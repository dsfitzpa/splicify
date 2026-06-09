"""
WorkflowInput — the canonical input shape every cloning-workflow handler reads.

Goal: collapse the three-way mismatch where each handler currently re-derives
its inputs from a different bag (raw `message`, `seq_data['sequences']`,
`intent_result['<intent>']`, `canonical_request`, raw upload bytes). After
this scaffolding lands, every handler will accept a single `WorkflowInput` —
the same shape regardless of whether the request came down the predesign path
(`_execute_unified_predesign`) or the inventory-router path
(`target_from_inventory_router.route_from_uploads`).

Step 1 (this file): pure additive scaffolding. No handler change, no behavior
change. The dataclasses + projection builders are exercised by unit tests in
`backend/tests/test_workflow_input.py`.

Companion docs: PROJECT_SUMMARY.md §3, CLONING_WORKFLOWS.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, TYPE_CHECKING

from .canonical_request import AssemblyPlan, CloningRequest, Constraint

if TYPE_CHECKING:
    # Type-only imports to avoid runtime circulars. The actual classes are
    # imported lazily inside the builders that need them.
    from .predesign.part_resolver import ResolvedPart
    from .predesign.target_builder import TargetPlasmid
    from .target_from_inventory_router import PlasmidContext


# ---------------------------------------------------------------------------
# PlasmidView — unified annotated-plasmid wrapper
# ---------------------------------------------------------------------------

@dataclass
class PlasmidView:
    """
    Unified projection of every annotated plasmid the handler might see:
      - the assembled target (predesign-driven path)
      - the uploaded target (inventory-driven path)
      - one resolved part (predesign-driven path)
      - one inventory plasmid (inventory-driven path)
      - a KB-derived insert
      - a derived intermediate (PCR amplicon, gBlock, ...)

    Superset of router's `PlasmidContext` and predesign's `ResolvedPart` +
    annotation pair. Either upstream pipeline projects losslessly into this.

    Role vocabulary (free-form string; handler-specific values allowed):
      - "target"          assembled / uploaded plasmid being built or edited
      - "vector"          backbone distinct from target (restriction / sgRNA)
      - "insert"          part being cloned in
      - "donor"           Gateway donor (carries attL/attR or attP entry)
      - "destination"     Gateway destination vector
      - "template"        SDM source plasmid
      - "fragment"        unordered piece in a multi-fragment assembly
      - "promoter" / "cds" / "terminator" / "polyA" / "utr" / "tag" /
        "linker" / "selection_marker" / "origin" / "mcs"
                          finer typing for Golden Gate / Gateway resolution
      - "auxiliary"       inventory plasmid not selected by the chosen workflow
      - "unknown"         role couldn't be inferred

    Source vocabulary:
      - "uploaded"               from `file` or `inventory_files`
      - "kb_lookup"              resolved from pLannotate feature KB
      - "prompt_sequence"        literal DNA in user message
      - "predesign_assembled"    target built by TargetPlasmidBuilder
      - "derived"                constructed by the backend (PCR, synth)
    """

    name: str
    sequence: str
    role: str = "unknown"
    source: str = "uploaded"
    source_file: Optional[str] = None
    topology: str = "linear"
    length: int = 0
    gb_text: Optional[str] = None
    annotations: List[Dict[str, Any]] = field(default_factory=list)
    modules: List[Dict[str, Any]] = field(default_factory=list)
    interactions: List[Dict[str, Any]] = field(default_factory=list)
    cloning_features: Any = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.sequence and self.length == 0:
            self.length = len(self.sequence)

    @classmethod
    def from_plasmid_context(
        cls,
        ctx: "PlasmidContext",
        role: str = "unknown",
        source: str = "uploaded",
        source_file: Optional[str] = None,
    ) -> "PlasmidView":
        """Project a router-side PlasmidContext (target_from_inventory_router).

        PlasmidContext fields: name, sequence, gb_text, annotations, modules,
        cloning. Topology is not on PlasmidContext — the router minimal-GB
        wrapper assumes circular; we preserve that default.
        """
        return cls(
            name=ctx.name,
            sequence=ctx.sequence,
            role=role,
            source=source,
            source_file=source_file,
            topology="circular",
            length=len(ctx.sequence),
            gb_text=ctx.gb_text,
            annotations=list(ctx.annotations or []),
            modules=list(ctx.modules or []),
            interactions=[],
            cloning_features=ctx.cloning,
            metadata={},
        )

    @classmethod
    def from_resolved_part(
        cls,
        part: "ResolvedPart",
        annotation: Optional[Dict[str, Any]] = None,
        role: Optional[str] = None,
    ) -> "PlasmidView":
        """Project a predesign ResolvedPart (+ optional annotation pair).

        `annotation` is one entry from `predesign_context['part_annotations']`,
        i.e. {name, length, modules, interactions}. When None, modules and
        interactions are left empty.
        """
        return cls(
            name=part.name,
            sequence=part.sequence,
            role=role or part.role or "unknown",
            source=_source_from_resolved_origin(part.origin),
            source_file=part.source_detail,
            topology="linear",
            length=part.length,
            gb_text=None,
            annotations=list(part.features or []),
            modules=list((annotation or {}).get("modules") or []),
            interactions=list((annotation or {}).get("interactions") or []),
            cloning_features=None,
            metadata={
                "canonical_id": part.canonical_id,
                "description": part.description,
                "confidence": part.confidence,
            },
        )

    @classmethod
    def from_target_plasmid(
        cls,
        target: "TargetPlasmid",
        target_modules: Optional[List[Dict[str, Any]]] = None,
        target_interactions: Optional[List[Dict[str, Any]]] = None,
        name: str = "assembled_target",
    ) -> "PlasmidView":
        """Project a predesign TargetPlasmid into a PlasmidView with role='target'."""
        return cls(
            name=name,
            sequence=target.sequence,
            role="target",
            source="predesign_assembled",
            source_file=None,
            topology=target.topology or "circular",
            length=target.length,
            gb_text=None,
            annotations=list(target.features or []),
            modules=list(target_modules or []),
            interactions=list(target_interactions or []),
            cloning_features=None,
            metadata=dict(target.metadata or {}),
        )


def _source_from_resolved_origin(origin: str) -> str:
    """Map predesign ResolvedPart.origin → PlasmidView.source vocabulary."""
    if origin == "knowledge_base":
        return "kb_lookup"
    if origin == "user_file":
        return "uploaded"
    if origin == "user_provided":
        return "prompt_sequence"
    return "derived"


# ---------------------------------------------------------------------------
# WorkflowInput — top-level canonical input
# ---------------------------------------------------------------------------

@dataclass
class WorkflowInput:
    """
    Canonical input for every per-intent cloning workflow handler.

    Built once per request by `build_workflow_input` and consumed in place of
    raw `message` / `seq_data` / `intent_result['<intent>']` / `canonical_request`
    / per-handler upload re-parsing.
    """

    intent: str
    workflow: str
    session_id: str

    target: Optional[PlasmidView] = None
    parts: List[PlasmidView] = field(default_factory=list)
    vector: Optional[PlasmidView] = None

    plasmid_spec: Optional[Dict[str, Any]] = None
    spec_diff: Optional[Dict[str, Any]] = None

    constraints: List[Constraint] = field(default_factory=list)
    assembly: AssemblyPlan = field(default_factory=AssemblyPlan)

    workflow_args: Dict[str, Any] = field(default_factory=dict)

    # Provenance bag for audit lines + downstream LLM explanation.
    # Known keys: "routing_audit", "predesign_context", "normalizer_notes".
    provenance: Dict[str, Any] = field(default_factory=dict)

    def parts_by_role(self, role: str) -> List[PlasmidView]:
        return [p for p in self.parts if p.role == role]

    def part_by_name(self, name: str) -> Optional[PlasmidView]:
        for p in self.parts:
            if p.name == name:
                return p
        return None

    def has_predesign(self) -> bool:
        return self.plasmid_spec is not None

    def has_routing_audit(self) -> bool:
        return bool(self.provenance.get("routing_audit"))


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def build_from_predesign(
    *,
    canonical_request: Optional[CloningRequest],
    intent: str,
    session_id: str,
    resolved_parts: List["ResolvedPart"],
    target_plasmid: Optional["TargetPlasmid"],
    part_annotations: Optional[List[Dict[str, Any]]] = None,
    target_modules: Optional[List[Dict[str, Any]]] = None,
    target_interactions: Optional[List[Dict[str, Any]]] = None,
    plasmid_spec: Optional[Dict[str, Any]] = None,
    spec_diff: Optional[Dict[str, Any]] = None,
    predesign_context: Optional[Dict[str, Any]] = None,
) -> WorkflowInput:
    """Project predesign internals + canonical_request → WorkflowInput.

    Caller (chat.py / _execute_unified_predesign) is responsible for surfacing
    `resolved_parts` and `target_plasmid` from the pipeline; today they live
    only inside `_execute_unified_predesign` (chat.py:723) — promoting them
    is step 2 of the migration.
    """
    ann_by_name: Dict[str, Dict[str, Any]] = {}
    for ann in part_annotations or []:
        n = ann.get("name")
        if n:
            ann_by_name[n] = ann

    parts: List[PlasmidView] = []
    for rp in resolved_parts:
        parts.append(PlasmidView.from_resolved_part(rp, ann_by_name.get(rp.name)))

    target_view: Optional[PlasmidView] = None
    if target_plasmid is not None:
        target_view = PlasmidView.from_target_plasmid(
            target_plasmid,
            target_modules=target_modules,
            target_interactions=target_interactions,
        )

    workflow_args = _merge_canonical_params(canonical_request, intent)

    constraints: List[Constraint] = list((canonical_request.constraints if canonical_request else []) or [])
    assembly = (canonical_request.assembly if canonical_request else None) or AssemblyPlan()

    return WorkflowInput(
        intent=intent,
        workflow=intent,
        session_id=session_id,
        target=target_view,
        parts=parts,
        vector=None,
        plasmid_spec=plasmid_spec,
        spec_diff=spec_diff,
        constraints=constraints,
        assembly=assembly,
        workflow_args=workflow_args,
        provenance={
            "predesign_context": predesign_context,
            "normalizer_notes": list((canonical_request.normalizer_notes if canonical_request else []) or []),
        },
    )


def build_from_router(
    *,
    canonical_request: Optional[CloningRequest],
    session_id: str,
    routing_audit_payload: Dict[str, Any],
    target_ctx: "PlasmidContext",
    inventory_ctxs: List["PlasmidContext"],
    target_source_file: Optional[str] = None,
    inventory_source_files: Optional[List[Optional[str]]] = None,
) -> WorkflowInput:
    """Project router output (chosen + contexts) → WorkflowInput.

    The router's `chosen_handler_args` carry cross-references to inventory
    plasmids by name (donor_name, cargo_source, vector_source, template_source,
    fragments[*].source). We use those to assign roles when projecting the
    inventory list into `parts`.

    Multisite Gateway will eventually carry an explicit per-donor role on the
    router payload (separate router change). Until that lands, `donors[]`
    are tagged as role="donor" without per-position differentiation.
    """
    chosen_intent = routing_audit_payload.get("chosen_intent") or ""
    chosen_args: Dict[str, Any] = dict(routing_audit_payload.get("chosen_handler_args") or {})
    gateway_variant = chosen_args.get("gateway_variant")
    workflow = chosen_intent + ((":" + gateway_variant) if gateway_variant else "")

    role_by_name = _role_assignments_from_handler_args(chosen_intent, chosen_args)

    inventory_source_files = inventory_source_files or [None] * len(inventory_ctxs)

    parts: List[PlasmidView] = []
    vector_view: Optional[PlasmidView] = None
    for ctx, src_file in zip(inventory_ctxs, inventory_source_files):
        role = role_by_name.get(ctx.name, "auxiliary")
        view = PlasmidView.from_plasmid_context(
            ctx, role=role, source="uploaded", source_file=src_file,
        )
        if role == "vector":
            vector_view = view
        parts.append(view)

    target_view = PlasmidView.from_plasmid_context(
        target_ctx, role="target", source="uploaded", source_file=target_source_file,
    )

    merged_args = _merge_canonical_params(canonical_request, chosen_intent)
    merged_args.update(chosen_args)

    constraints: List[Constraint] = list((canonical_request.constraints if canonical_request else []) or [])
    assembly = (canonical_request.assembly if canonical_request else None) or AssemblyPlan()

    return WorkflowInput(
        intent=chosen_intent,
        workflow=workflow,
        session_id=session_id,
        target=target_view,
        parts=parts,
        vector=vector_view,
        plasmid_spec=None,
        spec_diff=None,
        constraints=constraints,
        assembly=assembly,
        workflow_args=merged_args,
        provenance={
            "routing_audit": routing_audit_payload,
            "normalizer_notes": list((canonical_request.normalizer_notes if canonical_request else []) or []),
        },
    )


def build_from_single_upload(
    *,
    canonical_request: Optional[CloningRequest],
    intent: str,
    session_id: str,
    target_view: PlasmidView,
) -> WorkflowInput:
    """Minimal builder for the single-upload case (annotate_gb, plasmid_design
    text-only, sdm with one file). The caller is responsible for running the
    annotation pipeline on the upload and constructing the PlasmidView."""
    return WorkflowInput(
        intent=intent,
        workflow=intent,
        session_id=session_id,
        target=target_view,
        parts=[],
        vector=None,
        plasmid_spec=None,
        spec_diff=None,
        constraints=list((canonical_request.constraints if canonical_request else []) or []),
        assembly=(canonical_request.assembly if canonical_request else None) or AssemblyPlan(),
        workflow_args=_merge_canonical_params(canonical_request, intent),
        provenance={
            "normalizer_notes": list((canonical_request.normalizer_notes if canonical_request else []) or []),
        },
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_INTENT_TO_PARAM_FIELD: Dict[str, str] = {
    "gibson_design": "gibson_params",
    "golden_gate_primer_design": "golden_gate_params",
    "restriction_cloning": "restriction_params",
    "gateway_cloning": "gateway_params",
    "sdm_design": "sdm_params",
    "sgrna_golden_gate": "sgrna_params",
}


def _merge_canonical_params(
    canonical_request: Optional[CloningRequest], intent: str
) -> Dict[str, Any]:
    if canonical_request is None:
        return {}
    field_name = _INTENT_TO_PARAM_FIELD.get(intent)
    if not field_name:
        return {}
    raw = getattr(canonical_request, field_name, None)
    return dict(raw) if isinstance(raw, dict) else {}


def _role_assignments_from_handler_args(
    chosen_intent: str, args: Dict[str, Any]
) -> Dict[str, str]:
    """Inspect the router's `chosen_handler_args` and emit a name → role map
    that `build_from_router` applies when projecting inventory into `parts`.
    """
    out: Dict[str, str] = {}

    if chosen_intent == "gateway_cloning":
        if args.get("gateway_variant") == "multisite":
            for d in args.get("donors") or []:
                if d:
                    out[d] = "donor"
            dest = args.get("destination")
            if dest:
                out[dest] = "destination"
        else:
            donor = args.get("donor_name")
            if donor:
                out[donor] = "donor"
            cargo = args.get("cargo_source")
            if cargo:
                out[cargo] = "insert"

    elif chosen_intent == "sgrna_golden_gate":
        v = args.get("vector_source")
        if v:
            out[v] = "vector"

    elif chosen_intent == "sdm_design":
        t = args.get("template_source")
        if t:
            out[t] = "template"

    elif chosen_intent in ("restriction_cloning", "golden_gate_primer_design"):
        for f in args.get("fragments") or []:
            if isinstance(f, dict):
                src = f.get("source")
                if src:
                    out[src] = "fragment"
            elif isinstance(f, str):
                out[f] = "fragment"

    elif chosen_intent == "gibson_design":
        for f in args.get("fragments") or []:
            if isinstance(f, dict):
                src = f.get("source")
                if src:
                    out[src] = "fragment"

    elif chosen_intent == "plasmid_design":
        for a in args.get("inventory_anchors") or []:
            if isinstance(a, dict):
                src = a.get("source")
                if src:
                    out[src] = "fragment"

    return out


__all__ = [
    "PlasmidView",
    "WorkflowInput",
    "build_from_predesign",
    "build_from_router",
    "build_from_single_upload",
]
