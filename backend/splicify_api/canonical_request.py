"""
Canonical cloning-request schema.

The single data structure that every workflow reads, produced by
`request_normalizer.normalize_request` from the raw message + uploaded files +
per-intent parser output.

Design goals:
- One schema, many workflows. Gibson, Golden Gate, Restriction, SDM, Gateway
  handlers should all consume `CloningRequest` without re-parsing the raw
  message.
- Non-destructive to the existing intent_result dict. Normalization enriches;
  it never requires existing callers to change their input shape.
- No third-party deps (pydantic isn't installed on the VPS). Plain dataclasses
  + to_dict/from_dict for JSON logging and LLM round-trips.

Slot inventory (see README in module for discussion):
- parts:        the biological pieces the user wants assembled/modified
- vector:       the backbone / target plasmid when distinct from the parts
- enzymes:      user-specified restriction / type-IIs enzymes
- mutations:    for SDM — structured list, supports multi-site requests
- constraints:  negative/positive constraints ("avoid internal BsaI",
                "keep native Kozak", "preserve reading frame")
- assembly:     requested assembly method + topology, when the user stated it
- pcr_params / gibson_params / golden_gate_params / restriction_params /
  gateway_params / sdm_params: preserved pass-through of the per-intent LLM
  parse so handlers can transition gradually.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# Leaf types
# ---------------------------------------------------------------------------

@dataclass
class Part:
    """
    A piece of DNA referenced by the user.

    Exactly one of `sequence` or `name` is required; both are allowed. When only
    `name` is present, the handler is expected to resolve it via the pLannotate
    feature KB. When `sequence` is present, the handler may skip KB lookup.

    `role` captures the grammatical/functional role the user assigned:
      - "insert"    — the thing being cloned in
      - "backbone"  — the receiving vector (also see top-level `vector`)
      - "fragment"  — an unordered piece in a multi-fragment assembly
      - "promoter" / "cds" / "terminator" / "polyA" / "utr" / "tag" /
        "linker" / "selection_marker" / "origin" / "mcs" — finer typing
        mostly for Golden Gate / Gateway module resolution
      - "unknown"   — role couldn't be inferred

    `source` says where this Part came from, for audit/debug:
      - "prompt_sequence" — literal DNA in the user message
      - "prompt_name"     — a name in the user message, resolved via KB
      - "target_file"     — the uploaded target .gb
      - "inventory_file"  — one of the uploaded inventory .gb files
      - "kb_lookup"       — resolved from pLannotate KB by name
      - "derived"         — constructed by the backend (e.g. PCR amplicon)
    """
    name: Optional[str] = None
    sequence: Optional[str] = None
    role: str = "unknown"
    source: str = "prompt_name"
    source_file: Optional[str] = None
    # Ordering hint when the user specified a positional order (e.g.
    # "EF1a + eGFP + bGH polyA" → 0, 1, 2). None means order-independent.
    order: Optional[int] = None
    # Free-form bag for handler-specific extras (kb_sseqid, feature_type,
    # tm hint, etc.) so downstream code doesn't need new fields to experiment.
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Vector:
    """
    The backbone/target plasmid, when it has a distinct role from the parts
    list. Restriction / Gateway / SDM all need this; Gibson's inventory version
    does too.
    """
    name: Optional[str] = None
    sequence: Optional[str] = None
    source: str = "prompt_name"
    source_file: Optional[str] = None
    is_circular: bool = True
    # e.g. "pUC19", "pDONR221", "lentiCRISPR v2"
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Enzyme:
    """
    A restriction / type-IIs enzyme named in the request. Kept as its own type
    so handlers can carry buffer/heat-inactivation/star-activity metadata
    without stringly-typed gymnastics.
    """
    name: str
    # "type_ii" | "type_iis" | "nicking" | "homing" | "unknown"
    kind: str = "unknown"
    # User-expressed preference strength: "required", "preferred", "avoid"
    preference: str = "required"


@dataclass
class Mutation:
    """
    One site-directed mutagenesis request. A multi-site request becomes a list
    of these; each handler decides whether it can combine them in one primer
    pair (Q5 can, within a window) or needs sequential reactions.
    """
    # "substitution" | "deletion" | "insertion"
    mutation_type: str
    # "codon" | "feature" | "position" | "sequence" | "truncation"
    target_method: str
    target_feature_name: Optional[str] = None
    target_position_start: Optional[int] = None
    target_position_end: Optional[int] = None
    codon_position: Optional[int] = None
    codon_from: Optional[str] = None
    codon_to: Optional[str] = None
    # "N" | "C" | None — for N/C-terminal tag insertions and truncations
    terminus: Optional[str] = None
    # For direct-sequence mode or inserted payload sequence
    old_sequence: Optional[str] = None
    new_sequence: Optional[str] = None
    # For truncations: how many residues / bp to remove from the terminus
    truncation_length: Optional[int] = None
    description: Optional[str] = None


@dataclass
class Constraint:
    """
    A positive or negative constraint the user expressed that handlers should
    respect. Generic on purpose — the `kind` vocabulary can grow without
    schema changes.

    Known kinds:
      - "avoid_internal_site"    target="BsaI"
      - "preserve_feature"       target="Kozak" / "CMV promoter"
      - "preserve_reading_frame" target=feature_name or None
      - "junction_at"            target=sequence motif or feature_name
      - "require_domestication"  target=enzyme_name
      - "use_standard"           target="moclo" | "loop" | "mobius"
      - "max_fragments"          target=<int as str>
      - "selection"              target="ampR" | "kanR" | ...
    """
    kind: str
    target: Optional[str] = None
    # "positive" (must) | "negative" (must not) | "preferred" | "dispreferred"
    polarity: str = "positive"
    notes: Optional[str] = None


@dataclass
class AssemblyPlan:
    """
    User-stated assembly intent. May be None when the router infers it from
    the intent label alone.
    """
    # "gibson" | "golden_gate" | "restriction" | "gateway" | "sdm" | "pcr" |
    # "inv_gib" | "plasmid_design" | "annotate"
    method: Optional[str] = None
    topology: Optional[str] = None  # "linear" | "circular"
    # Ordered list of sub-methods when the request is composite, e.g.
    # ["golden_gate", "sdm"] → assemble, then mutate. Kept for future use.
    steps: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Top-level request
# ---------------------------------------------------------------------------

@dataclass
class CloningRequest:
    """
    Canonical representation of a user's cloning request.

    Populated by `request_normalizer.normalize_request`. Handlers read from
    this object (progressively — handlers that haven't migrated still read
    `intent_result`).
    """
    intent: str
    parts: List[Part] = field(default_factory=list)
    vector: Optional[Vector] = None
    enzymes: List[Enzyme] = field(default_factory=list)
    mutations: List[Mutation] = field(default_factory=list)
    constraints: List[Constraint] = field(default_factory=list)
    assembly: AssemblyPlan = field(default_factory=AssemblyPlan)

    # Pass-through of the per-intent LLM parse. Handlers in transition read
    # from here; once migrated they read from `parts`/`vector`/etc. directly.
    gibson_params: Optional[Dict[str, Any]] = None
    pcr_params: Optional[Dict[str, Any]] = None
    sdm_params: Optional[Dict[str, Any]] = None
    sgrna_params: Optional[Dict[str, Any]] = None
    golden_gate_params: Optional[Dict[str, Any]] = None
    restriction_params: Optional[Dict[str, Any]] = None
    gateway_params: Optional[Dict[str, Any]] = None

    # Provenance / diagnostics
    confidence: float = 1.0
    normalizer_notes: List[str] = field(default_factory=list)
    needs_clarification: bool = False
    clarification_question: Optional[str] = None
    # Verbatim original message, redacted DNA form, for audit/LLM round-trip.
    raw_message: Optional[str] = None
    redacted_message: Optional[str] = None

    # ------------------------------------------------------------------
    # Convenience lookups used by handlers
    # ------------------------------------------------------------------

    def parts_by_role(self, role: str) -> List[Part]:
        return [p for p in self.parts if p.role == role]

    def ordered_fragments(self) -> List[Part]:
        frags = [p for p in self.parts if p.role in ("fragment", "insert", "promoter", "cds", "terminator", "polyA", "utr", "tag", "linker")]
        frags.sort(key=lambda p: (p.order if p.order is not None else 10_000))
        return frags

    def constraint(self, kind: str) -> Optional[Constraint]:
        for c in self.constraints:
            if c.kind == kind:
                return c
        return None

    def has_constraint(self, kind: str, target: Optional[str] = None) -> bool:
        for c in self.constraints:
            if c.kind != kind:
                continue
            if target is None or c.target == target:
                return True
        return False

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CloningRequest":
        parts = [Part(**p) for p in data.get("parts", [])]
        vector = Vector(**data["vector"]) if data.get("vector") else None
        enzymes = [Enzyme(**e) for e in data.get("enzymes", [])]
        mutations = [Mutation(**m) for m in data.get("mutations", [])]
        constraints = [Constraint(**c) for c in data.get("constraints", [])]
        assembly = AssemblyPlan(**data.get("assembly", {})) if data.get("assembly") else AssemblyPlan()
        return cls(
            intent=data.get("intent", "unknown"),
            parts=parts,
            vector=vector,
            enzymes=enzymes,
            mutations=mutations,
            constraints=constraints,
            assembly=assembly,
            gibson_params=data.get("gibson_params"),
            pcr_params=data.get("pcr_params"),
            sdm_params=data.get("sdm_params"),
            sgrna_params=data.get("sgrna_params"),
            golden_gate_params=data.get("golden_gate_params"),
            restriction_params=data.get("restriction_params"),
            gateway_params=data.get("gateway_params"),
            confidence=data.get("confidence", 1.0),
            normalizer_notes=data.get("normalizer_notes", []) or [],
            needs_clarification=data.get("needs_clarification", False),
            clarification_question=data.get("clarification_question"),
            raw_message=data.get("raw_message"),
            redacted_message=data.get("redacted_message"),
        )
