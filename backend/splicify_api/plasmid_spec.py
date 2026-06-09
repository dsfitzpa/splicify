"""
PlasmidSpec — abstract representation of the plasmid the user wants.

Built from the intent + KB-resolved features + uploaded sequences. Carried
through the unified predesign pipeline so every stage (PartResolver,
TargetPlasmidBuilder, CloningRouter, target_from_inventory_router) compares
candidate targets against the same spec.

This is the missing layer between "user prompt" and "assembled target": it is
the answer to "what do we believe the user wants on the final plasmid?" and is
what the (future) LLM orchestrator will reason over.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# Module / role inference from the prompt
# ---------------------------------------------------------------------------
_HOST_KEYWORDS = {
    "mammalian": ("mammalian", "human cell", "hek293", "cho", "u2os"),
    "lentiviral": ("lentivir", "lentiv", "lvv", "ltr", "wpre", "rre"),
    "aav": ("aav", "adeno-associated"),
    "bacterial": ("bacterial", "e. coli", "e.coli", "ecoli", "iptg", "t7"),
    "yeast": ("yeast", "saccharomyces", "pichia"),
    "plant": ("plant", "agrobacterium", "t-dna", "arabidopsis"),
}

_MODULE_ROLE_KEYWORDS = {
    "promoter":         ("promoter", "promoters"),
    "polyA_signal":     ("polya", "poly-a", "poly a", "polyadenylation"),
    "terminator":       ("terminator",),
    "selection_marker": ("selection", "puromycin", "puror", "blasticidin", "neor",
                          "hygror", "ampicillin", "kanamycin"),
    "ori":              ("origin of replication", "ori "),
    "nls":              ("nls", "nuclear localization"),
    "tag":              ("flag tag", "ha tag", "myc tag", "v5 tag", "his tag",
                          "his-tag", "6xhis", "strep tag"),
    "reporter":         ("egfp", "gfp", "mcherry", "mtagrfp", "mscarlet",
                          "luciferase", "rluc", "fluc"),
    "nuclease":         ("cas9", "dcas9", "cas12", "cas13", "talen", "zinc finger"),
    "guide_cassette":   ("sgrna", "grna", "guide rna"),
    "wpre":             ("wpre",),
    "lentiviral":       ("ltr", "5' ltr", "3' ltr", "rre", "psi", "cppt"),
}


@dataclass
class SpecModule:
    """A module the user wants on the final plasmid.

    `origin` is "described" (named in the prompt), "uploaded" (present in an
    input file), or "inferred" (added by the spec builder for biological
    completeness, e.g. ori + AmpR for mammalian expression).
    """
    name: str
    role: Optional[str] = None
    origin: str = "described"  # described | uploaded | inferred
    sequence: Optional[str] = None
    canonical_id: Optional[str] = None
    feature_id: Optional[str] = None
    length: Optional[int] = None


@dataclass
class PlasmidSpec:
    modules_required: List[SpecModule] = field(default_factory=list)
    modules_present:  List[SpecModule] = field(default_factory=list)
    assembly_hint:    Optional[str] = None  # gateway/gibson/golden_gate/...
    host_hint:        Optional[str] = None  # mammalian/lentiviral/bacterial/...
    topology:         str = "circular"
    description:      str = ""

    def required_role_set(self) -> set:
        return {m.role for m in self.modules_required if m.role}

    def to_dict(self) -> Dict[str, Any]:
        return {
            "modules_required": [m.__dict__ for m in self.modules_required],
            "modules_present":  [m.__dict__ for m in self.modules_present],
            "assembly_hint": self.assembly_hint,
            "host_hint": self.host_hint,
            "topology": self.topology,
            "description": self.description,
        }


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------
def build_plasmid_spec(
    *,
    message: str,
    intent_result: Dict[str, Any],
    uploaded_modules: Optional[List[Dict[str, Any]]] = None,
) -> PlasmidSpec:
    """Build a PlasmidSpec from intent + uploaded modules.

    `uploaded_modules` is the flat list of module annotations from any uploaded
    .gb files (already produced by the annotation pipeline). Each entry should
    have at least {name, role|module_type, sequence?}.
    """
    msg = (message or "").lower()
    intent = (intent_result or {}).get("intent")

    spec = PlasmidSpec(
        description=message or "",
        topology="linear" if "linear" in msg else "circular",
        assembly_hint=_assembly_hint_from_intent(intent),
        host_hint=_host_hint_from_text(msg),
    )

    # Modules the user explicitly named (KB-resolved if possible)
    kb = (intent_result or {}).get("kb_resolved") or {}
    for cand in kb.get("candidates") or []:
        spec.modules_required.append(SpecModule(
            name=cand.get("name") or "",
            role=cand.get("feature_type"),
            origin="described",
        ))
    # Replace with KB-identified entries when available (carry sequence/feature_id)
    by_name = {m.name.lower(): m for m in spec.modules_required}
    for ident in kb.get("identified") or []:
        nm = (ident.get("query") or ident.get("name") or "").strip()
        if not nm:
            continue
        m = by_name.get(nm.lower())
        if m is None:
            m = SpecModule(name=nm, origin="described")
            spec.modules_required.append(m)
            by_name[nm.lower()] = m
        m.role = ident.get("feature_type") or m.role
        m.sequence = ident.get("sequence") or m.sequence
        m.canonical_id = ident.get("name") or ident.get("feature_id") or m.canonical_id
        m.feature_id = ident.get("feature_id") or m.feature_id
        m.length = ident.get("length") or m.length

    # Modules already present in any uploaded sequence
    for mod in uploaded_modules or []:
        spec.modules_present.append(SpecModule(
            name=str(mod.get("name") or mod.get("module_type") or ""),
            role=str(mod.get("role") or mod.get("module_type") or "") or None,
            origin="uploaded",
            sequence=mod.get("sequence"),
            length=mod.get("length") or (
                len(mod["sequence"]) if mod.get("sequence") else None
            ),
        ))

    # Role keywords mentioned without a specific KB part
    for role, kws in _MODULE_ROLE_KEYWORDS.items():
        if any(kw in msg for kw in kws):
            already = any(
                (m.role or "").lower() == role
                for m in spec.modules_required + spec.modules_present
            )
            if not already:
                spec.modules_required.append(SpecModule(
                    name=role, role=role, origin="described",
                ))

    return spec


def _assembly_hint_from_intent(intent: Optional[str]) -> Optional[str]:
    return {
        "gibson_design": "gibson",
        "gateway_cloning": "gateway",
        "golden_gate_primer_design": "golden_gate",
        "sgrna_golden_gate": "sgrna_golden_gate",
        "restriction_cloning": "restriction",
        "sdm_design": "sdm",
        "plasmid_design": "describe",
        "annotate_gb": None,
        "unknown": None,
    }.get(intent or "")


def _host_hint_from_text(msg: str) -> Optional[str]:
    for host, kws in _HOST_KEYWORDS.items():
        if any(kw in msg for kw in kws):
            return host
    return None


# ---------------------------------------------------------------------------
# Spec ↔ assembled target diff
# ---------------------------------------------------------------------------
def diff_spec_against_target(
    spec: PlasmidSpec,
    target_modules: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    """Compare what the spec asks for to what is actually on the assembled target.

    Returns:
        {
          "satisfied":  [{name, role}],   # required & present in target
          "missing":    [{name, role}],   # required & absent in target
          "unexpected": [{name, role}],   # in target but not requested
        }
    """
    target_role_to_names: Dict[str, List[str]] = {}
    target_names: List[str] = []
    for m in target_modules or []:
        nm = str(m.get("name") or m.get("module_type") or "")
        rl = str(m.get("role") or m.get("module_type") or "") or None
        target_names.append(nm.lower())
        if rl:
            target_role_to_names.setdefault(rl.lower(), []).append(nm.lower())

    satisfied: List[Dict[str, Any]] = []
    missing:   List[Dict[str, Any]] = []
    seen_target_idxs: set = set()

    for req in spec.modules_required:
        nm = (req.name or "").lower()
        rl = (req.role or "").lower() if req.role else None
        ok = False
        if nm and nm in target_names:
            seen_target_idxs.add(target_names.index(nm))
            ok = True
        elif rl and target_role_to_names.get(rl):
            seen_target_idxs.add(target_names.index(target_role_to_names[rl][0]))
            ok = True
        (satisfied if ok else missing).append({"name": req.name, "role": req.role})

    unexpected: List[Dict[str, Any]] = []
    for i, m in enumerate(target_modules or []):
        if i in seen_target_idxs:
            continue
        unexpected.append({
            "name": m.get("name") or m.get("module_type"),
            "role": m.get("role") or m.get("module_type"),
        })

    return {"satisfied": satisfied, "missing": missing, "unexpected": unexpected}
