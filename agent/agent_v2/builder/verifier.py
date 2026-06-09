"""Verifier: takes a VirtualConstruct + IntentSpec and emits a
structured verdict with diagnostics the builder can act on.

The verifier runs the materialized construct sequence through
annotate_llm_cached (same pipeline used everywhere else) so the
checks are run against the canonical annotation output the rest of
the system already trusts."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

from agent_v2.builder.intent_spec import IntentSpec
from agent_v2.builder.virtual_construct import VirtualConstruct, assess_methods


# Diagnostic codes the builder knows how to act on.
DIAG_MODULE_MISSING       = "module_missing"
DIAG_MODULE_PARTIAL       = "module_partial"
DIAG_INTERACTION_MISSING  = "interaction_missing"
DIAG_ORIENTATION_WRONG    = "orientation_wrong"
DIAG_ORDER_WRONG          = "order_wrong"
DIAG_JUNCTION_INCOMPAT    = "junction_incompatible"
DIAG_ROLE_MISSING         = "role_missing"
DIAG_FORBIDDEN_FEATURE    = "forbidden_feature_present"
DIAG_NO_METHOD            = "no_assembly_method_feasible"


@dataclass
class Diagnostic:
    code: str
    detail: str
    slot_idx: Optional[int] = None
    suggested_action: Optional[str] = None     # human-readable
    structured_action: Optional[dict[str, Any]] = None   # builder uses this

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "detail": self.detail,
            "slot_idx": self.slot_idx,
            "suggested_action": self.suggested_action,
            "structured_action": self.structured_action,
        }


@dataclass
class VerifierResult:
    passed: bool
    diagnostics: list[Diagnostic] = field(default_factory=list)
    annotation_summary: dict[str, Any] = field(default_factory=dict)
    method_pick: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "diagnostics": [d.to_dict() for d in self.diagnostics],
            "annotation_summary": self.annotation_summary,
            "method_pick": self.method_pick,
        }


# --- expected-shape mapping: function → required modules/interactions/roles ---
FUNCTION_TO_CHECKS = {
    "expression": {
        "modules": ["mammalian_pol2_expression_cassette"],
        "interactions": [],
        "roles": ["promoter", "cds", "polya"],
    },
    "sgrna_cloning": {
        "modules": ["guide_expression_cassette"],
        "interactions": [],
        "roles": ["promoter", "scaffold", "stuffer"],
    },
    "lentiviral_vector": {
        "modules": ["lentiviral_payload"],
        "interactions": ["lentiviral_three_module"],
        "roles": ["promoter", "cds", "polya", "ltr"],
    },
    "crispr_knockout": {
        "modules": ["guide_expression_cassette", "mammalian_pol2_expression_cassette"],
        "interactions": [],
        "roles": ["promoter", "cds", "scaffold"],
    },
    "cloning_vector": {
        "modules": ["bacterial_selection_cassette"],
        "interactions": [],
        "roles": ["selection_marker", "origin", "polylinker"],
    },
}


async def verify(
    construct: VirtualConstruct,
    intent: IntentSpec,
    *,
    annotate_fn=None,
) -> VerifierResult:
    """Run every deterministic primitive against the materialized
    construct + the rule-based detector's annotation output."""
    if annotate_fn is None:
        from splicify_api.annotation_cache import annotate_llm_cached
        annotate_fn = annotate_llm_cached

    diagnostics: list[Diagnostic] = []
    sequence = construct.materialize()

    if not sequence:
        diagnostics.append(Diagnostic(
            code="empty_construct",
            detail="VirtualConstruct has no parts; nothing to verify.",
        ))
        return VerifierResult(passed=False, diagnostics=diagnostics)

    try:
        ann = await annotate_fn(sequence, circular=intent.topology == "circular")
    except Exception as e:
        diagnostics.append(Diagnostic(
            code="annotation_failed",
            detail=f"{type(e).__name__}: {e}",
        ))
        return VerifierResult(passed=False, diagnostics=diagnostics)

    modules = ann.get("modules") or []
    interactions = ann.get("interactions") or []
    module_types = {m.get("module_type") for m in modules}
    interaction_types = {ix.get("type") or ix.get("interaction_type") for ix in interactions}

    # 1. Required modules.
    expected = FUNCTION_TO_CHECKS.get(intent.function, {})
    required_modules = list(intent.required_modules or expected.get("modules", []))
    for mt in required_modules:
        if mt in module_types:
            continue
        # partial check: is any module containing the same role chain present?
        partial = any(
            mt.split("_")[0] in (m.get("module_type") or "")
            for m in modules
        )
        diagnostics.append(Diagnostic(
            code=DIAG_MODULE_PARTIAL if partial else DIAG_MODULE_MISSING,
            detail=f"Required module {mt!r} not detected. Modules present: {sorted(module_types)}",
            structured_action={
                "action": "request_more_parts" if not partial else "reorder_or_reorient",
                "missing_module": mt,
            },
        ))

    # 2. Required interactions.
    required_interactions = list(intent.required_interactions or expected.get("interactions", []))
    for it in required_interactions:
        if it in interaction_types:
            continue
        diagnostics.append(Diagnostic(
            code=DIAG_INTERACTION_MISSING,
            detail=f"Required interaction {it!r} not detected. Present: {sorted(t for t in interaction_types if t)}",
            structured_action={"action": "reorder_or_reorient", "missing_interaction": it},
        ))

    # 3. Required roles in the slot composition.
    role_counts: dict[str, int] = {}
    for s in construct.slots:
        role_counts[s.part.role] = role_counts.get(s.part.role, 0) + 1
    required_roles = list(intent.required_roles or expected.get("roles", []))
    for role in required_roles:
        if role_counts.get(role, 0) == 0:
            diagnostics.append(Diagnostic(
                code=DIAG_ROLE_MISSING,
                detail=f"No slot carries role {role!r}. Roles present: {role_counts}",
                structured_action={"action": "request_more_parts", "missing_role": role},
            ))

    # 4. Orientation check — promoter must precede CDS on the same
    # strand; polyA must follow CDS on the same strand. We walk
    # adjacent role triples.
    for i, slot in enumerate(construct.slots):
        if slot.part.role != "cds":
            continue
        before = construct.slots[i - 1] if i > 0 else None
        after = construct.slots[i + 1] if i + 1 < len(construct.slots) else None
        # When circular, allow wrap-around for the "before" / "after" of edges.
        if construct.topology == "circular":
            before = before or construct.slots[-1]
            after = after or construct.slots[0]
        if before and before.part.role == "promoter" and before.orientation != slot.orientation:
            diagnostics.append(Diagnostic(
                code=DIAG_ORIENTATION_WRONG, slot_idx=i - 1 if i > 0 else len(construct.slots) - 1,
                detail=f"Promoter at slot {i-1} ({before.part.name}) is on opposite strand "
                       f"from its CDS at slot {i} ({slot.part.name}).",
                structured_action={"action": "swap_orientation", "slot_idx": i - 1},
            ))
        if after and after.part.role == "polya" and after.orientation != slot.orientation:
            diagnostics.append(Diagnostic(
                code=DIAG_ORIENTATION_WRONG, slot_idx=i + 1 if i + 1 < len(construct.slots) else 0,
                detail=f"polyA at slot {i+1} ({after.part.name}) is on opposite strand "
                       f"from its CDS at slot {i} ({slot.part.name}).",
                structured_action={"action": "swap_orientation", "slot_idx": i + 1},
            ))

    # 5. Order check — promoter → CDS → polyA pattern around every CDS.
    for i, slot in enumerate(construct.slots):
        if slot.part.role != "cds":
            continue
        before = construct.slots[(i - 1) % len(construct.slots)] if construct.slots else None
        after = construct.slots[(i + 1) % len(construct.slots)] if construct.slots else None
        if before and after:
            if before.part.role not in ("promoter", "kozak", "tag") and \
               after.part.role not in ("polya", "tag", "linker", "ltr"):
                diagnostics.append(Diagnostic(
                    code=DIAG_ORDER_WRONG, slot_idx=i,
                    detail=f"CDS at slot {i} ({slot.part.name}) is not flanked by "
                           f"a 5' promoter/kozak/tag or a 3' polyA/tag/linker/ltr. "
                           f"Surrounding: {before.part.role!r} / {after.part.role!r}.",
                    structured_action={"action": "reorder", "around_slot": i},
                ))

    # 6. Forbidden features.
    for fname in intent.forbidden_features:
        for s in construct.slots:
            if fname.lower() in (s.part.name or "").lower():
                diagnostics.append(Diagnostic(
                    code=DIAG_FORBIDDEN_FEATURE,
                    detail=f"Forbidden feature {fname!r} present at slot containing {s.part.name!r}.",
                    structured_action={"action": "remove_part", "part_id": s.part.part_id},
                ))

    # 7. Assembly-method feasibility — at least one method must be
    # feasible given the junction profile.
    method = assess_methods(construct)
    if not method.feasible:
        diagnostics.append(Diagnostic(
            code=DIAG_NO_METHOD,
            detail=(
                "No assembly method works for this part set without adding new sequence. "
                f"Rejected: {method.rejected}"
            ),
            structured_action={"action": "request_more_parts",
                                "hint": "need parts with shared junction homology, "
                                        "Type IIs sites, or restriction sites"},
        ))

    summary = {
        "annotations_returned": len(ann.get("annotations") or []),
        "modules_present": sorted(t for t in module_types if t),
        "interactions_present": sorted(t for t in interaction_types if t),
        "construct_bp": len(sequence),
        "n_slots": len(construct.slots),
    }
    return VerifierResult(
        passed=not diagnostics,
        diagnostics=diagnostics,
        annotation_summary=summary,
        method_pick=method.pick,
    )
