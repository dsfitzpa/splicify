"""Builder loop: take a PartSet + IntentSpec, propose an ordered
VirtualConstruct, hand it to the verifier, apply mechanical
modifications based on the verdict, retry up to N times. No new
sequence is ever added — only reorder / reorient / swap parts."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

from agent_v2.builder.intent_spec import IntentSpec
from agent_v2.builder.part_set import PartSet, Part
from agent_v2.builder.verifier import (
    VerifierResult, Diagnostic, verify,
    DIAG_ORIENTATION_WRONG, DIAG_ORDER_WRONG,
    DIAG_MODULE_MISSING, DIAG_MODULE_PARTIAL,
    DIAG_ROLE_MISSING, DIAG_FORBIDDEN_FEATURE,
    DIAG_INTERACTION_MISSING, DIAG_JUNCTION_INCOMPAT, DIAG_NO_METHOD,
)
from agent_v2.builder.virtual_construct import (
    Slot, VirtualConstruct, assess_methods, MethodAssessment,
)


# Canonical role order for an expression cassette. The builder tries
# this layout first, then permutations if verification fails.
DEFAULT_ORDER = [
    "ltr",                  # 5' LTR (lentiviral)
    "lentiviral_cis",       # ψ, RRE, cPPT
    "enhancer",
    "promoter",
    "kozak",
    "tag",                  # N-terminal tag
    "cds",
    "polya",
    "wpre",
    "scaffold",             # gRNA scaffold (POL3)
    "stuffer",
    "selection_marker",
    "att_site",
    "polylinker",
    "origin",
]


@dataclass
class BuildJournalEntry:
    iteration: int
    action: str
    detail: str
    construct_snapshot: dict[str, Any] = field(default_factory=dict)
    verifier_result: dict[str, Any] = field(default_factory=dict)


@dataclass
class BuildResult:
    success: bool
    final_construct: Optional[VirtualConstruct] = None
    final_verification: Optional[VerifierResult] = None
    method_pick: Optional[str] = None
    method_assessment: Optional[MethodAssessment] = None
    journal: list[BuildJournalEntry] = field(default_factory=list)
    unresolved_diagnostics: list[Diagnostic] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "success": self.success,
            "method_pick": self.method_pick,
            "method_assessment": (self.method_assessment.to_dict()
                                   if self.method_assessment else None),
            "final_construct": (self.final_construct.to_dict()
                                 if self.final_construct else None),
            "final_verification": (self.final_verification.to_dict()
                                    if self.final_verification else None),
            "journal": [{"iteration": e.iteration, "action": e.action,
                          "detail": e.detail} for e in self.journal],
            "unresolved_diagnostics": [d.to_dict() for d in self.unresolved_diagnostics],
        }


def _initial_layout(part_set: PartSet) -> VirtualConstruct:
    """Pick a starting order using DEFAULT_ORDER for known roles, then
    everything else by source-plasmid order."""
    by_role: dict[str, list[Part]] = {}
    for p in part_set.parts:
        by_role.setdefault(p.role, []).append(p)
    ordered: list[Slot] = []
    for role in DEFAULT_ORDER:
        for p in by_role.pop(role, []):
            ordered.append(Slot(part=p, orientation=p.source_strand))
    # Anything with an unknown role gets appended at the end in
    # source order.
    leftover: list[Part] = []
    for role_parts in by_role.values():
        leftover.extend(role_parts)
    leftover.sort(key=lambda p: (p.source_plasmid_id, p.source_start))
    for p in leftover:
        ordered.append(Slot(part=p, orientation=p.source_strand))
    return VirtualConstruct(slots=ordered, topology="circular")


def _apply_action(construct: VirtualConstruct, diag: Diagnostic,
                  part_set: PartSet) -> tuple[bool, str]:
    """Apply ONE mechanical modification to the construct based on a
    diagnostic. Returns (applied, detail). No new sequence added."""
    a = diag.structured_action or {}
    action = a.get("action")
    if not action:
        return False, "no structured_action on diagnostic"

    if action == "swap_orientation":
        idx = a.get("slot_idx")
        if idx is None or idx >= len(construct.slots):
            return False, f"swap_orientation: bad slot_idx {idx}"
        construct.slots[idx].orientation *= -1
        return True, f"flipped slot {idx} ({construct.slots[idx].part.name})"

    if action == "reorder":
        # Promote a CDS into a window flanked by a promoter (left) and
        # a polyA (right). If we can find a promoter not adjacent to a
        # CDS, move it next to the closest CDS.
        cds_idx = a.get("around_slot")
        if cds_idx is None or cds_idx >= len(construct.slots):
            return False, f"reorder: bad slot_idx {cds_idx}"
        slots = construct.slots
        # Find a promoter slot to move next to cds_idx
        prom_idx = next((i for i, s in enumerate(slots) if s.part.role == "promoter"
                          and abs(i - cds_idx) > 1), None)
        polya_idx = next((i for i, s in enumerate(slots) if s.part.role == "polya"
                          and abs(i - cds_idx) > 1), None)
        moved = []
        if prom_idx is not None:
            prom_slot = slots.pop(prom_idx)
            new_pos = cds_idx if prom_idx > cds_idx else cds_idx - 1
            slots.insert(max(0, new_pos), prom_slot)
            moved.append(f"promoter→{new_pos}")
        if polya_idx is not None and polya_idx < len(slots):
            polya_slot = slots.pop(polya_idx)
            cds_in = next((i for i, s in enumerate(slots) if s.part.role == "cds"), cds_idx)
            slots.insert(min(len(slots), cds_in + 1), polya_slot)
            moved.append(f"polya→{cds_in + 1}")
        if not moved:
            return False, "reorder: no movable promoter/polya found"
        return True, "; ".join(moved)

    if action == "reorder_or_reorient":
        # The simpler tactic: try flipping each non-matching slot's
        # orientation toward the construct's majority strand.
        if not construct.slots:
            return False, "no slots to reorient"
        majority = max(set(s.orientation for s in construct.slots),
                        key=lambda o: sum(1 for s in construct.slots if s.orientation == o))
        flipped = 0
        for s in construct.slots:
            if s.orientation != majority and s.part.role in {"promoter", "polya", "cds", "ltr"}:
                s.orientation = majority
                flipped += 1
        if flipped:
            return True, f"flipped {flipped} slots to strand {majority}"
        return False, "no flipping possible"

    if action == "remove_part":
        pid = a.get("part_id")
        if pid is None:
            return False, "remove_part: missing part_id"
        before = len(construct.slots)
        construct.slots = [s for s in construct.slots if s.part.part_id != pid]
        if len(construct.slots) == before:
            return False, f"remove_part: part_id {pid} not found"
        return True, f"removed part {pid}"

    if action == "request_more_parts":
        # Try to pull a candidate from part_set.candidates that fills
        # the missing role / module.
        wanted_role = a.get("missing_role")
        wanted_module = a.get("missing_module")
        target = None
        if wanted_role:
            target = next((c for c in part_set.candidates if c.role == wanted_role), None)
        if not target and wanted_module:
            target = next((c for c in part_set.candidates
                           if (c.name or "").lower().find(wanted_module.lower()) >= 0), None)
        if not target:
            return False, (
                f"request_more_parts: no candidate matches "
                f"role={wanted_role!r} module={wanted_module!r}"
            )
        construct.slots.append(Slot(part=target, orientation=target.source_strand))
        # Move it into the default-role-order position.
        return True, f"added candidate part {target.name!r} (role={target.role})"

    return False, f"unknown action: {action}"


async def build(
    part_set: PartSet,
    intent: IntentSpec,
    *,
    max_iters: int = 5,
    annotate_fn=None,
) -> BuildResult:
    """Construct + verify + revise loop. Returns the final state
    (success=True if verification passed, False otherwise) along
    with the BuildJournal."""
    construct = _initial_layout(part_set)
    journal: list[BuildJournalEntry] = []
    journal.append(BuildJournalEntry(
        iteration=0, action="initial_layout",
        detail=f"laid out {len(construct.slots)} slots via DEFAULT_ORDER",
        construct_snapshot=construct.to_dict(),
    ))
    last_result: Optional[VerifierResult] = None

    for it in range(1, max_iters + 1):
        method = assess_methods(construct)
        # Hard abort: zero parts means nothing to verify.
        if not construct.slots:
            return BuildResult(
                success=False, final_construct=construct,
                final_verification=None, method_pick=method.pick,
                method_assessment=method, journal=journal,
                unresolved_diagnostics=[Diagnostic(
                    code="empty_construct",
                    detail="No slots left after revisions.",
                )],
            )
        last_result = await verify(construct, intent, annotate_fn=annotate_fn)
        snapshot = construct.to_dict()
        snapshot["method"] = method.to_dict()
        journal.append(BuildJournalEntry(
            iteration=it, action="verify",
            detail=(f"verifier {'PASS' if last_result.passed else 'FAIL'} — "
                    f"{len(last_result.diagnostics)} diagnostic(s); "
                    f"method_pick={method.pick}"),
            construct_snapshot=snapshot,
            verifier_result=last_result.to_dict(),
        ))
        if last_result.passed and method.pick is not None:
            return BuildResult(
                success=True, final_construct=construct,
                final_verification=last_result,
                method_pick=method.pick, method_assessment=method,
                journal=journal,
            )

        # Apply at most one fix per iteration so the journal stays
        # readable and we don't compound bad moves.
        applied_any = False
        for d in last_result.diagnostics:
            ok, detail = _apply_action(construct, d, part_set)
            journal.append(BuildJournalEntry(
                iteration=it, action="modify",
                detail=f"{d.code}: {'applied' if ok else 'skipped'} — {detail}",
            ))
            if ok:
                applied_any = True
                break
        if not applied_any:
            # No mechanical fix available — bail with the unresolved
            # diagnostics so the caller can hand back to the
            # interpreter for a different PartSet.
            return BuildResult(
                success=False, final_construct=construct,
                final_verification=last_result,
                method_pick=method.pick, method_assessment=method,
                journal=journal,
                unresolved_diagnostics=last_result.diagnostics,
            )

    # Out of iterations.
    method = assess_methods(construct)
    return BuildResult(
        success=False, final_construct=construct,
        final_verification=last_result,
        method_pick=method.pick, method_assessment=method,
        journal=journal,
        unresolved_diagnostics=(last_result.diagnostics if last_result else []),
    )
