"""
RestrictionOperator — Phase 2 evaluator for restriction/ligation cloning.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Set, Tuple

from ..utils import reverse_complement
from .build_plan import BuildStep, OperatorMetrics, RestrictionBuildPlan, RestrictionJunctionPlan
from .junction import Junction, build_junctions
from .lab_profile import DEFAULT_LAB_PROFILE, LabProfile
from .re_database import RE_DATABASE, RestrictionEnzyme, are_compatible, compute_scar


_PREFERRED_ENGINEERED_ORDER = [
    "EcoRI", "XhoI", "BamHI", "HindIII", "NheI", "AgeI", "MluI", "KpnI", "BglII", "NotI",
]


class RestrictionOperator:
    """Evaluate feasibility and expected effort for restriction cloning."""

    def __init__(self, lab_profile: Optional[LabProfile] = None) -> None:
        self.lab = lab_profile or DEFAULT_LAB_PROFILE

    def evaluate(self, modules: List[dict], topology: str = "circular") -> RestrictionBuildPlan:
        plan = RestrictionBuildPlan(assembly_topology=topology)
        real_modules = [m for m in modules if (m.get("sequence") or "").strip()]
        plan.fragment_count = len(real_modules)

        if len(real_modules) < 2:
            plan.feasible = False
            plan.infeasibility_reasons.append("Restriction cloning requires at least 2 sequenced modules")
            return plan

        junctions = build_junctions(real_modules, topology=topology)
        used_enzymes: Set[str] = set()
        engineered_primers: List[Dict[str, object]] = []

        for j in junctions:
            left_mod = real_modules[j.left_module_index]
            right_mod = real_modules[j.right_module_index]

            left_native = self._native_site_candidates(j, side="left")
            right_native = self._native_site_candidates(j, side="right")
            conflicts = self._collect_conflicts(left_mod, right_mod, left_native, right_native)

            pair = self._select_pair(left_native, right_native, used_enzymes)
            strategy = "native_sites"

            if pair is None:
                pair = self._pick_engineered_pair(used_enzymes)
                strategy = "engineered_sites"
                if pair is None:
                    plan.feasible = False
                    plan.infeasibility_reasons.append(
                        f"J{j.junction_index}: unable to assign a compatible enzyme pair"
                    )
                    continue
                engineered_primers.extend(self._engineered_primers(j, left_mod, right_mod, pair))

            left_enz, right_enz = pair
            used_enzymes.update([left_enz.name, right_enz.name])
            scar = compute_scar(left_enz, right_enz)

            warnings: List[str] = []
            if scar is None:
                plan.feasible = False
                plan.infeasibility_reasons.append(
                    f"J{j.junction_index}: incompatible overhangs ({left_enz.name}/{right_enz.name})"
                )
            elif j.reading_frame_continuation and scar and (len(scar) % 3 != 0):
                warnings.append("Scar length not frame-safe for coding fusion")
            if j.regulatory_boundary and scar:
                warnings.append("Scar at regulatory boundary may alter expression")

            plan.junction_plans.append(
                RestrictionJunctionPlan(
                    junction_index=j.junction_index,
                    left_module_name=j.left_module_name,
                    right_module_name=j.right_module_name,
                    left_enzyme=left_enz.name,
                    right_enzyme=right_enz.name,
                    strategy=strategy,
                    scar_sequence=scar,
                    internal_conflicts=conflicts,
                    warnings=warnings,
                )
            )
            for warning in warnings:
                plan.warnings.append(f"J{j.junction_index}: {warning}")
            for conflict in conflicts:
                plan.warnings.append(f"J{j.junction_index}: {conflict}")

        plan.engineered_primer_table = engineered_primers
        plan.bom = self._build_bom(plan)
        plan.steps = self._build_steps(plan)
        plan.metrics = self._compute_metrics(real_modules, plan, len(junctions))
        plan.summary = self._build_summary(plan)

        return plan

    def _native_site_candidates(self, junction: Junction, side: str) -> List[RestrictionEnzyme]:
        candidates: List[RestrictionEnzyme] = []
        boundary = len(junction.left_flank)
        for name, positions in junction.internal_restriction_sites.items():
            enzyme = RE_DATABASE.get(name)
            if not enzyme:
                continue
            if side == "left" and any(pos < boundary for pos in positions):
                candidates.append(enzyme)
            if side == "right" and any(pos >= boundary for pos in positions):
                candidates.append(enzyme)
        # keep deterministic ordering by preferred list first
        order = {name: i for i, name in enumerate(_PREFERRED_ENGINEERED_ORDER)}
        return sorted(candidates, key=lambda e: order.get(e.name, 999))

    def _collect_conflicts(
        self,
        left_mod: dict,
        right_mod: dict,
        left_native: List[RestrictionEnzyme],
        right_native: List[RestrictionEnzyme],
    ) -> List[str]:
        conflicts: List[str] = []
        for enzyme in left_native:
            hits = self._count_sites(left_mod.get("sequence", ""), enzyme)
            if hits > 1:
                conflicts.append(f"{enzyme.name} appears {hits}x in left fragment")
        for enzyme in right_native:
            hits = self._count_sites(right_mod.get("sequence", ""), enzyme)
            if hits > 1:
                conflicts.append(f"{enzyme.name} appears {hits}x in right fragment")
        return conflicts

    def _select_pair(
        self,
        left_native: List[RestrictionEnzyme],
        right_native: List[RestrictionEnzyme],
        used_enzymes: Set[str],
    ) -> Optional[Tuple[RestrictionEnzyme, RestrictionEnzyme]]:
        # Prefer pairs with no reuse and compatible overhangs.
        best: Optional[Tuple[RestrictionEnzyme, RestrictionEnzyme]] = None
        best_score = -999
        for a in left_native:
            for b in right_native:
                if not are_compatible(a, b):
                    continue
                score = 0
                if a.name not in used_enzymes:
                    score += 3
                if b.name not in used_enzymes:
                    score += 3
                if a.name == b.name:
                    score -= 1
                score -= int(a.star_activity) + int(b.star_activity)
                if score > best_score:
                    best = (a, b)
                    best_score = score
        return best

    def _pick_engineered_pair(
        self,
        used_enzymes: Set[str],
    ) -> Optional[Tuple[RestrictionEnzyme, RestrictionEnzyme]]:
        pool = [RE_DATABASE[name] for name in _PREFERRED_ENGINEERED_ORDER if name in RE_DATABASE]
        best: Optional[Tuple[RestrictionEnzyme, RestrictionEnzyme]] = None
        best_score = -999
        for a in pool:
            for b in pool:
                if a.name == b.name:
                    continue
                if not are_compatible(a, b):
                    continue
                score = 0
                if a.name not in used_enzymes:
                    score += 2
                if b.name not in used_enzymes:
                    score += 2
                if a.buffer == b.buffer:
                    score += 1
                if score > best_score:
                    best = (a, b)
                    best_score = score
        return best

    def _engineered_primers(
        self,
        junction: Junction,
        left_mod: dict,
        right_mod: dict,
        pair: Tuple[RestrictionEnzyme, RestrictionEnzyme],
    ) -> List[Dict[str, object]]:
        left_enz, right_enz = pair
        left_seq = (left_mod.get("sequence") or "").upper()
        right_seq = (right_mod.get("sequence") or "").upper()

        left_anneal = reverse_complement(left_seq[-22:]) if left_seq else "N" * 22
        right_anneal = right_seq[:22] if right_seq else "N" * 22

        return [
            {
                "primer_name": f"J{junction.junction_index}_{junction.left_module_name}_REV_eng",
                "sequence": f"{left_enz.recognition_seq}{left_anneal}",
                "purpose": "engineered_re_site",
                "enzyme": left_enz.name,
            },
            {
                "primer_name": f"J{junction.junction_index}_{junction.right_module_name}_FWD_eng",
                "sequence": f"{right_enz.recognition_seq}{right_anneal}",
                "purpose": "engineered_re_site",
                "enzyme": right_enz.name,
            },
        ]

    def _count_sites(self, seq: str, enzyme: RestrictionEnzyme) -> int:
        if not seq:
            return 0
        s = seq.upper()
        pattern = enzyme.recognition_seq.upper()
        rc = reverse_complement(pattern)
        hits = len([m.start() for m in re.finditer(pattern, s)])
        if rc != pattern:
            hits += len([m.start() for m in re.finditer(rc, s)])
        return hits

    def _build_bom(self, plan: RestrictionBuildPlan) -> List[str]:
        unique_enzymes = sorted({jp.left_enzyme for jp in plan.junction_plans} | {jp.right_enzyme for jp in plan.junction_plans})
        bom = [
            "Restriction digest buffer",
            "T4 DNA ligase + ligase buffer",
            "Competent cells",
            "LB agar + selection antibiotic",
        ]
        if unique_enzymes:
            bom.append("Restriction enzymes: " + ", ".join(unique_enzymes))
        if plan.engineered_primer_table:
            bom.append(f"PCR primers for engineered sites ({len(plan.engineered_primer_table)} total)")
        return bom

    def _build_steps(self, plan: RestrictionBuildPlan) -> List[BuildStep]:
        steps: List[BuildStep] = []
        n = 1
        if plan.engineered_primer_table:
            steps.append(BuildStep(
                step_number=n,
                step_type="pcr",
                description="PCR amplify fragments with engineered restriction-site tails",
                materials=["Q5 polymerase", "engineered primers"],
                estimated_hours=2.0,
                estimated_days=0.0,
            ))
            n += 1
        steps.append(BuildStep(
            step_number=n,
            step_type="digestion",
            description="Digest vector/insert fragments with assigned restriction enzymes",
            materials=["restriction enzymes", "digest buffer"],
            estimated_hours=1.5,
            estimated_days=0.0,
        ))
        n += 1
        steps.append(BuildStep(
            step_number=n,
            step_type="gel_purification",
            description="Gel purify digested fragments",
            materials=["agarose gel", "gel extraction kit"],
            estimated_hours=1.0,
            estimated_days=0.0,
        ))
        n += 1
        steps.append(BuildStep(
            step_number=n,
            step_type="ligation",
            description="Ligate compatible overhang fragments",
            materials=["T4 DNA ligase"],
            estimated_hours=0.8,
            estimated_days=0.0,
        ))
        n += 1
        steps.append(BuildStep(
            step_number=n,
            step_type="transformation",
            description="Transform ligation product and plate colonies",
            materials=["competent cells", "selection plates"],
            estimated_hours=1.0,
            estimated_days=1.0,
        ))
        n += 1
        steps.append(BuildStep(
            step_number=n,
            step_type="colony_screening",
            description="Screen colonies by colony PCR or digest",
            materials=["screening primers"],
            estimated_hours=1.5,
            estimated_days=0.0,
        ))
        n += 1
        steps.append(BuildStep(
            step_number=n,
            step_type="sequencing",
            description="Sequence-verify final construct",
            materials=["Sanger sequencing"],
            estimated_hours=0.3,
            estimated_days=1.0,
        ))
        return steps

    def _compute_metrics(self, modules: List[dict], plan: RestrictionBuildPlan, junction_count: int) -> OperatorMetrics:
        metrics = OperatorMetrics()
        metrics.primer_count = len(plan.engineered_primer_table)
        metrics.pcr_count = 1 if plan.engineered_primer_table else 0
        metrics.gel_count = 1
        metrics.assembly_count = 1
        metrics.transformation_count = 1
        metrics.miniprep_count = 1
        metrics.sequencing_count = max(1, min(3, junction_count))

        # Length-scaled primer cost: each engineered primer carries the
        # actual oligo length when available, falling back to a 35 mer
        # estimate (RE-site primers run ~30-40 mer with the cut site +
        # flanks). $0.24/bp aligns with IDT 25 nmol desalted.
        primer_lens = [len(p.get("sequence", "")) or 35
                       for p in plan.engineered_primer_table]
        metrics.primer_cost_usd = sum(self.lab.primer_cost(L) for L in primer_lens)
        metrics.pcr_cost_usd = metrics.pcr_count * self.lab.pcr_rxn_cost_usd
        metrics.gel_cost_usd = metrics.gel_count * self.lab.gel_lane_cost_usd

        unique_enzymes = {jp.left_enzyme for jp in plan.junction_plans} | {jp.right_enzyme for jp in plan.junction_plans}
        digest_cost = sum(RE_DATABASE[e].cost_per_rxn_usd for e in unique_enzymes if e in RE_DATABASE)
        ligation_cost = 3.5
        metrics.assembly_cost_usd = digest_cost + ligation_cost

        metrics.transformation_cost_usd = self.lab.transformation_cost_usd
        metrics.miniprep_cost_usd = self.lab.miniprep_cost_usd
        # ONT sequencing: one full-plasmid read covers all junctions; cap
        # to 1 read /construct regardless of junction_count.
        metrics.sequencing_count = self.lab.sequencing_reads_per_construct
        metrics.sequencing_cost_usd = metrics.sequencing_count * self.lab.sequencing_cost_usd

        metrics.total_cost_usd = (
            metrics.primer_cost_usd
            + metrics.pcr_cost_usd
            + metrics.gel_cost_usd
            + metrics.assembly_cost_usd
            + metrics.transformation_cost_usd
            + metrics.miniprep_cost_usd
            + metrics.sequencing_cost_usd
        )

        metrics.total_labor_hours = 6.0 + (2.0 if plan.engineered_primer_table else 0.0)
        metrics.total_calendar_days = 3.0 + (1.0 if plan.engineered_primer_table else 0.0)

        conflict_count = sum(len(j.internal_conflicts) for j in plan.junction_plans)
        engineered_junctions = sum(1 for j in plan.junction_plans if j.strategy == "engineered_sites")
        metrics.pcr_risk = 0.10 + (0.08 if engineered_junctions else 0.0)
        metrics.assembly_risk = 0.10 + 0.03 * conflict_count + 0.04 * engineered_junctions
        metrics.overall_risk_score = min(1.0, (metrics.pcr_risk + metrics.assembly_risk) / 2.0)

        if conflict_count:
            metrics.risk_flags.append(f"{conflict_count} internal restriction-site conflicts")
        if engineered_junctions:
            metrics.risk_flags.append(f"{engineered_junctions} junction(s) require engineered sites")

        return metrics

    def _build_summary(self, plan: RestrictionBuildPlan) -> str:
        if not plan.junction_plans:
            return "No restriction junction assignments were generated"
        engineered = sum(1 for j in plan.junction_plans if j.strategy == "engineered_sites")
        return (
            f"Restriction cloning plan across {len(plan.junction_plans)} junctions; "
            f"{engineered} use engineered PCR-tail sites; "
            f"estimated ${plan.metrics.total_cost_usd:.2f} and {plan.metrics.total_calendar_days:.1f} days."
        )
