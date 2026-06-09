"""
GoldenGateOperator — Phase 4 evaluator for Type IIS modular assembly.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Set

from ..utils import reverse_complement
from .build_plan import BuildStep, GoldenGateBuildPlan, GoldenGateJunctionPlan, OperatorMetrics
from .junction import build_junctions
from .lab_profile import DEFAULT_LAB_PROFILE, LabProfile


_TYPE_IIS = {
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "BbsI": "GAAGAC",
}

# Practical overhang panel with generally strong ligation behavior.
# Values are coarse fidelity scores (0-1) inspired by Potapov-style trends.
_OVERHANG_PANEL = [
    ("AATG", 0.96), ("GCTT", 0.95), ("GGAG", 0.94), ("CGCT", 0.94),
    ("TCGA", 0.93), ("CAGC", 0.92), ("ATCC", 0.92), ("GTCA", 0.91),
    ("AGGT", 0.90), ("TGCC", 0.90), ("CTTG", 0.89), ("GATG", 0.88),
]

# Avoid palindromes and simple homopolymers in 4bp overhangs.
_BAD_OVERHANGS = {
    "AAAA", "TTTT", "CCCC", "GGGG",
    "ATAT", "TATA", "CGCG", "GCGC",
}


class GoldenGateOperator:
    """Evaluate native or engineered Golden Gate assembly plan."""

    def __init__(self, lab_profile: Optional[LabProfile] = None, domestication_threshold: int = 6) -> None:
        self.lab = lab_profile or DEFAULT_LAB_PROFILE
        self.domestication_threshold = domestication_threshold

    def evaluate(self, modules: List[dict], topology: str = "circular") -> GoldenGateBuildPlan:
        plan = GoldenGateBuildPlan(assembly_topology=topology)
        real_modules = [m for m in modules if (m.get("sequence") or "").strip()]
        plan.fragment_count = len(real_modules)

        if len(real_modules) < 2:
            plan.feasible = False
            plan.infeasibility_reasons.append("Golden Gate requires at least 2 sequenced modules")
            return plan

        junctions = build_junctions(real_modules, topology=topology)
        enzyme = self._pick_enzyme(real_modules)
        plan.enzyme = enzyme

        domestication_burden = self._count_internal_sites(real_modules, _TYPE_IIS[enzyme])
        plan.domestication_burden = domestication_burden

        native_possible = self._native_mode_possible(junctions, enzyme)
        if native_possible:
            plan.strategy = "native_sites"
            plan.junction_plans = self._build_native_junctions(junctions, enzyme)
        else:
            plan.strategy = "engineered_sites"
            plan.junction_plans = self._build_engineered_junctions(junctions, enzyme)
            plan.primer_table = self._build_engineered_primers(real_modules, plan.junction_plans, enzyme)

        if domestication_burden > self.domestication_threshold:
            plan.warnings.append(
                f"High domestication burden: {domestication_burden} internal {enzyme} site(s)"
            )
            # keep feasible but raise risk/cost hard

        if len(plan.junction_plans) != len(junctions):
            plan.feasible = False
            plan.infeasibility_reasons.append("Unable to assign unique 4bp overhangs to all junctions")

        plan.bom = self._build_bom(plan)
        plan.steps = self._build_steps(plan)
        plan.metrics = self._compute_metrics(plan, len(junctions))
        plan.summary = self._build_summary(plan)
        return plan

    def _pick_enzyme(self, modules: List[dict]) -> str:
        # Prefer enzyme with fewer internal sites across all fragments.
        counts = {
            name: self._count_internal_sites(modules, site)
            for name, site in _TYPE_IIS.items()
        }
        return min(counts, key=counts.get)

    def _count_internal_sites(self, modules: List[dict], site: str) -> int:
        total = 0
        rc = reverse_complement(site)
        for mod in modules:
            seq = (mod.get("sequence") or "").upper()
            total += len(re.findall(site, seq))
            if rc != site:
                total += len(re.findall(rc, seq))
        return total

    def _native_mode_possible(self, junctions: List[object], enzyme: str) -> bool:
        for j in junctions:
            if enzyme not in j.internal_restriction_sites:
                return False
        return True

    def _build_native_junctions(self, junctions: List[object], enzyme: str) -> List[GoldenGateJunctionPlan]:
        plans: List[GoldenGateJunctionPlan] = []
        used: Set[str] = set()
        for j in junctions:
            overhang = self._choose_overhang(used)
            if overhang is None:
                return plans
            used.add(overhang[0])
            plans.append(
                GoldenGateJunctionPlan(
                    junction_index=j.junction_index,
                    left_module_name=j.left_module_name,
                    right_module_name=j.right_module_name,
                    overhang_4bp=overhang[0],
                    enzyme=enzyme,
                    strategy="native_sites",
                    overhang_fidelity_score=overhang[1],
                    warnings=[],
                )
            )
        return plans

    def _build_engineered_junctions(self, junctions: List[object], enzyme: str) -> List[GoldenGateJunctionPlan]:
        plans: List[GoldenGateJunctionPlan] = []
        used: Set[str] = set()
        for j in junctions:
            overhang = self._choose_overhang(used)
            if overhang is None:
                return plans
            used.add(overhang[0])
            warnings = []
            if j.reading_frame_continuation and overhang[0] in {"AATG", "GGAG"}:
                warnings.append("Check frame continuity at coding junction")
            if j.regulatory_boundary:
                warnings.append("Verify boundary overhang does not alter UTR context")
            plans.append(
                GoldenGateJunctionPlan(
                    junction_index=j.junction_index,
                    left_module_name=j.left_module_name,
                    right_module_name=j.right_module_name,
                    overhang_4bp=overhang[0],
                    enzyme=enzyme,
                    strategy="engineered_sites",
                    overhang_fidelity_score=overhang[1],
                    warnings=warnings,
                )
            )
        return plans

    def _choose_overhang(self, used: Set[str]) -> Optional[tuple[str, float]]:
        for oh, score in _OVERHANG_PANEL:
            if oh in used or oh in _BAD_OVERHANGS:
                continue
            if reverse_complement(oh) in used:
                continue
            return (oh, score)
        return None

    def _build_engineered_primers(
        self,
        modules: List[dict],
        junction_plans: List[GoldenGateJunctionPlan],
        enzyme: str,
    ) -> List[Dict[str, object]]:
        site = _TYPE_IIS[enzyme]
        primers: List[Dict[str, object]] = []
        for jp in junction_plans:
            left_seq = (modules[jp.junction_index].get("sequence") or "").upper()
            right_seq = (modules[(jp.junction_index + 1) % len(modules)].get("sequence") or "").upper()
            left_anneal = reverse_complement(left_seq[-22:]) if left_seq else "N" * 22
            right_anneal = right_seq[:22] if right_seq else "N" * 22

            primers.append({
                "primer_name": f"GG_J{jp.junction_index}_{jp.left_module_name}_REV",
                "sequence": f"{site}{jp.overhang_4bp}{left_anneal}",
                "enzyme": enzyme,
                "overhang_4bp": jp.overhang_4bp,
                "purpose": "golden_gate_tail",
            })
            primers.append({
                "primer_name": f"GG_J{jp.junction_index}_{jp.right_module_name}_FWD",
                "sequence": f"{site}{jp.overhang_4bp}{right_anneal}",
                "enzyme": enzyme,
                "overhang_4bp": jp.overhang_4bp,
                "purpose": "golden_gate_tail",
            })
        return primers

    def _build_bom(self, plan: GoldenGateBuildPlan) -> List[str]:
        bom = [
            f"{plan.enzyme} Type IIS enzyme",
            "T4 DNA ligase",
            "Golden Gate cycling buffer",
            "Competent cells",
            "Selection plates",
            "Sanger sequencing",
        ]
        if plan.primer_table:
            bom.append(f"Engineered Type IIS primers ({len(plan.primer_table)} total)")
        return bom

    def _build_steps(self, plan: GoldenGateBuildPlan) -> List[BuildStep]:
        steps: List[BuildStep] = []
        n = 1
        if plan.primer_table:
            steps.append(BuildStep(
                step_number=n,
                step_type="pcr",
                description="PCR amplify modules with Type IIS tails + overhangs",
                materials=["Q5 polymerase", "Golden Gate primers"],
                estimated_hours=2.0,
                estimated_days=0.0,
            ))
            n += 1

        steps.extend([
            BuildStep(n, "digestion_ligation", "One-pot Golden Gate cycling digest/ligation", [plan.enzyme, "T4 ligase"], 1.5, 0.0),
            BuildStep(n + 1, "transformation", "Transform assembly product", ["competent cells"], 0.8, 1.0),
            BuildStep(n + 2, "colony_screening", "Screen colonies", ["colony PCR reagents"], 1.2, 0.0),
            BuildStep(n + 3, "sequencing", "Sequence verify junctions", ["Sanger sequencing"], 0.3, 1.0),
        ])
        return steps

    def _compute_metrics(self, plan: GoldenGateBuildPlan, junction_count: int) -> OperatorMetrics:
        metrics = OperatorMetrics()
        metrics.primer_count = len(plan.primer_table)
        metrics.pcr_count = 1 if plan.primer_table else 0
        metrics.assembly_count = 1
        metrics.transformation_count = 1
        metrics.sequencing_count = max(1, min(4, junction_count))

        # Length-scaled primer cost. Golden Gate primers carry the BsaI
        # site + 4-nt overhang + anneal: typical 35-45 mer.
        primer_lens = []
        for p in plan.primer_table:
            seq = p.get("sequence") or p.get("forward_primer") or p.get("reverse_primer") or ""
            primer_lens.append(len(seq) or 40)
        if len(primer_lens) < metrics.primer_count:
            primer_lens += [40] * (metrics.primer_count - len(primer_lens))
        metrics.primer_cost_usd = sum(self.lab.primer_cost(L) for L in primer_lens)

        metrics.pcr_cost_usd = metrics.pcr_count * self.lab.pcr_rxn_cost_usd
        # NEB Golden Gate Assembly Mix per rxn (lab catalog: neb_gg_assembly).
        metrics.assembly_cost_usd = self.lab.catalog["neb_gg_assembly"].cost_per_rxn
        metrics.transformation_cost_usd = self.lab.transformation_cost_usd
        # ONT: one read covers the whole plasmid.
        metrics.sequencing_count = self.lab.sequencing_reads_per_construct
        metrics.sequencing_cost_usd = metrics.sequencing_count * self.lab.sequencing_cost_usd

        metrics.total_cost_usd = (
            metrics.primer_cost_usd
            + metrics.pcr_cost_usd
            + metrics.assembly_cost_usd
            + metrics.transformation_cost_usd
            + metrics.sequencing_cost_usd
        )

        metrics.total_labor_hours = 4.2 + (1.8 if plan.primer_table else 0.0)
        metrics.total_calendar_days = 2.2 + (1.0 if plan.primer_table else 0.0)

        mean_fidelity = 0.0
        if plan.junction_plans:
            mean_fidelity = sum(jp.overhang_fidelity_score for jp in plan.junction_plans) / len(plan.junction_plans)

        metrics.pcr_risk = 0.08 + (0.08 if plan.primer_table else 0.0)
        metrics.assembly_risk = max(0.05, 0.20 - (mean_fidelity * 0.12) + (plan.domestication_burden * 0.015))
        metrics.overall_risk_score = min(1.0, (metrics.pcr_risk + metrics.assembly_risk) / 2.0)

        if plan.domestication_burden:
            metrics.risk_flags.append(f"{plan.domestication_burden} internal {plan.enzyme} site(s) to domesticate")
        if mean_fidelity < 0.90:
            metrics.risk_flags.append("Low average overhang fidelity")

        return metrics

    def _build_summary(self, plan: GoldenGateBuildPlan) -> str:
        return (
            f"Golden Gate ({plan.enzyme}) {plan.strategy} plan across {len(plan.junction_plans)} junctions; "
            f"domestication burden {plan.domestication_burden}; "
            f"estimated ${plan.metrics.total_cost_usd:.2f} and {plan.metrics.total_calendar_days:.1f} days."
        )
