"""
Splicify Cloning Operator Library
==================================
Evaluates cloning strategies for a resolved plasmid design.

Architecture:
- Junction: atomic boundary between two adjacent modules (foundation)
- Operator: "compiler backend" that takes a ModuleGraph → CandidatePlan
  - GibsonOperator  (Phase 1, implemented)
  - RestrictionOperator  (Phase 2, implemented)
  - DiffRouter / SDMOperator  (Phase 3, implemented)
  - GoldenGateOperator  (Phase 4, implemented)

Usage:
    from splicify_api.cloning import GibsonOperator, build_junctions, LabProfile

    junctions = build_junctions(resolved_modules, topology="circular")
    op = GibsonOperator()
    plan = op.evaluate(resolved_modules, topology="circular")
    print(plan.metrics.total_cost_usd)
"""
from .junction import Junction, build_junctions
from .lab_profile import LabProfile, DEFAULT_LAB_PROFILE
from .build_plan import (
    FragmentSource, FragmentSourceType,
    OverlapDesign, BuildStep, OperatorMetrics, GibsonBuildPlan,
    RestrictionJunctionPlan, RestrictionBuildPlan,
    SDMPrimerDesign, SDMBuildPlan,
    GoldenGateJunctionPlan, GoldenGateBuildPlan,
)
from .gibson_operator import GibsonOperator
from .restriction_operator import RestrictionOperator
from .diff_router import DiffRouter, EditType, ModuleChange
from .sdm_operator import SDMOperator
from .golden_gate_operator import GoldenGateOperator
from .re_database import RestrictionEnzyme, RE_DATABASE, COMPATIBLE_ENZYME_PAIRS, compute_scar

__all__ = [
    "Junction", "build_junctions",
    "LabProfile", "DEFAULT_LAB_PROFILE",
    "FragmentSource", "FragmentSourceType",
    "OverlapDesign", "BuildStep", "OperatorMetrics", "GibsonBuildPlan",
    "RestrictionJunctionPlan", "RestrictionBuildPlan",
    "SDMPrimerDesign", "SDMBuildPlan",
    "GoldenGateJunctionPlan", "GoldenGateBuildPlan",
    "GibsonOperator",
    "RestrictionOperator",
    "DiffRouter", "EditType", "ModuleChange",
    "SDMOperator",
    "GoldenGateOperator",
    "RestrictionEnzyme", "RE_DATABASE", "COMPATIBLE_ENZYME_PAIRS", "compute_scar",
]
