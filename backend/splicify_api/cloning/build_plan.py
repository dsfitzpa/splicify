"""
Shared data structures for cloning build plans.

All operators produce plans built from these types, enabling apples-to-apples
comparison across Gibson, restriction, Golden Gate, etc.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class FragmentSourceType(str, Enum):
    """How a module's sequence will be obtained for assembly."""
    PCR_LIBRARY  = "pcr_library"    # PCR from module DB source plasmid (in-house)
    PCR_INVENTORY = "pcr_inventory"  # PCR from user-uploaded inventory plasmid
    PCR_ADDGENE  = "pcr_addgene"    # Order from Addgene, then PCR amplify
    SYNTHESIS    = "synthesis"       # Synthesize de novo (IDT, Twist, Azenta, etc.)
    EXISTING     = "existing"        # Use directly — linearized plasmid, no PCR
    UNKNOWN      = "unknown"         # Source not determinable from available data


@dataclass
class FragmentSource:
    """Sourcing information for a single module in the assembly."""
    module_index: int
    module_name: str
    source_type: FragmentSourceType

    # PCR-specific fields
    template_plasmid: Optional[str] = None   # filename / Addgene accession
    template_source_group: Optional[str] = None
    expected_amplicon_bp: int = 0
    pcr_difficulty: str = "easy"             # "easy" | "moderate" | "difficult"
    pcr_difficulty_reasons: List[str] = field(default_factory=list)

    # Synthesis-specific fields
    synthesis_bp: int = 0
    synthesis_tier: str = "standard"         # "standard" | "complex"
    synthesis_cost_estimate_usd: float = 0.0
    synthesis_days: float = 0.0

    notes: List[str] = field(default_factory=list)


@dataclass
class OverlapDesign:
    """
    Overlap/homology region and associated primer pair for a single junction
    in Gibson/HiFi assembly.

    Convention:
    - overlap_sequence: last overlap_len bp of the LEFT module's sequence.
      This appears naturally at the 3' end of the LEFT fragment's PCR product
      AND is added as a 5' tail on the RIGHT fragment's forward primer.
    - forward_primer: primer for the RIGHT module (5'-[overlap_seq]-[anneal]-3')
    - reverse_primer: primer for the LEFT module  (5'-RC[last anneal_len of left]-3')
    """
    junction_index: int
    left_module_name: str
    right_module_name: str

    # Overlap region
    overlap_sequence: str
    overlap_length: int
    overlap_tm: float
    overlap_gc: float

    # Primer sequences
    forward_primer: str        # for RIGHT module — contains overlap tail
    reverse_primer: str        # for LEFT module  — anneals to left module end

    # Annealing portions only (for Tm reporting)
    forward_primer_anneal: str
    reverse_primer_anneal: str
    forward_anneal_tm: float
    reverse_anneal_tm: float

    # Quality metrics
    uniqueness_score: float = 100.0  # 0-100; <50 = misassembly risk
    hairpin_dg_fwd: float = 0.0      # kcal/mol; < -3 = flag
    hairpin_dg_rev: float = 0.0
    self_dimer_dg_fwd: float = 0.0   # kcal/mol; < -6 = flag
    self_dimer_dg_rev: float = 0.0
    quality_score: float = 100.0     # 0-100 aggregate (from OverlapScorer)
    warnings: List[str] = field(default_factory=list)


@dataclass
class BuildStep:
    """A single step in the wet-lab protocol."""
    step_number: int
    step_type: str        # pcr | gel_purification | assembly | transformation |
                          # colony_screening | miniprep | sequencing | synthesis_wait
    description: str
    materials: List[str] = field(default_factory=list)
    estimated_hours: float = 0.0   # hands-on hours
    estimated_days: float = 0.0    # calendar days (including wait)


@dataclass
class OperatorMetrics:
    """
    Aggregate cost, time, labor, and risk for a complete build plan.
    Used for ranking across operator options.
    """
    # Counts
    primer_count: int = 0
    pcr_count: int = 0
    gel_count: int = 0
    assembly_count: int = 1
    transformation_count: int = 1
    miniprep_count: int = 0
    sequencing_count: int = 0

    # Cost breakdown (USD)
    primer_cost_usd: float = 0.0
    pcr_cost_usd: float = 0.0
    gel_cost_usd: float = 0.0
    assembly_cost_usd: float = 0.0
    transformation_cost_usd: float = 0.0
    miniprep_cost_usd: float = 0.0
    sequencing_cost_usd: float = 0.0
    synthesis_cost_usd: float = 0.0
    total_cost_usd: float = 0.0

    # Labor
    total_labor_hours: float = 0.0

    # Calendar time (critical path, days)
    total_calendar_days: float = 0.0

    # Risk (0.0 = low, 1.0 = high)
    pcr_risk: float = 0.0
    assembly_risk: float = 0.0
    overall_risk_score: float = 0.0
    risk_flags: List[str] = field(default_factory=list)


@dataclass
class GibsonBuildPlan:
    """
    Complete Gibson/HiFi Assembly build plan produced by GibsonOperator.evaluate().

    Contains everything needed to execute the cloning:
    - fragment_sources: how to obtain each PCR fragment
    - overlap_designs: primer pairs + quality metrics for each junction
    - primer_table: flat list for ordering (all primers in one place)
    - bom: bill of materials
    - steps: step-by-step protocol
    - metrics: cost/time/labor/risk summary
    """
    method: str = "gibson_hifi"
    feasible: bool = True
    infeasibility_reasons: List[str] = field(default_factory=list)

    fragment_count: int = 0
    assembly_topology: str = "circular"

    fragment_sources: List[FragmentSource] = field(default_factory=list)
    overlap_designs: List[OverlapDesign] = field(default_factory=list)

    # Flat primer table: one entry per primer (order-ready)
    # Keys: primer_name, sequence, tm_anneal, tm_overlap, length, gc_content,
    #       purpose, fragment, overlap_tail, annealing_portion,
    #       hairpin_dg, quality_score, warnings
    primer_table: List[Dict[str, Any]] = field(default_factory=list)

    bom: List[str] = field(default_factory=list)
    steps: List[BuildStep] = field(default_factory=list)

    metrics: OperatorMetrics = field(default_factory=OperatorMetrics)
    summary: str = ""
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to JSON-compatible dict."""
        return dataclasses.asdict(self)


@dataclass
class RestrictionJunctionPlan:
    """Restriction-cloning design details for a single junction."""
    junction_index: int
    left_module_name: str
    right_module_name: str
    left_enzyme: str
    right_enzyme: str
    strategy: str = "native_sites"   # native_sites | engineered_sites
    scar_sequence: Optional[str] = None
    internal_conflicts: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass
class RestrictionBuildPlan:
    """Complete restriction/ligation build plan produced by RestrictionOperator."""
    method: str = "restriction_cloning"
    feasible: bool = True
    infeasibility_reasons: List[str] = field(default_factory=list)

    fragment_count: int = 0
    assembly_topology: str = "circular"

    junction_plans: List[RestrictionJunctionPlan] = field(default_factory=list)
    engineered_primer_table: List[Dict[str, Any]] = field(default_factory=list)
    bom: List[str] = field(default_factory=list)
    steps: List[BuildStep] = field(default_factory=list)

    metrics: OperatorMetrics = field(default_factory=OperatorMetrics)
    summary: str = ""
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)


@dataclass
class SDMPrimerDesign:
    """Primer pair for a single Q5-style SDM edit."""
    template_name: str
    edit_start: int
    edit_end: int
    old_sequence: str
    new_sequence: str
    forward_primer: str
    reverse_primer: str
    forward_tm: float
    reverse_tm: float
    warnings: List[str] = field(default_factory=list)


@dataclass
class SDMBuildPlan:
    """Complete site-directed mutagenesis plan produced by SDMOperator."""
    method: str = "q5_sdm"
    feasible: bool = True
    infeasibility_reasons: List[str] = field(default_factory=list)

    template_name: str = ""
    primer_design: Optional[SDMPrimerDesign] = None
    bom: List[str] = field(default_factory=list)
    steps: List[BuildStep] = field(default_factory=list)

    metrics: OperatorMetrics = field(default_factory=OperatorMetrics)
    summary: str = ""
    warnings: List[str] = field(default_factory=list)
    
    # Additional fields for visualization
    mutated_sequence: str = ""
    mutation_type: str = ""  # deletion | insertion | substitution
    primer_strategy: str = ""  # back_to_back | single_primer | overlapping
    fwd_anneal_seq: str = ""
    rev_anneal_seq: str = ""
    overlap_tm: Optional[float] = None
    overlap_seq: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)


@dataclass
class GoldenGateJunctionPlan:
    """Golden Gate design details for a single junction."""
    junction_index: int
    left_module_name: str
    right_module_name: str
    overhang_4bp: str
    enzyme: str = "BsaI"
    strategy: str = "engineered_sites"  # native_sites | engineered_sites
    overhang_fidelity_score: float = 0.0
    warnings: List[str] = field(default_factory=list)


@dataclass
class GoldenGateBuildPlan:
    """Complete Golden Gate plan produced by GoldenGateOperator."""
    method: str = "golden_gate"
    feasible: bool = True
    infeasibility_reasons: List[str] = field(default_factory=list)

    fragment_count: int = 0
    assembly_topology: str = "circular"
    enzyme: str = "BsaI"
    strategy: str = "engineered_sites"
    domestication_burden: int = 0

    junction_plans: List[GoldenGateJunctionPlan] = field(default_factory=list)
    primer_table: List[Dict[str, Any]] = field(default_factory=list)
    bom: List[str] = field(default_factory=list)
    steps: List[BuildStep] = field(default_factory=list)

    metrics: OperatorMetrics = field(default_factory=OperatorMetrics)
    summary: str = ""
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)


# Import AttSiteMatch from gateway_sites
# (This will be imported when needed to avoid circular imports)


@dataclass
class GatewayJunctionPlan:
    """Gateway recombination design for a single junction."""
    junction_index: int
    left_module_name: str
    right_module_name: str
    left_att_site: str  # "attB1", "attL1", etc.
    right_att_site: str
    strategy: str  # "native_sites" | "pcr_add_sites"
    product_left_site: str  # "attL1" for BP reaction
    product_right_site: str
    reaction_type: str  # "BP" | "LR"
    ccdb_selection: bool = False
    warnings: List[str] = field(default_factory=list)


@dataclass
class GatewayBuildPlan:
    """Complete Gateway cloning build plan."""
    method: str = "gateway_cloning"
    feasible: bool = True
    infeasibility_reasons: List[str] = field(default_factory=list)

    fragment_count: int = 0
    assembly_topology: str = "circular"
    reaction_type: str = "BP"  # "BP" | "LR" | "MultiSite"

    junction_plans: List[GatewayJunctionPlan] = field(default_factory=list)
    primer_table: List[Dict[str, Any]] = field(default_factory=list)

    # Product information
    product_sequence: str = ""
    product_annotations: List[Dict[str, Any]] = field(default_factory=list)

    # Byproduct information
    byproduct_sequence: str = ""
    byproduct_description: str = ""

    bom: List[str] = field(default_factory=list)
    steps: List[BuildStep] = field(default_factory=list)
    metrics: OperatorMetrics = field(default_factory=OperatorMetrics)
    summary: str = ""
    warnings: List[str] = field(default_factory=list)

    # 2026-05-06: orientation-aware substrate classifications per input
    # module (one entry per compatible att pair found, with
    # kind='intermolecular'/'excision'/'inversion'). Surfaced by
    # GatewayOperator.evaluate.
    substrate_classifications: List[Dict[str, Any]] = field(default_factory=list)
    def to_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)
