"""
LabProfile — configurable cost and labor parameters for a wet lab.

Pricing model (2026-05 refresh):
  - Reagents are stored as bulk catalog entries (price + units per pack);
    per-reaction cost is derived. Each workflow declares the bundle of
    reagents it consumes per build so a single edit to the catalog
    propagates to every workflow.
  - Primers are length-scaled at $0.24/bp (25 nmol scale, IDT desalted).
  - Synthesis is tiered by length (0.5-1.8 / 1.8-3.2 / 3.2-5.0 kbp
    bands @ $0.07 / $0.08 / $0.09 per bp respectively).
  - Sequencing defaults to one ONT (long-read) read per construct at $15.

The legacy flat *_cost_usd attributes are still exposed for backwards
compatibility with operators that read them directly; they are derived
from the reagent catalog at construction time.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Tuple


@dataclass
class Reagent:
    """A bulk-pack reagent. cost_per_rxn = bulk_price_usd / bulk_rxns."""
    name: str
    bulk_price_usd: float
    bulk_rxns: int                          # reactions per pack
    usage_note: str = ""                    # human-readable usage detail
    catalog_no: str = ""                    # vendor SKU when known

    @property
    def cost_per_rxn(self) -> float:
        return self.bulk_price_usd / max(1, self.bulk_rxns)


# ---------------------------------------------------------------------------
# Bulk catalog — single source of truth for reagent prices.
# ---------------------------------------------------------------------------
# Numbers come from the user's 2026-02 cost worksheet and are derived from
# vendor catalog list prices for academic-lab-sized packs.
DEFAULT_CATALOG: Dict[str, Reagent] = {
    # Cells / plates
    "neb_10beta": Reagent(
        "NEB 10-beta competent E. coli", 282.00, 20,
        "1 transformation per pack rxn", "C3019H",
    ),
    "stbl3": Reagent(
        "Stbl3 competent E. coli (lentiviral)", 250.00, 20,
        "1 transformation; required for LTR-containing payloads",
    ),
    "lb_agar_plate": Reagent(
        "LB-agar selection plate (poured)", 51.00, 20,
        "1 plate per transformation",
    ),
    # Gels + ladders
    "ultrapure_agarose_g": Reagent(
        "UltraPure agarose (per gram)", 824.65, 500,
        "0.6 g / 60 mL gel = 1.0 g per construct (rounded up)",
    ),
    "sybr_safe_ul": Reagent(
        "SYBR Safe DNA Gel Stain (per µL)", 113.00, 400,
        "6 µL / 60 mL gel",
    ),
    "ladder_lane": Reagent(
        "1 kb Plus DNA Ladder (per lane)", 268.00, 1000,
        "1 lane per gel",
    ),
    # PCR
    "neb_q5_rxn": Reagent(
        "NEB Q5 Hot Start (per rxn)", 898.00, 500,
        "1 rxn / 25 µL PCR",
    ),
    # Gibson
    "nebuilder_hifi": Reagent(
        "NEBuilder HiFi DNA Assembly Master Mix (E2621X)", 2774.00, 250,
        "1 rxn / 20 µL assembly",
    ),
    # Gateway
    "bp_clonase_ii": Reagent(
        "Gateway BP Clonase II Enzyme Mix", 700.00, 20,
        "1 rxn / 10 µL recombination", "11789-020",
    ),
    "lr_clonase_ii": Reagent(
        "Gateway LR Clonase II Enzyme Mix", 700.00, 20,
        "1 rxn / 10 µL recombination", "11791-020",
    ),
    # Restriction / Golden Gate / sgRNA
    "neb_re_pair_rxn": Reagent(
        "NEB restriction enzyme pair (per double digest)", 350.00, 100,
        "20 U each enzyme + buffer for 1 µg DNA",
    ),
    "t4_ligase_rxn": Reagent(
        "NEB T4 DNA Ligase + 10X buffer (per rxn)", 200.00, 100,
        "1 rxn / 20 µL ligation",
    ),
    "neb_gg_assembly": Reagent(
        "NEB Golden Gate Assembly Mix (BsaI-HF v2; E1601)", 312.00, 20,
        "1 rxn / 25 µL one-pot Golden Gate",
    ),
    "esp3i_rxn": Reagent(
        "NEB Esp3I (BsmBI v2) (per rxn)", 300.00, 250,
        "1 rxn / 25 µL Type IIs digest (sgRNA-into-lentiCRISPR)",
    ),
    # SDM
    "q5_sdm_kit": Reagent(
        "NEB Q5 Site-Directed Mutagenesis Kit (E0554S)", 250.00, 10,
        "1 rxn = PCR + KLD + transformation control", "E0554S",
    ),
    # Sequencing — ONT default (Plasmidsaurus / Primordium / similar)
    "ont_seq_read": Reagent(
        "ONT plasmid sequencing read (Plasmidsaurus-class)", 15.00, 1,
        "1 read = full coverage of one ~10 kb plasmid",
    ),
}


# ---------------------------------------------------------------------------
# Primer + synthesis pricing functions
# ---------------------------------------------------------------------------
PRIMER_COST_PER_BP_USD: float = 0.24

# Tiered gene-synthesis pricing. Each entry is (max_kbp_inclusive, $/bp).
# Bands taken from the user's 2026-05 vendor quote; <500 bp uses the
# lowest tier and >5000 bp uses the highest (extrapolation).
SYNTHESIS_TIERS: List[Tuple[float, float]] = [
    (1800.0, 0.07),
    (3200.0, 0.08),
    (5000.0, 0.09),
]


def primer_cost(length_bp: int) -> float:
    """Cost of a single ordered primer at standard 25 nmol scale."""
    return max(0, length_bp) * PRIMER_COST_PER_BP_USD


def synthesis_cost_per_bp(length_bp: int) -> float:
    """Resolve the per-bp synthesis tier for a fragment of `length_bp` bp."""
    for upper_bp, rate in SYNTHESIS_TIERS:
        if length_bp <= upper_bp:
            return rate
    return SYNTHESIS_TIERS[-1][1]


def synthesis_cost(length_bp: int) -> float:
    """Total synthesis cost for one ordered fragment."""
    return length_bp * synthesis_cost_per_bp(length_bp)


# ---------------------------------------------------------------------------
# Per-workflow reagent bundles
# ---------------------------------------------------------------------------
# Each bundle is a list of (catalog_key, qty_rxns_per_construct) tuples. The
# qty multiplier is an integer (or float) to handle reagents that scale with
# fragment count — e.g. Gibson uses 1× Q5 per fragment so the operator
# overrides q5 quantity from its plan; the bundle declares the BASELINE.
WORKFLOW_REAGENT_BUNDLES: Dict[str, List[Tuple[str, float]]] = {
    "gibson": [
        ("neb_10beta", 1),
        ("lb_agar_plate", 1),
        ("nebuilder_hifi", 1),
        ("ultrapure_agarose_g", 0.6),
        ("sybr_safe_ul", 6),
        ("ladder_lane", 1),
    ],
    "gateway": [
        ("neb_10beta", 1),
        ("lb_agar_plate", 1),
        ("bp_clonase_ii", 1),    # BP for entry; LR added separately when both fire
        ("ultrapure_agarose_g", 0.6),
        ("sybr_safe_ul", 6),
        ("ladder_lane", 1),
    ],
    "restriction": [
        ("neb_10beta", 1),
        ("lb_agar_plate", 1),
        ("neb_re_pair_rxn", 1),
        ("t4_ligase_rxn", 1),
        ("ultrapure_agarose_g", 0.6),
        ("sybr_safe_ul", 6),
        ("ladder_lane", 1),
    ],
    "golden_gate": [
        ("neb_10beta", 1),
        ("lb_agar_plate", 1),
        ("neb_gg_assembly", 1),
        ("ultrapure_agarose_g", 0.6),
        ("sybr_safe_ul", 6),
        ("ladder_lane", 1),
    ],
    "sgrna_golden_gate": [
        ("stbl3", 1),
        ("lb_agar_plate", 1),
        ("esp3i_rxn", 1),
        ("t4_ligase_rxn", 1),
        ("ultrapure_agarose_g", 0.6),
        ("sybr_safe_ul", 6),
        ("ladder_lane", 1),
    ],
    "sdm": [
        ("neb_10beta", 1),
        ("lb_agar_plate", 1),
        ("q5_sdm_kit", 1),       # kit bundles enzyme + KLD + buffers
    ],
}


def reagent_lines(workflow: str, catalog: Dict[str, Reagent] = DEFAULT_CATALOG) -> List[Tuple[Reagent, float]]:
    """Resolve a workflow's reagent bundle into (Reagent, multiplier) pairs."""
    bundle = WORKFLOW_REAGENT_BUNDLES.get(workflow, [])
    out: List[Tuple[Reagent, float]] = []
    for key, qty in bundle:
        r = catalog.get(key)
        if r is None:
            continue
        out.append((r, float(qty)))
    return out


def workflow_reagent_cost_per_construct(workflow: str, catalog: Dict[str, Reagent] = DEFAULT_CATALOG) -> float:
    """Per-construct reagent cost for a workflow (excludes primers/synth/seq)."""
    return sum(r.cost_per_rxn * qty for r, qty in reagent_lines(workflow, catalog))



# ---------------------------------------------------------------------------
# Quick per-construct cost estimator for the predesign router
# ---------------------------------------------------------------------------
def estimate_workflow_cost(
    workflow: str,
    n_parts: int,
    *,
    catalog: Dict[str, Reagent] = DEFAULT_CATALOG,
    avg_primer_len: int = 35,
    primers_per_part: int = 2,
    avg_synthesis_bp: int = 0,
) -> float:
    """Approximate total per-construct cost without running the full operator.

    Covers primers (length-scaled), per-fragment Q5 PCR, the workflow reagent
    bundle (cells/plate/assembly/gel-share), and one ONT sequencing read.
    Used by predesign/cloning_router._evaluate_workflow to surface a
    realistic cost on workflow_trace.txt before any operator has actually
    designed primers.
    """
    primer_total = primers_per_part * n_parts * primer_cost(avg_primer_len)
    pcr_total = n_parts * catalog["neb_q5_rxn"].cost_per_rxn
    reagent_total = workflow_reagent_cost_per_construct(workflow, catalog)
    seq_total = catalog["ont_seq_read"].cost_per_rxn
    synth_total = synthesis_cost(avg_synthesis_bp) if avg_synthesis_bp > 0 else 0.0
    return primer_total + pcr_total + reagent_total + seq_total + synth_total


# ---------------------------------------------------------------------------
# LabProfile — main configurable container
# ---------------------------------------------------------------------------
@dataclass
class LabProfile:
    """
    All cost, labor, and time parameters used by cloning operators.

    Cost units: USD.  Labor units: hours hands-on.  Time units: calendar days.
    """

    # Reagent catalog — copy on construction so per-instance overrides are safe.
    catalog: Dict[str, Reagent] = field(default_factory=lambda: dict(DEFAULT_CATALOG))

    # ------------------------------------------------------------------
    # Pricing helpers (instance methods so subclasses can override)
    # ------------------------------------------------------------------
    def primer_cost(self, length_bp: int) -> float:
        return primer_cost(length_bp)

    def synthesis_cost(self, length_bp: int) -> float:
        return synthesis_cost(length_bp)

    def synthesis_cost_per_bp(self, length_bp: int) -> float:
        return synthesis_cost_per_bp(length_bp)

    def reagent_lines(self, workflow: str) -> List[Tuple[Reagent, float]]:
        return reagent_lines(workflow, self.catalog)

    def workflow_reagent_cost_per_construct(self, workflow: str) -> float:
        return workflow_reagent_cost_per_construct(workflow, self.catalog)

    # ------------------------------------------------------------------
    # Sequencing — ONT (long-read) default, 1 read / construct
    # ------------------------------------------------------------------
    @property
    def sequencing_cost_usd(self) -> float:
        return self.catalog["ont_seq_read"].cost_per_rxn

    sequencing_reads_per_construct: int = 1

    # ------------------------------------------------------------------
    # Backwards-compat flat fields — derived from the catalog so any
    # legacy callers (operators that read self.lab.pcr_rxn_cost_usd
    # directly) keep working without changes.
    # ------------------------------------------------------------------
    @property
    def primer_cost_usd(self) -> float:
        # Avg 30 mer Gibson primer; operators with primer-length info
        # should call self.primer_cost(len) directly for accuracy.
        return self.primer_cost(30)

    @property
    def pcr_rxn_cost_usd(self) -> float:
        return self.catalog["neb_q5_rxn"].cost_per_rxn

    @property
    def gel_lane_cost_usd(self) -> float:
        # Per-lane gel cost = agarose + stain + ladder share, derived
        # from the 60-mL-gel recipe (10 lanes per gel).
        agarose = self.catalog["ultrapure_agarose_g"].cost_per_rxn * 0.6
        stain = self.catalog["sybr_safe_ul"].cost_per_rxn * 6
        ladder = self.catalog["ladder_lane"].cost_per_rxn
        return (agarose + stain + ladder) / 10  # 10 lanes / gel

    @property
    def hifi_assembly_cost_usd(self) -> float:
        return self.catalog["nebuilder_hifi"].cost_per_rxn

    @property
    def transformation_cost_usd(self) -> float:
        return self.catalog["neb_10beta"].cost_per_rxn + self.catalog["lb_agar_plate"].cost_per_rxn

    @property
    def miniprep_cost_usd(self) -> float:
        # ONT sequencing services include miniprep upstream, so the
        # explicit miniprep line is folded into the seq workflow.
        return 0.0

    # ------------------------------------------------------------------
    # Synthesis legacy flats (still referenced by gibson_operator
    # decision logic). Resolve to the tier endpoints.
    # ------------------------------------------------------------------
    @property
    def synthesis_cost_per_bp_simple(self) -> float:
        return SYNTHESIS_TIERS[0][1]

    @property
    def synthesis_cost_per_bp_complex(self) -> float:
        return SYNTHESIS_TIERS[-1][1]

    synthesis_threshold_bp: int = 2000
    synthesis_complex_gc_threshold: float = 0.70
    synthesis_lead_time_days: float = 7.0

    # ------------------------------------------------------------------
    # Labor (hands-on hours per step) — unchanged from prior version.
    # ------------------------------------------------------------------
    pcr_labor_hours: float = 0.5
    gel_labor_hours: float = 0.5
    assembly_labor_hours: float = 0.5
    transformation_labor_hours: float = 1.0
    colony_pcr_labor_hours: float = 1.0
    miniprep_labor_hours: float = 0.5
    sequencing_labor_hours: float = 0.25
    primer_ordering_labor_hours: float = 0.25

    # ------------------------------------------------------------------
    # Calendar time (days per step, including wait times)
    # ------------------------------------------------------------------
    pcr_days: float = 0.25
    gel_days: float = 0.25
    assembly_days: float = 0.1
    transformation_days: float = 0.5
    overnight_incubation_days: float = 0.5
    colony_selection_days: float = 0.5
    miniprep_days: float = 0.25
    sequencing_turnaround_days: float = 1.5

    # ------------------------------------------------------------------
    # Screening assumptions
    # ------------------------------------------------------------------
    colonies_to_screen: int = 6
    minipreps_per_construct: int = 2
    sequencing_reads_per_junction: int = 0   # ONT covers the whole plasmid
    sequencing_reads_overhead: int = 1


DEFAULT_LAB_PROFILE = LabProfile()
