"""VirtualConstruct: parts arranged in order + orientation. Materializes
to a sequence that's piped through annotate_llm_cached for verification."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Optional

from agent_v2.builder.part_set import Part, _rc


# Type IIs enzymes the builder recognises for Golden Gate compatibility
# checks. Same list as the annotation pipeline uses for POL3-GG-01.
TYPE_IIS_SITES = {
    "BsmBI": "CGTCTC",
    "Esp3I": "CGTCTC",
    "BbsI":  "GAAGAC",
    "BsaI":  "GGTCTC",
    "SapI":  "GCTCTTC",
    "AarI":  "CACCTGC",
}

# Common Type II restriction enzymes used in restriction cloning.
COMMON_TYPE_II = {
    "EcoRI":  "GAATTC",
    "BamHI":  "GGATCC",
    "HindIII":"AAGCTT",
    "XhoI":   "CTCGAG",
    "NotI":   "GCGGCCGC",
    "SalI":   "GTCGAC",
    "PstI":   "CTGCAG",
    "SacI":   "GAGCTC",
    "NheI":   "GCTAGC",
    "SpeI":   "ACTAGT",
    "XbaI":   "TCTAGA",
    "KpnI":   "GGTACC",
    "EcoRV":  "GATATC",
    "PmeI":   "GTTTAAAC",
    "MluI":   "ACGCGT",
}


# Gateway att sites for Gateway-method detection (BP / LR reactions).
GATEWAY_SITES = {
    "attB1": "ACAAGTTTGTACAAAAAAGCAGGCT",
    "attB2": "ACCACTTTGTACAAGAAAGCTGGGT",
    "attL1": "CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCT",
    "attL2": "ACCCAGCTTTCTTGTACAAAGTTGGCATTATAAGAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTG",
    "attR1": "ACAAGTTTGTACAAAAAAGCTGAACG",
    "attR2": "ACCACTTTGTACAAGAAAGCTGGGTC",
}


# Minimum homology for Gibson assembly (NEB recommends 15-25 bp).
GIBSON_MIN_OVERLAP = 15


@dataclass
class Slot:
    """One occupied position in the virtual construct."""
    part: Part
    orientation: int = 1            # 1 (5'→3') or -1 (reverse complement)

    @property
    def body(self) -> str:
        return self.part.body(self.orientation)

    @property
    def upstream(self) -> str:
        # The junction visible at the 5' end of this slot after orientation.
        if self.orientation == 1:
            return self.part.upstream_junction
        return _rc(self.part.downstream_junction)

    @property
    def downstream(self) -> str:
        if self.orientation == 1:
            return self.part.downstream_junction
        return _rc(self.part.upstream_junction)


@dataclass
class JunctionAnalysis:
    """Compatibility findings for the seam between slot_i and slot_{i+1}."""
    left_idx: int
    right_idx: int
    gibson_overlap_bp: int = 0          # max homology between left.downstream and right.upstream
    type_iis_enzyme: Optional[str] = None
    type_iis_overhang: Optional[str] = None
    type_ii_enzyme: Optional[str] = None
    type_ii_site: Optional[str] = None
    gateway_pair: Optional[tuple[str, str]] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "left_idx": self.left_idx,
            "right_idx": self.right_idx,
            "gibson_overlap_bp": self.gibson_overlap_bp,
            "type_iis_enzyme": self.type_iis_enzyme,
            "type_iis_overhang": self.type_iis_overhang,
            "type_ii_enzyme": self.type_ii_enzyme,
            "type_ii_site": self.type_ii_site,
            "gateway_pair": list(self.gateway_pair) if self.gateway_pair else None,
        }


@dataclass
class VirtualConstruct:
    slots: list[Slot] = field(default_factory=list)
    topology: str = "circular"

    def materialize(self) -> str:
        """Concatenate body sequences in current order/orientation.
        Note: the junctions are NOT included in the final construct —
        they're only used for compatibility analysis. The body
        sequences themselves are the parts that get assembled."""
        return "".join(s.body for s in self.slots)

    def junctions(self) -> list[JunctionAnalysis]:
        """Pair-wise junction analysis for adjacent slots (wrapping
        around when topology is circular)."""
        n = len(self.slots)
        if n < 2:
            return []
        pairs = [(i, (i + 1) % n) for i in range(n)] if self.topology == "circular" \
                else [(i, i + 1) for i in range(n - 1)]
        out: list[JunctionAnalysis] = []
        for i, j in pairs:
            left = self.slots[i]
            right = self.slots[j]
            j_analysis = analyze_junction(left, right, idx_pair=(i, j))
            out.append(j_analysis)
        return out

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_slots": len(self.slots),
            "topology": self.topology,
            "slots": [{"name": s.part.name, "role": s.part.role,
                        "orientation": s.orientation, "length_bp": s.part.length_bp,
                        "part_id": s.part.part_id} for s in self.slots],
            "junctions": [j.to_dict() for j in self.junctions()],
            "total_bp": sum(s.part.length_bp for s in self.slots),
        }


def _longest_common_suffix_prefix(a: str, b: str, min_len: int = 1) -> int:
    """Length of the longest suffix of `a` that is also a prefix of `b`."""
    if not a or not b:
        return 0
    max_check = min(len(a), len(b), 50)
    for k in range(max_check, min_len - 1, -1):
        if a[-k:] == b[:k]:
            return k
    return 0


def analyze_junction(left: Slot, right: Slot, *, idx_pair: tuple[int, int]) -> JunctionAnalysis:
    """Inspect the seam between left.downstream and right.upstream and
    record what assembly methods could bridge it WITHOUT adding new
    sequence."""
    out = JunctionAnalysis(left_idx=idx_pair[0], right_idx=idx_pair[1])

    # Gibson: longest exact homology between left's 3' tail and right's
    # 5' head. The 50 bp junction sequences are designed for exactly
    # this comparison.
    out.gibson_overlap_bp = _longest_common_suffix_prefix(
        left.downstream, right.upstream, min_len=GIBSON_MIN_OVERLAP,
    )

    # Golden Gate: a Type IIs recognition site must sit in left's
    # downstream junction OR right's upstream junction (typically both,
    # cutting outward to leave matching 4 bp overhangs).
    for enz, site in TYPE_IIS_SITES.items():
        rc = _rc(site)
        if (site in left.downstream or rc in left.downstream) and \
           (site in right.upstream or rc in right.upstream):
            out.type_iis_enzyme = enz
            # Best-effort overhang extraction — Type IIs sites cut 1 nt
            # downstream with a 4 bp offset; we just record that the
            # enzyme matches, leaving overhang computation to the
            # downstream cloning workflow.
            out.type_iis_overhang = "<computed by cloning workflow>"
            break

    # Restriction (Type II): same site present in both halves of the
    # junction, palindromic so the same string finds itself.
    for enz, site in COMMON_TYPE_II.items():
        if site in left.downstream and site in right.upstream:
            out.type_ii_enzyme = enz
            out.type_ii_site = site
            break

    # Gateway: matching attB / attL / attR pair.
    for name_a, site_a in GATEWAY_SITES.items():
        if site_a in left.downstream:
            for name_b, site_b in GATEWAY_SITES.items():
                if site_b in right.upstream and name_a[-1] == name_b[-1]:
                    out.gateway_pair = (name_a, name_b)
                    break
            if out.gateway_pair:
                break

    return out


@dataclass
class MethodAssessment:
    """The set of assembly methods feasible for a given
    VirtualConstruct, ranked by preference."""
    feasible: list[str] = field(default_factory=list)
    rejected: dict[str, str] = field(default_factory=dict)        # method → reason
    junction_findings: list[JunctionAnalysis] = field(default_factory=list)
    pick: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "feasible": list(self.feasible),
            "rejected": dict(self.rejected),
            "pick": self.pick,
            "junctions": [j.to_dict() for j in self.junction_findings],
        }


def assess_methods(construct: VirtualConstruct) -> MethodAssessment:
    """Walk every junction and decide which assembly method(s) work
    without adding new sequence.

    Rules:
      - **gibson**: every junction has ≥GIBSON_MIN_OVERLAP bp Gibson homology.
      - **golden_gate**: every junction has a recognised Type IIs site
        on both sides (same enzyme).
      - **restriction**: every junction shares a Type II site on both
        sides (same enzyme); same enzyme used for every junction.
      - **gateway**: every junction has a recognised att pair.
      - **single-part**: 1 slot → site-directed mutagenesis (SDM).
      - **sgrna_gg**: 3-slot pattern (pol3_promoter / stuffer / scaffold)
        with Type IIs flanks — special case of golden_gate.
    """
    out = MethodAssessment()
    junctions = construct.junctions()
    out.junction_findings = junctions

    # Single-part SDM
    if len(construct.slots) == 1:
        out.feasible.append("sdm")

    if junctions:
        if all(j.gibson_overlap_bp >= GIBSON_MIN_OVERLAP for j in junctions):
            out.feasible.append("gibson")
        else:
            failing = [j for j in junctions if j.gibson_overlap_bp < GIBSON_MIN_OVERLAP]
            out.rejected["gibson"] = (
                f"{len(failing)}/{len(junctions)} junctions below {GIBSON_MIN_OVERLAP} bp homology"
            )

        gg_enzymes = [j.type_iis_enzyme for j in junctions if j.type_iis_enzyme]
        if len(gg_enzymes) == len(junctions) and len(set(gg_enzymes)) == 1:
            out.feasible.append("golden_gate")
        elif gg_enzymes:
            out.rejected["golden_gate"] = (
                f"only {len(gg_enzymes)}/{len(junctions)} junctions have a Type IIs site "
                f"or enzymes don't match across junctions"
            )

        rest_enzymes = [j.type_ii_enzyme for j in junctions if j.type_ii_enzyme]
        if len(rest_enzymes) == len(junctions) and len(set(rest_enzymes)) == 1:
            out.feasible.append("restriction")
        elif rest_enzymes:
            out.rejected["restriction"] = (
                f"only {len(rest_enzymes)}/{len(junctions)} junctions share a Type II site"
            )

        gw_pairs = [j.gateway_pair for j in junctions if j.gateway_pair]
        if len(gw_pairs) == len(junctions):
            out.feasible.append("gateway")
        elif gw_pairs:
            out.rejected["gateway"] = (
                f"only {len(gw_pairs)}/{len(junctions)} junctions have a Gateway pair"
            )

    # sgRNA Golden Gate special case: 3 slots, pol3_promoter / stuffer
    # / scaffold roles, GG feasible.
    roles_sequence = [s.part.role for s in construct.slots]
    if "golden_gate" in out.feasible and \
       any("pol3" in r for r in roles_sequence) and "stuffer" in roles_sequence \
       and any("scaffold" in r for r in roles_sequence):
        out.feasible.insert(0, "sgrna_gg")

    # Pick preference: gibson (most flexible) > golden_gate > restriction > gateway > sdm > sgrna_gg first if it matched
    preference = ["sgrna_gg", "gibson", "golden_gate", "restriction", "gateway", "sdm"]
    for m in preference:
        if m in out.feasible:
            out.pick = m
            break
    return out
