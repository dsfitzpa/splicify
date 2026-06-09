"""
Restriction enzyme metadata and compatibility helpers.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Union


@dataclass(frozen=True)
class RestrictionEnzyme:
    name: str
    recognition_seq: str
    overhang_seq: str
    overhang_type: str  # 5prime | 3prime | blunt
    left_of_cut: int
    right_of_cut: int
    dam_sensitive: bool = False
    dcm_sensitive: bool = False
    buffer: str = "CutSmart"
    cost_per_rxn_usd: float = 1.2
    star_activity: bool = False


RE_DATABASE: Dict[str, RestrictionEnzyme] = {
    "EcoRI": RestrictionEnzyme("EcoRI", "GAATTC", "AATT", "5prime", 1, 5),
    "BamHI": RestrictionEnzyme("BamHI", "GGATCC", "GATC", "5prime", 1, 5, dam_sensitive=True),
    "HindIII": RestrictionEnzyme("HindIII", "AAGCTT", "AGCT", "5prime", 1, 5),
    "NotI": RestrictionEnzyme("NotI", "GCGGCCGC", "GGCC", "5prime", 2, 6),
    "XhoI": RestrictionEnzyme("XhoI", "CTCGAG", "TCGA", "5prime", 1, 5),
    "NheI": RestrictionEnzyme("NheI", "GCTAGC", "CTAG", "5prime", 1, 5, dam_sensitive=True),
    "XbaI": RestrictionEnzyme("XbaI", "TCTAGA", "CTAG", "5prime", 1, 5),
    "SpeI": RestrictionEnzyme("SpeI", "ACTAGT", "CTAG", "5prime", 1, 5),
    "SalI": RestrictionEnzyme("SalI", "GTCGAC", "TCGA", "5prime", 1, 5),
    "KpnI": RestrictionEnzyme("KpnI", "GGTACC", "GTAC", "3prime", 5, 1),
    "NcoI": RestrictionEnzyme("NcoI", "CCATGG", "CATG", "5prime", 1, 5, dcm_sensitive=True),
    "BglII": RestrictionEnzyme("BglII", "AGATCT", "GATC", "5prime", 1, 5),
    "AgeI": RestrictionEnzyme("AgeI", "ACCGGT", "CCGG", "5prime", 1, 5),
    "PstI": RestrictionEnzyme("PstI", "CTGCAG", "TGCA", "3prime", 5, 1),
    "MluI": RestrictionEnzyme("MluI", "ACGCGT", "CGCG", "5prime", 1, 5),
    "ClaI": RestrictionEnzyme("ClaI", "ATCGAT", "CGAT", "5prime", 2, 4, dam_sensitive=True),
    "PacI": RestrictionEnzyme("PacI", "TTAATTAA", "TAAT", "blunt", 4, 4),
    "AscI": RestrictionEnzyme("AscI", "GGCGCGCC", "CGCG", "5prime", 2, 6),
    "SacI": RestrictionEnzyme("SacI", "GAGCTC", "AGCT", "5prime", 1, 5),
    "AvrII": RestrictionEnzyme("AvrII", "CCTAGG", "CTAG", "5prime", 1, 5),
}


# --- NEB-COMMON-ENZYMES-EXPANSION-2026-04-19 ---
# Commonly-used NEB Type II enzymes added 2026-04-19 (pure ACGT only;
# degenerate-base enzymes handled by a future IUPAC-aware scanner pass).
RE_DATABASE.update({
    "EcoRV":  RestrictionEnzyme("EcoRV",  "GATATC",   "",     "blunt",  3, 3),
    "SmaI":   RestrictionEnzyme("SmaI",   "CCCGGG",   "",     "blunt",  3, 3),
    "XmaI":   RestrictionEnzyme("XmaI",   "CCCGGG",   "CCGG", "5prime", 1, 5),
    "StuI":   RestrictionEnzyme("StuI",   "AGGCCT",   "",     "blunt",  3, 3),
    "NdeI":   RestrictionEnzyme("NdeI",   "CATATG",   "TA",   "5prime", 2, 4),
    "SphI":   RestrictionEnzyme("SphI",   "GCATGC",   "CATG", "3prime", 5, 1),
    "ApaI":   RestrictionEnzyme("ApaI",   "GGGCCC",   "GGCC", "3prime", 5, 1),
    "SacII":  RestrictionEnzyme("SacII",  "CCGCGG",   "GC",   "3prime", 4, 2),
    "MfeI":   RestrictionEnzyme("MfeI",   "CAATTG",   "AATT", "5prime", 1, 5),
    "SwaI":   RestrictionEnzyme("SwaI",   "ATTTAAAT", "",     "blunt",  4, 4),
    "FseI":   RestrictionEnzyme("FseI",   "GGCCGGCC", "CCGG", "3prime", 6, 2),
    "PmeI":   RestrictionEnzyme("PmeI",   "GTTTAAAC", "",     "blunt",  4, 4),
    "BstBI":  RestrictionEnzyme("BstBI",  "TTCGAA",   "CG",   "5prime", 2, 4),
    "BsiWI":  RestrictionEnzyme("BsiWI",  "CGTACG",   "GTAC", "5prime", 1, 5),
    "NsiI":   RestrictionEnzyme("NsiI",   "ATGCAT",   "TGCA", "3prime", 5, 1),
    "AflII":  RestrictionEnzyme("AflII",  "CTTAAG",   "TTAA", "5prime", 1, 5),
    "SnaBI":  RestrictionEnzyme("SnaBI",  "TACGTA",   "",     "blunt",  3, 3),
    "AfeI":   RestrictionEnzyme("AfeI",   "AGCGCT",   "",     "blunt",  3, 3),
    "EagI":   RestrictionEnzyme("EagI",   "CGGCCG",   "GGCC", "5prime", 1, 5),
    "MscI":   RestrictionEnzyme("MscI",   "TGGCCA",   "",     "blunt",  3, 3),
    "NruI":   RestrictionEnzyme("NruI",   "TCGCGA",   "",     "blunt",  3, 3),
    "PvuI":   RestrictionEnzyme("PvuI",   "CGATCG",   "AT",   "3prime", 4, 2),
    "PvuII":  RestrictionEnzyme("PvuII",  "CAGCTG",   "",     "blunt",  3, 3),
    "ScaI":   RestrictionEnzyme("ScaI",   "AGTACT",   "",     "blunt",  3, 3),
    "HpaI":   RestrictionEnzyme("HpaI",   "GTTAAC",   "",     "blunt",  3, 3),
    "ZraI":   RestrictionEnzyme("ZraI",   "GACGTC",   "",     "blunt",  3, 3),
    "AatII":  RestrictionEnzyme("AatII",  "GACGTC",   "ACGT", "3prime", 5, 1),
    "BsrGI":  RestrictionEnzyme("BsrGI",  "TGTACA",   "GTAC", "5prime", 1, 5),
    "DraI":   RestrictionEnzyme("DraI",   "TTTAAA",   "",     "blunt",  3, 3),
    "FspI":   RestrictionEnzyme("FspI",   "TGCGCA",   "",     "blunt",  3, 3),
    "MseI":   RestrictionEnzyme("MseI",   "TTAA",     "TA",   "5prime", 1, 3),
    "PmlI":   RestrictionEnzyme("PmlI",   "CACGTG",   "",     "blunt",  3, 3),
    "SbfI":   RestrictionEnzyme("SbfI",   "CCTGCAGG", "TGCA", "3prime", 6, 2),
    "PciI":   RestrictionEnzyme("PciI",   "ACATGT",   "CATG", "5prime", 1, 5),
    "BssHII": RestrictionEnzyme("BssHII", "GCGCGC",   "CGCG", "5prime", 1, 5),
})
# Extend compatible overhang groups to include the newly added matching families.

# Known compatible overhang families used in hybrid ligations.
_COMPATIBLE_OVERHANG_GROUPS = {
    "CTAG": {"XbaI", "SpeI", "NheI", "AvrII"},
    "GTAC": {"KpnI", "BsrGI", "BsiWI"},
    "CATG": {"NcoI", "SphI", "PciI"},
    "CCGG": {"AgeI", "XmaI", "FseI", "EagI"},
    "AATT": {"EcoRI", "MfeI"},
    "AGCT": {"HindIII", "SacI"},
    "TGCA": {"PstI", "NsiI", "SbfI"},
    "ACGT": {"AatII"},
    "TA": {"NdeI", "MseI"},
    "CGCG": {"MluI", "BssHII", "AscI"},
    "AT": {"ClaI", "PvuI"},
    "TCGA": {"XhoI", "SalI"},
    "GATC": {"BamHI", "BglII"},
}

COMPATIBLE_ENZYME_PAIRS = {
    frozenset({a, b})
    for members in _COMPATIBLE_OVERHANG_GROUPS.values()
    for a in members
    for b in members
    if a != b
}


def get_enzyme(name: str) -> Optional[RestrictionEnzyme]:
    return RE_DATABASE.get(name)


def are_compatible(a: Union[str, RestrictionEnzyme], b: Union[str, RestrictionEnzyme]) -> bool:
    enz_a = RE_DATABASE.get(a, a) if isinstance(a, str) else a
    enz_b = RE_DATABASE.get(b, b) if isinstance(b, str) else b
    if not isinstance(enz_a, RestrictionEnzyme) or not isinstance(enz_b, RestrictionEnzyme):
        return False
    if enz_a.name == enz_b.name:
        return True
    if enz_a.overhang_type == "blunt" and enz_b.overhang_type == "blunt":
        return True
    if enz_a.overhang_type != enz_b.overhang_type:
        return False
    return enz_a.overhang_seq == enz_b.overhang_seq


def compute_scar(
    a: Union[str, RestrictionEnzyme],
    b: Union[str, RestrictionEnzyme],
) -> Optional[str]:
    """
    Return the ligation scar implied by two enzymes.

    For compatible sticky ends this is the shared overhang sequence.
    For blunt ends the scar is empty.
    If enzymes are incompatible, returns None.
    """
    enz_a = RE_DATABASE.get(a, a) if isinstance(a, str) else a
    enz_b = RE_DATABASE.get(b, b) if isinstance(b, str) else b
    if not isinstance(enz_a, RestrictionEnzyme) or not isinstance(enz_b, RestrictionEnzyme):
        return None

    if enz_a.overhang_type == "blunt" and enz_b.overhang_type == "blunt":
        return ""
    if not are_compatible(enz_a, enz_b):
        return None
    return enz_a.overhang_seq
