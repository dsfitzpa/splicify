"""
Gateway Cloning att site definitions based on GateWayMine consensus sequences.

This module provides:
- All standard Gateway att site sequences (attB, attP, attL, attR for sites 1-5)
- Detection functions with exact and fuzzy matching
- Orthogonality validation based on 7bp core sequences
- Site type parsing utilities

Gateway Biology:
- BP Reaction: attB + attP → attL + attR
  Example: attB1 + attP1 → attL1 + attR1
- LR Reaction: attL + attR → attB + attP
  Example: attL1 + attR1 → attB1 + attP1
- Orthogonality: The 7bp core sequence determines specificity
  - Site 1: GTACAAA
  - Site 2: GTACAAG
  - Site 3: GTATAAT
  - Site 4: GTATAGA
  - Site 5: GTATACA

Consensus sequences from GateWayMine analysis (github.com/manulera/GateWayMine)
Using simplified consensus (most common bases, IUPAC codes resolved):
- M (A or C) → A
- R (A or G) → A
- K (G or T) → T
- Y (C or T) → T
- W (A or T) → A
- S (C or G) → C
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import List, Optional, Tuple
from Bio.Seq import Seq


# Gateway att site consensus sequences (simplified from IUPAC ambiguity codes)
# Structure: LEFT_FLANK + CORE_7bp + RIGHT_FLANK

# Left flanking sequence (long, ~80bp) - found in attP and attL
LEFT_LONG = "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAGCAATGCTTTTTATAATGCCAACTT"

# Right flanking sequence (long, ~100bp) - found in attP and attR
RIGHT_LONG = "TGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACAATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"

# Left flanking sequence (short, ~7bp) - found in attB and attR
LEFT_SHORT = "CAAATTT"

# Right flanking sequence (short, ~7bp) - found in attB and attL
RIGHT_SHORT = "AAATTTG"

# Core 7bp sequences (determine specificity)
CORE_SEQUENCES = {
    "1": "GTACAAA",
    "2": "GTACAAG",
    "3": "GTATAAT",
    "4": "GTATAGA",
    "5": "GTATACA",
}

# Build att site sequences
# attP sites: LONG_LEFT + CORE + LONG_RIGHT
# attB sites: SHORT_LEFT + CORE + SHORT_RIGHT
# attL sites: LONG_LEFT + CORE + SHORT_RIGHT
# attR sites: SHORT_LEFT + CORE + LONG_RIGHT

GATEWAY_ATT_SITES = {}

# Invitrogen pDONR vectors use slightly different sequences than GateWayMine consensus
# Add actual Invitrogen sequences first (these are the most common)
INVITROGEN_LEFT_LONG = "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTT"
INVITROGEN_RIGHT_LONG = "GAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA"
INVITROGEN_JUNCTION = "AAAGCTG"  # Different from GateWayMine "AAATTTG"

# Build Invitrogen-specific att sites (most common in actual use)
# Build Gateway att site sequences
# Use EXACT sequences from Invitrogen Gateway manual for attB sites (25 bp)
# attB1 and attB2 must be exactly as specified in the manual

# attB sites (25 bp total - from Gateway Technology manual page 24)
GATEWAY_ATT_SITES["attB1"] = "acaagtttgtacaaaaaagcaggct"  # 25 bp - Correct sequence from Gateway manual
GATEWAY_ATT_SITES["attB2"] = "accactttgtacaagaaagctgggt"  # 25 bp - Correct sequence from Gateway manual

# --- MULTISITE-GATEWAY-B3-B5-EXPANSION-2026-04-19 ---
# MultiSite Gateway Pro attB3/4/5 - 21 bp each (Invitrogen MultiSite Gateway manual).
# These are the canonical attB sequences downstream labs use for 4-fragment assembly.
GATEWAY_ATT_SITES["attB3"] = "caactttgtataataaagttg"  # 21 bp - MultiSite Gateway manual
GATEWAY_ATT_SITES["attB4"] = "caactttgtatagaaaagttg"  # 21 bp - MultiSite Gateway manual
GATEWAY_ATT_SITES["attB5"] = "caactttgtatacaaaagttg"  # 21 bp - MultiSite Gateway manual


# Build other att sites using flanking sequences
for num, core in CORE_SEQUENCES.items():
    # Invitrogen attP: longer left flank + core + AAAGCTG junction + shorter right flank
    GATEWAY_ATT_SITES[f"attP{num}"] = INVITROGEN_LEFT_LONG + core + INVITROGEN_JUNCTION + INVITROGEN_RIGHT_LONG
    # attL from BP reaction: left long + core + short right
    GATEWAY_ATT_SITES[f"attL{num}"] = INVITROGEN_LEFT_LONG + core + RIGHT_SHORT
    # attR from BP reaction: short left + core + junction + long right
    GATEWAY_ATT_SITES[f"attR{num}"] = LEFT_SHORT + core + INVITROGEN_JUNCTION + INVITROGEN_RIGHT_LONG



# pDEST-style attR variants — attB-left + attP-right-truncated structure
# used in Invitrogen destination vectors (pLenti6/V5-DEST, pcDNA3.1-DEST, etc.)
# Canonical length ~126 bp, starting with the attB1/attB2 first-18-bp left arm.
GATEWAY_ATT_SITES["attR1_pDEST"] = (
    # First 18 bp of attB1 (left arm through core)
    "ACAAGTTTGTACAAAAAA"
    # Junction (4 bp — truncated AAAGCTG)
    "GCTG"
    # Truncated INVITROGEN_RIGHT_LONG (first 104 bp)
    "AACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGG"
)
GATEWAY_ATT_SITES["attR2_pDEST"] = (
    # First 18 bp of attB2 left arm (site-2 core GTACAAG)
    "ACCACTTTGTACAAGAAAG"
    # Junction (3 bp remainder)
    "CTG"
    # Truncated INVITROGEN_RIGHT_LONG (first 104 bp)
    "AACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGG"
)

# Add actual pDONR223/pDONR221 sequences (variant from Invitrogen catalog plasmids)
# These differ from the consensus in the recombination regions
PDONR_ATT_SITES = {
    "attP1_pDONR": "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA",
    "attP2_pDONR": "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTCTTATAATGCCAACTTTGTACAAGAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA",
}

# Add pDONR variants to the main dictionary
GATEWAY_ATT_SITES["attP1_pDONR"] = PDONR_ATT_SITES["attP1_pDONR"]
GATEWAY_ATT_SITES["attP2_pDONR"] = PDONR_ATT_SITES["attP2_pDONR"]


# Alternative representations using full consensus with IUPAC codes (for reference)
# These can be used for more permissive matching
GATEWAY_ATT_CONSENSUS = {
    "attP1": "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTGTACAAAAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATRCTGTAAAACACAACATATSCAGTCAYWWTGAATCAACTACTTAGATGGTATTAGTGACCTGTA",
    "attP2": "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTGTACAAGAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATRCTGTAAAACACAACATATSCAGTCAYWWTGAATCAACTACTTAGATGGTATTAGTGACCTGTA",
    "attB1": "CMASTWTGTACAAAAAAGYWG",
    "attB2": "CMASTWTGTACAAGAAAGYWG",
    "attL1": "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTGTACAAAAAAGYWG",
    "attL2": "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTTTGTACAAGAAAGYWG",
    "attR1": "CMASTTTGTACAAAAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATRCTGTAAAACACAACATATSCAGTCAYWWTGAATCAACTACTTAGATGGTATTAGTGACCTGTA",
    "attR2": "CMASTTTGTACAAGAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATRCTGTAAAACACAACATATSCAGTCAYWWTGAATCAACTACTTAGATGGTATTAGTGACCTGTA",
}


# ccdB gene sequence (starts with ATG, ~300 bp)
# Used for selection in Gateway destination vectors
CCDB_GENE = "ATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATA"


@dataclass
class AttSiteMatch:
    """Detected Gateway att site in a sequence."""
    site_type: str  # "attB1", "attP2", "attL1", etc.
    start: int
    end: int
    sequence: str
    core_sequence: str  # The 7bp core (GTACAAA, GTACAAG, etc.)
    match_quality: str  # "exact" | "fuzzy_1_mismatch" | "fuzzy_2_mismatches"
    strand: int  # 1 (forward) or -1 (reverse)

    def __str__(self) -> str:
        strand_str = "+" if self.strand == 1 else "-"
        return (
            f"{self.site_type} ({strand_str}) at {self.start}-{self.end} "
            f"[{self.match_quality}, core: {self.core_sequence}]"
        )


# 2026-05-06: orientation-aware Gateway substrate classifier — shared
# between the cloning operator, predesign router, and auto-router so the
# three see the same excision/inversion/intermolecular verdict for the
# same att pair. Mirrors the geometry rules in
# rule_based_module_detector._detect_gateway_orientation_modules.

_GATEWAY_COMPAT = {("B", "P"), ("P", "B"), ("L", "R"), ("R", "L")}


def classify_gateway_substrate(att_sites):
    """Pair compatible att sites (attB↔attP and attL↔attR, matching numeric
    suffix) and classify each pair's geometry.

    Args:
        att_sites: list of AttSiteMatch (or any object with .site_type,
                   .start, .end, .strand attributes).

    Returns:
        list[dict] — one entry per compatible pair detected, with keys:
          - kind: 'intermolecular' | 'excision' | 'inversion'
          - left_site, right_site: the two AttSiteMatch objects
          - pair_label: e.g. 'attL1↔attR1'
          - cargo_start, cargo_end: span between the inner edges of the pair
          - notes: human-readable explanation of the verdict
    """
    if not att_sites or len(att_sites) < 2:
        return []
    parsed = []
    for s in att_sites:
        kind, num = parse_att_site_type(s.site_type)
        if kind in ("?", ""):
            continue
        parsed.append((s, kind.upper(), num or "0"))
    parsed.sort(key=lambda t: t[0].start)

    out = []
    used = set()
    for i, (sa, ka, na) in enumerate(parsed):
        if id(sa) in used:
            continue
        for sb, kb, nb in parsed[i + 1:]:
            if id(sb) in used:
                continue
            if (ka, kb) not in _GATEWAY_COMPAT:
                continue
            if na != nb:
                continue
            ls, rs = sa.strand, sb.strand
            if ls == rs:
                kind = "excision"
                notes = ("Compatible att pair on the same strand → BP/LR "
                         "clonase excises the intervening cargo "
                         "(intramolecular deletion).")
            elif ls > 0 and rs < 0:
                kind = "intermolecular"
                notes = ("Compatible att pair on opposite strands pointing "
                         "INWARD → standard substrate for intermolecular "
                         "BP/LR recombination with a compatible vector.")
            else:
                kind = "inversion"
                notes = ("Compatible att pair on opposite strands pointing "
                         "OUTWARD → BP/LR clonase inverts the intervening "
                         "cargo (intramolecular inversion).")
            cargo_start = min(sa.end, sb.end)
            cargo_end = max(sa.start, sb.start)
            out.append({
                "kind": kind,
                "left_site": sa,
                "right_site": sb,
                "pair_label": f"att{ka}{na}↔att{kb}{nb}",
                "cargo_start": cargo_start,
                "cargo_end": cargo_end,
                "notes": notes,
            })
            used.add(id(sa))
            used.add(id(sb))
            break
    return out


def extract_core_sequence(att_site_type: str) -> str:
    """
    Get the 7bp core sequence for an att site type.

    Args:
        att_site_type: e.g., "attB1", "attP2", "attL3"

    Returns:
        7bp core sequence (e.g., "GTACAAA" for site 1)
    """
    match = re.match(r"att[BPLR]([12345])", att_site_type)
    if not match:
        return "???????"

    site_number = match.group(1)
    return CORE_SEQUENCES.get(site_number, "???????")


def parse_att_site_type(site_name: str) -> Tuple[str, str]:
    """
    Parse att site name into type and number.

    Args:
        site_name: e.g., "attB1", "attP2", "attL3"

    Returns:
        Tuple of (site_type, site_number), e.g., ("B", "1")
    """
    match = re.match(r"att([BPLR])([12345])(?:_.*)?$", site_name)
    if not match:
        return ("?", "?")

    return (match.group(1), match.group(2))


def validate_orthogonality(site1: AttSiteMatch, site2: AttSiteMatch) -> bool:
    """
    Validate that two att sites are orthogonal (compatible for recombination).

    Rules:
    1. attB can only recombine with attP (BP reaction)
    2. attL can only recombine with attR (LR reaction)
    3. Core sequences (7bp) must match
    4. Site numbers must match (attB1 with attP1, etc.)

    Args:
        site1: First att site
        site2: Second att site

    Returns:
        True if sites are compatible for recombination
    """
    type1, num1 = parse_att_site_type(site1.site_type)
    type2, num2 = parse_att_site_type(site2.site_type)

    # Check number matching
    if num1 != num2:
        return False

    # Check core sequence matching
    if site1.core_sequence != site2.core_sequence:
        return False

    # Check type compatibility
    # BP reaction: B + P → L + R
    if (type1 == "B" and type2 == "P") or (type1 == "P" and type2 == "B"):
        return True

    # LR reaction: L + R → B + P
    if (type1 == "L" and type2 == "R") or (type1 == "R" and type2 == "L"):
        return True

    return False


def get_recombination_products(site1_type: str, site2_type: str) -> Tuple[str, str]:
    """
    Get the product att site types from a recombination reaction.

    Args:
        site1_type: First att site type (e.g., "attB1")
        site2_type: Second att site type (e.g., "attP1")

    Returns:
        Tuple of (product1_type, product2_type), e.g., ("attL1", "attR1")
        Returns ("?", "?") if sites are incompatible
    """
    type1, num1 = parse_att_site_type(site1_type)
    type2, num2 = parse_att_site_type(site2_type)

    if num1 != num2:
        return ("?", "?")

    # BP reaction: attB + attP → attL + attR
    if set([type1, type2]) == {"B", "P"}:
        return (f"attL{num1}", f"attR{num1}")

    # LR reaction: attL + attR → attB + attP
    if set([type1, type2]) == {"L", "R"}:
        return (f"attB{num1}", f"attP{num1}")

    return ("?", "?")


def _fuzzy_match(pattern: str, sequence: str, max_mismatches: int) -> List[int]:
    """
    Find fuzzy matches of pattern in sequence allowing up to max_mismatches.

    Args:
        pattern: Pattern to search for
        sequence: Sequence to search in
        max_mismatches: Maximum number of mismatches allowed

    Returns:
        List of start positions where fuzzy matches were found
    """
    matches = []
    pattern_len = len(pattern)

    for i in range(len(sequence) - pattern_len + 1):
        substring = sequence[i:i + pattern_len]
        mismatches = sum(1 for a, b in zip(pattern, substring) if a != b)

        if mismatches <= max_mismatches:
            matches.append(i)

    return matches


def scan_att_sites(
    sequence: str,
    fuzzy_threshold: int = 2,
    search_attB_only: bool = False
) -> List[AttSiteMatch]:
    """
    Scan a DNA sequence for Gateway att sites.

    Performs:
    1. Exact matching for all att sites (forward and reverse)
    2. Fuzzy matching for attB sites (allow 1-2 mismatches for degraded sites)

    Args:
        sequence: DNA sequence to scan
        fuzzy_threshold: Maximum mismatches for fuzzy matching (0-2)
        search_attB_only: If True, only search for attB sites (faster)

    Returns:
        List of AttSiteMatch objects, sorted by start position
    """
    matches = []
    sequence_upper = sequence.upper()

    # Prepare reverse complement
    try:
        seq_obj = Seq(sequence_upper)
        rc_sequence = str(seq_obj.reverse_complement())
    except Exception:
        rc_sequence = ""

    # Sites to search
    sites_to_search = {}
    if search_attB_only:
        sites_to_search = {k: v for k, v in GATEWAY_ATT_SITES.items() if k.startswith("attB")}
    else:
        sites_to_search = GATEWAY_ATT_SITES

    # Exact matching
    for site_name, site_seq in sites_to_search.items():
        site_seq_upper = site_seq.upper()
        core = extract_core_sequence(site_name)

        # Forward strand
        start = 0
        while True:
            pos = sequence_upper.find(site_seq_upper, start)
            if pos == -1:
                break

            matches.append(AttSiteMatch(
                site_type=site_name,
                start=pos,
                end=pos + len(site_seq_upper),
                sequence=site_seq_upper,
                core_sequence=core,
                match_quality="exact",
                strand=1
            ))
            start = pos + 1

        # Reverse strand
        if rc_sequence:
            start = 0
            while True:
                pos = rc_sequence.find(site_seq_upper, start)
                if pos == -1:
                    break

                # Convert position to forward strand coordinates
                fwd_start = len(sequence) - (pos + len(site_seq_upper))
                fwd_end = len(sequence) - pos

                matches.append(AttSiteMatch(
                    site_type=site_name,
                    start=fwd_start,
                    end=fwd_end,
                    sequence=site_seq_upper,
                    core_sequence=core,
                    match_quality="exact",
                    strand=-1
                ))
                start = pos + 1

    # Fuzzy matching for attB sites only (they're short enough to fuzzy match)
    if fuzzy_threshold > 0:
        attB_sites = {k: v for k, v in GATEWAY_ATT_SITES.items() if k.startswith("attB")}

        for site_name, site_seq in attB_sites.items():
            site_seq_upper = site_seq.upper()
            core = extract_core_sequence(site_name)

            # Forward strand fuzzy matching
            for max_mm in range(1, fuzzy_threshold + 1):
                fuzzy_positions = _fuzzy_match(site_seq_upper, sequence_upper, max_mm)

                for pos in fuzzy_positions:
                    # Check if we already have an exact match here
                    already_matched = any(
                        m.start == pos and m.strand == 1 and m.site_type == site_name
                        for m in matches
                    )

                    if not already_matched:
                        found_seq = sequence_upper[pos:pos + len(site_seq_upper)]

                        quality = f"fuzzy_{max_mm}_mismatch" if max_mm == 1 else f"fuzzy_{max_mm}_mismatches"

                        matches.append(AttSiteMatch(
                            site_type=site_name,
                            start=pos,
                            end=pos + len(site_seq_upper),
                            sequence=found_seq,
                            core_sequence=core,
                            match_quality=quality,
                            strand=1
                        ))

            # Reverse strand fuzzy matching
            if rc_sequence:
                for max_mm in range(1, fuzzy_threshold + 1):
                    fuzzy_positions = _fuzzy_match(site_seq_upper, rc_sequence, max_mm)

                    for pos in fuzzy_positions:
                        fwd_start = len(sequence) - (pos + len(site_seq_upper))
                        fwd_end = len(sequence) - pos

                        # Check if already matched
                        already_matched = any(
                            m.start == fwd_start and m.strand == -1 and m.site_type == site_name
                            for m in matches
                        )

                        if not already_matched:
                            found_seq = rc_sequence[pos:pos + len(site_seq_upper)]

                            quality = f"fuzzy_{max_mm}_mismatch" if max_mm == 1 else f"fuzzy_{max_mm}_mismatches"

                            matches.append(AttSiteMatch(
                                site_type=site_name,
                                start=fwd_start,
                                end=fwd_end,
                                sequence=found_seq,
                                core_sequence=core,
                                match_quality=quality,
                                strand=-1
                            ))

    # Sort by start position
    matches.sort(key=lambda m: m.start)

    return matches


def scan_for_ccdb(sequence: str) -> List[Tuple[int, int]]:
    """
    Scan sequence for ccdB gene (used for negative selection in Gateway).

    Args:
        sequence: DNA sequence to scan

    Returns:
        List of (start, end) positions where ccdB was found
    """
    matches = []
    sequence_upper = sequence.upper()

    # Search for ccdB gene
    start = 0
    while True:
        pos = sequence_upper.find(CCDB_GENE.upper(), start)
        if pos == -1:
            break

        matches.append((pos, pos + len(CCDB_GENE)))
        start = pos + 1

    return matches
