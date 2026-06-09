"""
sgRNA Golden Gate Oligo Designer

Core logic for designing annealed oligos for cloning sgRNAs into Type IIS
enzyme-based vectors like lentiCRISPR v2.

Type IIS enzymes (BsmBI, BbsI, BsaI) cut outside their recognition sequence,
leaving 4-bp sticky overhangs. For sgRNA cloning:
1. Digest vector with Type IIS enzyme to remove filler sequence
2. Design two complementary oligos that anneal to form ds-DNA with matching overhangs
3. Ligate annealed oligos into the digested vector
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ..utils import normalize_dna, reverse_complement

# Try to import primer3 for Tm calculation
try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False


# ---------------------------------------------------------------------------
# Type IIS Enzyme Database
# ---------------------------------------------------------------------------

# Type IIS enzymes cut outside their recognition sequence
# Format: {name: {recognition: forward_seq, rc: reverse_complement, cut_offset: (top, bottom)}}
# cut_offset is the number of bp downstream from the end of the recognition site

TYPE_IIS_ENZYMES: Dict[str, Dict[str, str]] = {
    "BsmBI": {
        "recognition": "CGTCTC",  # Cuts 1 nt downstream (top), 5 nt downstream (bottom)
        "rc": "GAGACG",
    },
    "BbsI": {
        "recognition": "GAAGAC",  # Same cut pattern
        "rc": "GTCTTC",
    },
    "BsaI": {
        "recognition": "GGTCTC",  # Same cut pattern
        "rc": "GAGACC",
    },
}


@dataclass
class TypeIISSite:
    """Represents a Type IIS restriction site found in a sequence."""
    enzyme: str
    position: int  # Start position of recognition site (0-based)
    orientation: str  # "forward" or "reverse"
    recognition_seq: str
    cut_position: int  # Where the top strand is cut (0-based)


@dataclass
class SgRNAOligoDesign:
    """Result of sgRNA oligo design for Golden Gate cloning."""
    forward_oligo: str  # 5'-[5p_overhang][effective_grna]-3'
    reverse_oligo: str  # 5'-[RC(3p_overhang)][RC(effective_grna)]-3'
    forward_tm: float
    reverse_tm: float
    grna_sequence: str  # input gRNA as supplied
    effective_grna: str  # gRNA actually placed in the cassette (may differ for Pol III)
    grna_was_modified: bool  # True if a G was prepended and the last base dropped
    grna_modification_note: str  # human-readable description of the modification
    enzyme: str
    five_prime_overhang: str  # vector top-strand 4-bp window at the upstream cut
    three_prime_overhang: str  # vector top-strand 4-bp window at the downstream cut
    insert_length: int
    annealed_product_display: str
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Core Functions
# ---------------------------------------------------------------------------

def calculate_tm(sequence: str, salt_conc: float = 50.0, dna_conc: float = 250.0) -> float:
    """Calculate Tm using primer3 if available, otherwise use basic formula."""
    seq = normalize_dna(sequence)
    if not seq or len(seq) < 2:
        return 0.0

    if HAS_PRIMER3:
        try:
            tm = primer3.calcTm(seq, mv_conc=salt_conc, dv_conc=0.0, dntp_conc=0.0, dna_conc=dna_conc)
            return float(tm)
        except Exception:
            pass

    # Fallback: basic nearest-neighbor approximation
    gc_count = seq.count("G") + seq.count("C")
    at_count = seq.count("A") + seq.count("T")
    if len(seq) < 14:
        return 2.0 * at_count + 4.0 * gc_count
    else:
        return 64.9 + 41.0 * (gc_count - 16.4) / len(seq)


def find_type_iis_sites(
    sequence: str,
    enzyme: str = "auto"
) -> List[TypeIISSite]:
    """
    Find Type IIS enzyme sites in a sequence.

    Args:
        sequence: DNA sequence to search
        enzyme: Enzyme name ("BsmBI", "BbsI", "BsaI") or "auto" to detect

    Returns:
        List of TypeIISSite objects sorted by position
    """
    seq_upper = normalize_dna(sequence)
    sites: List[TypeIISSite] = []

    enzymes_to_check = TYPE_IIS_ENZYMES.keys() if enzyme == "auto" else [enzyme]

    for enz_name in enzymes_to_check:
        if enz_name not in TYPE_IIS_ENZYMES:
            continue

        enz = TYPE_IIS_ENZYMES[enz_name]
        fwd_pattern = enz["recognition"]
        rev_pattern = enz["rc"]

        # Forward orientation: enzyme cuts downstream
        for m in re.finditer(fwd_pattern, seq_upper):
            sites.append(TypeIISSite(
                enzyme=enz_name,
                position=m.start(),
                orientation="forward",
                recognition_seq=fwd_pattern,
                cut_position=m.end() + 1,  # 1 nt downstream of recognition site end
            ))

        # Reverse orientation: enzyme cuts upstream (on the complement strand)
        for m in re.finditer(rev_pattern, seq_upper):
            sites.append(TypeIISSite(
                enzyme=enz_name,
                position=m.start(),
                orientation="reverse",
                recognition_seq=rev_pattern,
                cut_position=m.start() - 1,  # Cuts before the recognition site
            ))

    return sorted(sites, key=lambda x: x.position)


def find_cloning_site_pair(
    sites: List[TypeIISSite],
    enzyme: str
) -> Optional[Tuple[TypeIISSite, TypeIISSite]]:
    """
    Find a pair of Type IIS sites suitable for cloning.

    For standard sgRNA cloning cassettes:
    - First site: reverse orientation (GAGACG for BsmBI) - creates 5' overhang
    - Second site: forward orientation (CGTCTC for BsmBI) - creates 3' overhang

    Returns:
        Tuple of (upstream_site, downstream_site) or None if no valid pair found
    """
    enzyme_sites = [s for s in sites if s.enzyme == enzyme]

    # Look for reverse then forward pattern
    reverse_sites = [s for s in enzyme_sites if s.orientation == "reverse"]
    forward_sites = [s for s in enzyme_sites if s.orientation == "forward"]

    if not reverse_sites or not forward_sites:
        return None

    # Find pairs where reverse is upstream of forward
    for rev_site in reverse_sites:
        for fwd_site in forward_sites:
            if rev_site.position < fwd_site.position:
                return (rev_site, fwd_site)

    return None


def calculate_sticky_ends(
    sequence: str,
    upstream_site: TypeIISSite,
    downstream_site: TypeIISSite
) -> Tuple[str, str]:
    """
    Calculate the 4-bp sticky ends created by Type IIS enzyme digestion.

    For BsmBI cloning in lentiCRISPR v2:
    - 5' overhang from vector: typically CACC (after digestion)
    - 3' overhang from vector: typically GTTT (after digestion)

    Returns:
        Tuple of (five_prime_overhang, three_prime_overhang)
    """
    seq_upper = normalize_dna(sequence)

    # For reverse-oriented site (e.g., GAGACG):
    # The overhang is generated upstream of the recognition site
    # For BsmBI: cuts 5 nt before the recognition on bottom, 1 nt before on top
    # The 4-bp overhang is just before the recognition site
    five_prime_start = upstream_site.position - 4
    five_prime_overhang = seq_upper[five_prime_start:five_prime_start + 4] if five_prime_start >= 0 else "CACC"

    # For forward-oriented site (e.g., CGTCTC):
    # The overhang is generated downstream of the recognition site
    # The 4-bp overhang starts 1 nt after the recognition site end
    three_prime_start = downstream_site.position + len(downstream_site.recognition_seq) + 1
    three_prime_overhang = seq_upper[three_prime_start:three_prime_start + 4] if three_prime_start + 4 <= len(seq_upper) else "GTTT"

    return five_prime_overhang, three_prime_overhang


def adjust_grna_for_pol3(
    grna: str,
    prepend_g_for_pol3: bool,
    five_prime_overhang: str = "",
) -> Tuple[str, bool, str]:
    """For Pol III (U6 / H1 / 7SK) cassettes, transcription must initiate at
    a G. The user-supplied gRNA is preserved verbatim — we never trim its
    last base; we only PREPEND a G when neither the vector overhang nor the
    gRNA itself already supplies one.

    Rule (post 2026-05-07):
      - if `five_prime_overhang.endswith("G")` → vector supplies the initiating
        G; do NOT prepend, do NOT trim.
      - elif `grna.startswith("G")` → gRNA's own first base is the initiator;
        do NOT prepend, do NOT trim.
      - else → prepend a G; do NOT trim.

    Consequence: the transcribed spacer can be 20-22 nt depending on how many
    leading Gs the vector contributes (0, 1, or 2 typical) plus any prepend.
    Callers are expected to surface that length in the reply.

    Returns (effective_grna, was_modified, human_readable_note).
    """
    if not prepend_g_for_pol3:
        return grna, False, ""

    overhang_supplies_g = bool(five_prime_overhang) and five_prime_overhang.endswith("G")
    if overhang_supplies_g or grna.startswith("G"):
        return grna, False, ""

    effective = "G" + grna
    note = (
        f"gRNA did not start with G and vector sticky end "
        f"'{five_prime_overhang or '(none)'}' does not end with G; prepended G "
        f"so Pol III transcription initiates at G "
        f"({grna} → {effective}, {len(grna)} → {len(effective)} nt)."
    )
    return effective, True, note


def design_sgrna_oligos(
    grna_sequence: str,
    five_prime_overhang: str,   # vector top-strand 4-bp window at the upstream cut
    three_prime_overhang: str,  # vector top-strand 4-bp window at the downstream cut
    enzyme: str = "BsmBI",
    *,
    prepend_g_for_pol3: bool = True,
) -> SgRNAOligoDesign:
    """Design forward + reverse oligos for sgRNA Golden Gate cloning.

    The overhang convention is the **vector top-strand 4-bp window** at each
    cut (i.e. the 4 bases that re-form on the top strand after ligation),
    not the protruding ssDNA. With this convention:
      - forward oligo 5'-extension = `five_prime_overhang` AS-IS (the insert
        top strand restores those bases at the upstream cut).
      - reverse oligo 5'-extension = `RC(three_prime_overhang)` (the insert
        bottom strand carries the RC of the downstream top window).

    Pol III paths (`prepend_g_for_pol3=True`): if the gRNA does not start
    with G, prepend G and drop the last base — the trailing C on the reverse
    oligo is the complement of the prepended G.
    """
    grna = normalize_dna(grna_sequence)
    warnings: List[str] = []

    if len(grna) < 17:
        warnings.append(f"gRNA is short ({len(grna)} bp). Standard gRNAs are 17-25 bp.")
    if len(grna) > 25:
        warnings.append(f"gRNA is long ({len(grna)} bp). Standard gRNAs are 17-25 bp.")
    invalid_bases = set(grna) - set("ACGT")
    if invalid_bases:
        warnings.append(f"gRNA contains non-standard bases: {invalid_bases}")

    effective_grna, grna_modified, modification_note = adjust_grna_for_pol3(
        grna, prepend_g_for_pol3, five_prime_overhang=five_prime_overhang,
    )

    fwd_overhang = five_prime_overhang
    rev_overhang = reverse_complement(three_prime_overhang)

    forward_oligo = f"{fwd_overhang}{effective_grna}"
    reverse_oligo = f"{rev_overhang}{reverse_complement(effective_grna)}"

    forward_tm = calculate_tm(forward_oligo)
    reverse_tm = calculate_tm(reverse_oligo)

    if len(forward_oligo) > 60:
        warnings.append(f"Forward oligo is {len(forward_oligo)} bp - may be too long for standard synthesis.")
    if len(reverse_oligo) > 60:
        warnings.append(f"Reverse oligo is {len(reverse_oligo)} bp - may be too long for standard synthesis.")

    annealed_display = generate_annealed_display(
        five_prime_overhang=five_prime_overhang,
        three_prime_overhang=three_prime_overhang,
        forward_oligo=forward_oligo,
        reverse_oligo=reverse_oligo,
    )

    return SgRNAOligoDesign(
        forward_oligo=forward_oligo,
        reverse_oligo=reverse_oligo,
        forward_tm=round(forward_tm, 1),
        reverse_tm=round(reverse_tm, 1),
        grna_sequence=grna,
        effective_grna=effective_grna,
        grna_was_modified=grna_modified,
        grna_modification_note=modification_note,
        enzyme=enzyme,
        five_prime_overhang=five_prime_overhang,
        three_prime_overhang=three_prime_overhang,
        insert_length=len(effective_grna),
        annealed_product_display=annealed_display,
        warnings=warnings,
    )


def generate_annealed_display(
    *,
    five_prime_overhang: str,
    three_prime_overhang: str,
    forward_oligo: str,
    reverse_oligo: str,
) -> str:
    """Render the ligated insert region in monospace, showing both strands
    with their respective vector sticky ends.

    Layout (lentiCRISPR v2: top-strand windows CACC / GTTT, EMX1 gRNA):

        5'-----CACCGAGTCCGAGCAGAAGAAGAA----GTTT-3'   forward oligo + 3' vector sticky end
        3'-GTGG----CTCAGGCTCGTCTTCTTCTTCAAA-----5'   reverse oligo + 5' vector sticky end

    The TOP strand carries the forward oligo (which restores the upstream
    sticky end CACC) and ends with the vector's top-strand 3' sticky end
    (GTTT). The BOTTOM strand starts with the vector's bottom-strand 5'
    sticky end (GTGG = per-base complement of the top window CACC, drawn
    3'→5' left-to-right) and carries the reverse oligo drawn 3'→5'. Dashes
    mark recessed regions where the opposite strand's overhang protrudes.
    """
    _COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
    bottom_5p_sticky = five_prime_overhang.translate(_COMPLEMENT)
    top_3p_sticky = three_prime_overhang

    overhang_left = len(five_prime_overhang)
    overhang_right = len(three_prime_overhang)
    dash_left = "-" * overhang_left
    dash_right = "-" * overhang_right

    bottom_oligo_3to5 = reverse_oligo[::-1]

    top_line = (
        f"  5'-{dash_left}{forward_oligo}{dash_right}{top_3p_sticky}-3'"
        f"   forward oligo + 3' vector sticky end"
    )
    bot_line = (
        f"  3'-{bottom_5p_sticky}{dash_left}{bottom_oligo_3to5}{dash_right}-5'"
        f"   reverse oligo + 5' vector sticky end"
    )

    return "\n".join([top_line, bot_line])


def design_sgrna_oligos_from_vector(
    vector_sequence: str,
    grna_sequence: str,
    enzyme: str = "auto"
) -> SgRNAOligoDesign:
    """
    Design sgRNA oligos by analyzing the vector sequence to determine overhangs.

    Args:
        vector_sequence: Full vector sequence (GenBank format or raw DNA)
        grna_sequence: The gRNA target sequence
        enzyme: Type IIS enzyme ("BsmBI", "BbsI", "BsaI", or "auto")

    Returns:
        SgRNAOligoDesign with oligo sequences
    """
    # Parse GenBank if needed
    if "ORIGIN" in vector_sequence:
        # Extract sequence from GenBank format
        origin_idx = vector_sequence.find("ORIGIN")
        seq_section = vector_sequence[origin_idx:]
        # Remove ORIGIN header and // footer, extract just bases
        lines = seq_section.split("\n")[1:]  # Skip ORIGIN line
        seq_parts = []
        for line in lines:
            if line.strip() == "//":
                break
            # Remove line numbers and spaces
            seq_parts.append(re.sub(r"[\d\s]", "", line))
        vector_seq = "".join(seq_parts).upper()
    else:
        vector_seq = normalize_dna(vector_sequence)

    # Find Type IIS sites
    sites = find_type_iis_sites(vector_seq, enzyme)

    if not sites:
        # Fall back to standard lentiCRISPR v2 vector overhangs
        # (oligos will get CACC and AAAC as the reverse complement)
        return design_sgrna_oligos(
            grna_sequence=grna_sequence,
            five_prime_overhang="GGTG",  # Vector's 5' sticky end → oligo gets CACC
            three_prime_overhang="GTTT",  # Vector's 3' sticky end → oligo gets AAAC
            enzyme="BsmBI",
        )

    # Determine which enzyme was found (use the first one)
    detected_enzyme = sites[0].enzyme if enzyme == "auto" else enzyme

    # Find cloning site pair
    site_pair = find_cloning_site_pair(sites, detected_enzyme)

    if site_pair:
        upstream, downstream = site_pair
        five_prime, three_prime = calculate_sticky_ends(vector_seq, upstream, downstream)
    else:
        # Fall back to standard lentiCRISPR v2 vector overhangs
        five_prime = "GTGG"  # Vector's 5' sticky end → oligo gets CACC
        three_prime = "GTTT"  # Vector's 3' sticky end → oligo gets AAAC

    return design_sgrna_oligos(
        grna_sequence=grna_sequence,
        five_prime_overhang=five_prime,
        three_prime_overhang=three_prime,
        enzyme=detected_enzyme,
    )


# ---------------------------------------------------------------------------
# Helper: Load lentiCRISPR v2 vector
# ---------------------------------------------------------------------------

_LENTICRISPR_V2_PATH = (
    Path(__file__).resolve().parent.parent
    / "Module_Library_gb"
    / "CRISPR Plasmids"
    / "lentiCRISPR v2.gb"
)


def load_lenticrispr_v2() -> str:
    """Load the lentiCRISPR v2 GenBank file."""
    if _LENTICRISPR_V2_PATH.exists():
        return _LENTICRISPR_V2_PATH.read_text(encoding="utf-8")

    # Try alternative path
    alt_path = _LENTICRISPR_V2_PATH.parent / "lentiCRISPR v2 unannotated.gb"
    if alt_path.exists():
        return alt_path.read_text(encoding="utf-8")

    raise FileNotFoundError(f"lentiCRISPR v2 not found at {_LENTICRISPR_V2_PATH}")


def parse_genbank_sequence(genbank_text: str) -> str:
    """Extract the DNA sequence from GenBank format text."""
    origin_idx = genbank_text.find("ORIGIN")
    if origin_idx < 0:
        return normalize_dna(genbank_text)

    seq_section = genbank_text[origin_idx:]
    lines = seq_section.split("\n")[1:]
    seq_parts = []
    for line in lines:
        if line.strip() == "//":
            break
        seq_parts.append(re.sub(r"[\d\s]", "", line))
    return "".join(seq_parts).upper()


def parse_genbank_features(genbank_text: str) -> List[Dict]:
    """Extract features from GenBank format for visualization."""
    features = []

    lines = genbank_text.split("\n")
    in_features = False
    current_feature = None

    for line in lines:
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN"):
            if current_feature:
                features.append(current_feature)
            break

        if not in_features:
            continue

        # New feature line
        if line.startswith("     ") and not line.startswith("                     "):
            if current_feature:
                features.append(current_feature)

            parts = line.strip().split()
            if len(parts) >= 2:
                feat_type = parts[0]
                location = parts[1]

                # Parse location
                direction = -1 if "complement" in location else 1
                numbers = re.findall(r"(\d+)\.\.(\d+)", location)
                if numbers:
                    start = int(numbers[0][0]) - 1  # Convert to 0-based
                    end = int(numbers[0][1])
                    current_feature = {
                        "type": feat_type,
                        "start": start,
                        "end": end,
                        "direction": direction,
                        "qualifiers": {},
                    }

        # Qualifier line
        elif line.startswith("                     /") and current_feature:
            match = re.match(r'\s*/(\w+)="?([^"]*)"?', line)
            if match:
                current_feature["qualifiers"][match.group(1)] = match.group(2)

    return features


# ---------------------------------------------------------------------------
# In Silico Cloning: Assemble the final plasmid
# ---------------------------------------------------------------------------

@dataclass
class AssembledPlasmid:
    """Result of in silico sgRNA cloning assembly."""
    sequence: str
    total_length: int
    grna_sequence: str
    grna_name: str
    enzyme: str
    vector_name: str
    filler_removed_bp: int
    insert_length: int
    annotations: List[Dict]
    restriction_sites: List[Dict]
    ligation_junctions: List[Dict]
    # Predesign-derived metadata (post 2026-05-05). Empty / sentinel values
    # mean the assembler ran in legacy auto-detect mode rather than from a
    # predesign result.
    cassette_kind: str = ""
    promoter_name: str = ""
    promoter_end: int = -1
    upstream_site_pos: int = -1
    downstream_site_pos: int = -1
    five_prime_overhang: str = ""
    three_prime_overhang: str = ""
    insert_start: int = -1
    insert_end: int = -1
    prepend_g: bool = False
    cassette_validation: Dict = field(default_factory=dict)


def assemble_sgrna_plasmid(
    vector_sequence: str,
    grna_sequence: str,
    grna_name: str = "sgRNA",
    enzyme: str = "BsmBI",
    original_features: Optional[List[Dict]] = None,
    *,
    upstream_site_pos: Optional[int] = None,
    downstream_site_pos: Optional[int] = None,
    five_prime_overhang: Optional[str] = None,
    three_prime_overhang: Optional[str] = None,
    cassette_kind: str = "",
    promoter_name: str = "",
    promoter_end: int = -1,
    prepend_g: Optional[bool] = None,
    vector_name: str = "lentiCRISPR v2",
) -> AssembledPlasmid:
    """
    Perform in silico cloning: digest vector and ligate sgRNA insert.

    This simulates:
    1. BsmBI digestion of the vector (removes filler between sites)
    2. Ligation of annealed sgRNA oligos with matching sticky ends

    Args:
        vector_sequence: Full vector sequence (raw DNA)
        grna_sequence: The gRNA target sequence
        grna_name: Name for the gRNA (e.g., "EMX1")
        enzyme: Type IIS enzyme used
        original_features: Features from the original vector to remap

    Returns:
        AssembledPlasmid with the final sequence and annotations
    """
    seq = normalize_dna(vector_sequence)
    grna = normalize_dna(grna_sequence)

    # Resolve the Type IIs cut sites. The predesign step normally supplies
    # them; if not (legacy callers), fall back to the auto-detect path.
    if upstream_site_pos is None or downstream_site_pos is None:
        sites = find_type_iis_sites(seq, enzyme)
        site_pair = find_cloning_site_pair(sites, enzyme)
        if not site_pair:
            raise ValueError(f"Could not find valid {enzyme} site pair for cloning")
        upstream_site, downstream_site = site_pair
        rec_len = len(upstream_site.recognition_seq)
        # Top-strand cut for a forward site is 1 nt past the recognition end;
        # for a reverse site it is 5 nt before the recognition start.
        if upstream_site.orientation == "forward":
            upstream_cut_top = upstream_site.position + rec_len + 1
        else:
            upstream_cut_top = upstream_site.position - 5
        if downstream_site.orientation == "forward":
            downstream_cut_top = downstream_site.position + len(downstream_site.recognition_seq) + 1
        else:
            downstream_cut_top = downstream_site.position - 5
        five_oh = seq[upstream_cut_top:upstream_cut_top + 4]
        three_oh = seq[downstream_cut_top:downstream_cut_top + 4]
        upstream_site_pos = upstream_site.position
        downstream_site_pos = downstream_site.position
        five_prime_overhang = five_oh if five_prime_overhang is None else five_prime_overhang
        three_prime_overhang = three_oh if three_prime_overhang is None else three_prime_overhang

    # Top-strand cut positions (used to splice out the filler).
    # For BsmBI/BbsI/BsaI recognition (6 bp, 1/5 cut offset) we recover the
    # cut top-strand position by inspecting the bases directly, since the
    # caller has only handed us recognition-site starts.
    rec_len_up = 6
    rec_len_down = 6
    for enz_data in TYPE_IIS_ENZYMES.values():
        if seq[upstream_site_pos:upstream_site_pos + len(enz_data["recognition"])] == enz_data["recognition"]:
            rec_len_up = len(enz_data["recognition"])
            up_orientation = "forward"
            break
        if seq[upstream_site_pos:upstream_site_pos + len(enz_data["rc"])] == enz_data["rc"]:
            rec_len_up = len(enz_data["rc"])
            up_orientation = "reverse"
            break
    else:
        up_orientation = "reverse"  # lentiCRISPR convention
    for enz_data in TYPE_IIS_ENZYMES.values():
        if seq[downstream_site_pos:downstream_site_pos + len(enz_data["recognition"])] == enz_data["recognition"]:
            rec_len_down = len(enz_data["recognition"])
            down_orientation = "forward"
            break
        if seq[downstream_site_pos:downstream_site_pos + len(enz_data["rc"])] == enz_data["rc"]:
            rec_len_down = len(enz_data["rc"])
            down_orientation = "reverse"
            break
    else:
        down_orientation = "forward"

    upstream_cut_top = (
        upstream_site_pos + rec_len_up + 1
        if up_orientation == "forward"
        else upstream_site_pos - 5
    )
    downstream_cut_top = (
        downstream_site_pos + rec_len_down + 1
        if down_orientation == "forward"
        else downstream_site_pos - 5
    )

    # The 4-bp top-strand windows at each cut. They sit at positions
    # [upstream_cut_top, upstream_cut_top+4) and [downstream_cut_top,
    # downstream_cut_top+4). After ligation the overhangs are restored from
    # the insert oligos, so the assembled top strand is the vector top strand
    # with the inter-cut stretch replaced by the insert.
    if not five_prime_overhang:
        five_prime_overhang = seq[upstream_cut_top:upstream_cut_top + 4]
    if not three_prime_overhang:
        three_prime_overhang = seq[downstream_cut_top:downstream_cut_top + 4]

    # Insert: for Pol III cassettes, transcription initiates at G. The
    # initiating G may come from the vector's top-strand sticky-end window
    # (when it ends in G) or be prepended by the oligo. Either way the
    # gRNA's last base is dropped so the transcribed spacer length stays at
    # `len(grna)`. T7 fallback callers pass prepend_g=False explicitly.
    if prepend_g is None:
        prepend_g_request = True
    else:
        prepend_g_request = bool(prepend_g)
    insert_seq, prepend_g_eff, _ = adjust_grna_for_pol3(
        grna, prepend_g_request, five_prime_overhang=five_prime_overhang,
    )

    # Filler stretch that disappears on cleavage. Top-strand bases
    # [upstream_cut_top + 4, downstream_cut_top) are released with the
    # filler fragment; the 4-bp windows at each cut stay in the backbone.
    filler_bp = max(0, downstream_cut_top - (upstream_cut_top + 4))

    # Splice: keep top strand up to and including the upstream 4-bp window,
    # add the insert, then resume top strand from the start of the downstream
    # 4-bp window. This avoids the +1 bp duplication that the prior
    # five_prime + three_prime concatenation produced at the downstream end.
    upstream_keep_end = upstream_cut_top + 4  # exclusive end of left fragment in top coords
    downstream_keep_start = downstream_cut_top  # inclusive start of right fragment

    assembled_seq = seq[:upstream_keep_end] + insert_seq + seq[downstream_keep_start:]

    # Insert spans these top-strand coordinates in the assembled sequence.
    insert_position = upstream_keep_end
    delta = len(insert_seq) - (downstream_keep_start - upstream_keep_end)

    # Remap original features
    remapped_annotations = []
    if original_features:
        for feat in original_features:
            label = (
                feat.get("qualifiers", {}).get("label") or
                feat.get("qualifiers", {}).get("gene") or
                feat.get("type", "feature")
            )

            # Skip source and filler features
            if feat.get("type") == "source":
                continue
            if "filler" in str(label).lower():
                continue

            orig_start = feat.get("start", 0)
            orig_end = feat.get("end", 0)

            # Remap based on position relative to cloning site
            if orig_end <= upstream_keep_end:
                # Feature is entirely before cloning site - keep as is
                new_start = orig_start
                new_end = orig_end
            elif orig_start >= downstream_keep_start:
                # Feature is entirely after cloning site - shift by delta
                new_start = orig_start + delta
                new_end = orig_end + delta
            else:
                # Feature overlaps cloning site - skip (it was part of filler)
                continue

            remapped_annotations.append({
                "name": label,
                "start": new_start,
                "end": new_end,
                "direction": feat.get("direction", 1),
                "color": _get_feature_color(label),
            })

    # Add new annotations for the cloned insert
    grna_insert_start = insert_position
    grna_insert_end = insert_position + len(insert_seq)

    # Add gRNA annotation
    remapped_annotations.append({
        "name": f"{grna_name} gRNA",
        "start": grna_insert_start,
        "end": grna_insert_end,
        "direction": 1,
        "color": "#22C55E",  # Green for the gRNA
    })

    # Build ligation junction annotations (sticky ends). The 5' junction sits
    # in the 4 bp immediately UPSTREAM of the insert (the vector's existing
    # bases that re-form after ligation); the 3' junction in the 4 bp
    # immediately DOWNSTREAM. Both labels report the actual overhang sequence
    # from this vector — no hardcoded CACC/GTTT.
    ligation_junctions = [
        {
            "name": f"5' ligation junction ({five_prime_overhang})",
            "start": insert_position - 4,
            "end": insert_position,
            "direction": 1,
            "color": "#F97316",
            "junction_type": "5_prime",
            "overhang_sequence": five_prime_overhang,
        },
        {
            "name": f"3' ligation junction ({three_prime_overhang})",
            "start": grna_insert_end,
            "end": grna_insert_end + 4,
            "direction": 1,
            "color": "#F97316",
            "junction_type": "3_prime",
            "overhang_sequence": three_prime_overhang,
        },
    ]

    # Add junction annotations to main annotations
    for junc in ligation_junctions:
        remapped_annotations.append({
            "name": junc["name"],
            "start": junc["start"],
            "end": junc["end"],
            "direction": junc["direction"],
            "color": junc["color"],
        })

    # Note: In the assembled plasmid, the BsmBI sites are GONE (they were part of filler)
    restriction_sites = []  # No BsmBI sites remain after cloning

    return AssembledPlasmid(
        sequence=assembled_seq,
        total_length=len(assembled_seq),
        grna_sequence=grna,
        grna_name=grna_name,
        enzyme=enzyme,
        vector_name=vector_name,
        filler_removed_bp=filler_bp,
        insert_length=len(insert_seq),
        annotations=remapped_annotations,
        restriction_sites=restriction_sites,
        ligation_junctions=ligation_junctions,
        cassette_kind=cassette_kind,
        promoter_name=promoter_name,
        promoter_end=promoter_end,
        upstream_site_pos=upstream_site_pos,
        downstream_site_pos=downstream_site_pos,
        five_prime_overhang=five_prime_overhang,
        three_prime_overhang=three_prime_overhang,
        insert_start=grna_insert_start,
        insert_end=grna_insert_end,
        prepend_g=prepend_g_eff,
    )


def _get_feature_color(label: str) -> str:
    """Get a color for a feature based on its label."""
    label_lower = str(label).lower()

    color_map = {
        "u6": "#06B6D4",  # Cyan for U6 promoter
        "grna scaffold": "#10B981",  # Emerald for scaffold
        "cas9": "#3B82F6",  # Blue for Cas9
        "puror": "#F59E0B",  # Amber for selection
        "bleor": "#F59E0B",
        "ampr": "#F59E0B",
        "cmv": "#EC4899",  # Pink for CMV
        "ef-1": "#8B5CF6",  # Purple for EF1a
        "ltr": "#64748B",  # Slate for LTR
        "wpre": "#14B8A6",  # Teal for WPRE
        "ori": "#6366F1",  # Indigo for ori
        "psi": "#A855F7",  # Purple for Psi
        "rre": "#A855F7",
        "flag": "#EF4444",  # Red for tags
        "nls": "#EF4444",
        "p2a": "#F472B6",  # Pink for 2A
    }

    for key, color in color_map.items():
        if key in label_lower:
            return color

    return "#86A8D9"  # Default blue-gray
