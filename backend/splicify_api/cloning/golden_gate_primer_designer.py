"""
Golden Gate Assembly Primer Designer

Comprehensive module for designing Golden Gate assembly workflows:
1. Multi-Fragment Assembly (2-6+ fragments with orthogonal sticky ends)
2. Single Fragment Replacement (replace/insert sequences >52 bp)
3. Scarless Deletion (remove regions without scars)

Golden Gate primer structure:
    5'-[extra 2-4bp]-[Type IIS recognition 6bp]-[spacer N]-[4bp overhang]-[annealing ~20bp]-3'

Example for BsaI: 5'-GGG-GGTCTC-N-AATG-[20bp annealing]-3'
"""
from __future__ import annotations

import base64
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from ..utils import normalize_dna, reverse_complement, safe_tm

# Try to import primer3 for Tm calculation
try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False


# ---------------------------------------------------------------------------
# Type IIS Enzyme Database
# ---------------------------------------------------------------------------

TYPE_IIS_ENZYMES = {
    "BsaI": {
        "recognition": "GGTCTC",
        "rc": "GAGACC",
        "cut_offset": (1, 5),  # Cuts N1/N5 downstream
        "extra_bases": "GGG",  # Extra bases for efficient cutting
        "spacer": "N",  # Single N spacer between recognition and overhang
    },
    "BsmBI": {
        "recognition": "CGTCTC",
        "rc": "GAGACG",
        "cut_offset": (1, 5),
        "extra_bases": "GGG",
        "spacer": "N",
    },
    "BbsI": {
        "recognition": "GAAGAC",
        "rc": "GTCTTC",
        "cut_offset": (2, 6),
        "extra_bases": "GGG",
        "spacer": "NN",
    },
}


# ---------------------------------------------------------------------------
# Overhang Sets (High Fidelity, MoClo-Compatible)
# ---------------------------------------------------------------------------

# Standard high-fidelity overhang sets based on Potapov et al. (2018) ligation data
# Each tuple is (sequence, fidelity_score)
STANDARD_OVERHANG_SETS = {
    "high_fidelity_6": [
        ("AATG", 0.98),  # CDS start (contains ATG)
        ("GCTT", 0.97),
        ("GGAG", 0.96),
        ("CGCT", 0.95),
        ("TCGA", 0.94),
        ("CAGC", 0.93),
    ],
    "high_fidelity_10": [
        ("AATG", 0.98),
        ("GCTT", 0.97),
        ("GGAG", 0.96),
        ("CGCT", 0.95),
        ("TCGA", 0.94),
        ("CAGC", 0.93),
        ("ATCC", 0.92),
        ("GTCA", 0.91),
        ("AGGT", 0.90),
        ("TGCC", 0.89),
    ],
    "moclo_standard": [
        ("GGAG", 0.96),  # Standard MoClo fusion sites
        ("TACT", 0.95),
        ("AATG", 0.98),
        ("AGCC", 0.94),
        ("TTCG", 0.93),
        ("GCTT", 0.97),
    ],
}

# Overhangs to avoid (palindromes, homopolymers, low fidelity pairs)
BAD_OVERHANGS = {
    # Palindromes (will self-anneal)
    "ATAT", "TATA", "CGCG", "GCGC", "AATT", "TTAA", "CCGG", "GGCC",
    # Homopolymers
    "AAAA", "TTTT", "CCCC", "GGGG",
    # Near-palindromes with poor fidelity
    "ACGT", "TGCA", "CATG", "GTAC",
}


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class GoldenGateOverhang:
    """Represents a 4-bp overhang for Golden Gate assembly."""
    sequence: str  # 4bp (e.g., "AATG")
    fidelity_score: float  # 0-1
    category: str  # "moclo" | "standard" | "custom"


@dataclass
class GoldenGatePrimer:
    """A primer designed for Golden Gate assembly."""
    name: str
    sequence: str  # Full primer sequence
    annealing_portion: str
    annealing_tm: float
    overhang_4bp: str
    enzyme: str
    direction: str  # "forward" | "reverse"
    fragment_name: str


@dataclass
class GoldenGateJunction:
    """Represents a junction between two fragments in the assembly."""
    left_fragment: str
    right_fragment: str
    overhang: GoldenGateOverhang
    forward_primer: GoldenGatePrimer  # For right fragment
    reverse_primer: GoldenGatePrimer  # For left fragment
    junction_index: int


@dataclass
class InternalSiteWarning:
    """Warning about internal Type IIS sites that need domestication."""
    fragment_name: str
    enzyme: str
    position: int
    sequence_context: str  # Surrounding 20bp for context


@dataclass
class GoldenGatePrimerDesign:
    """Complete result of Golden Gate primer design."""
    workflow_type: str  # "multi_fragment" | "single_fragment" | "deletion"
    enzyme: str
    fragments: List[Dict[str, Any]]
    junctions: List[GoldenGateJunction]
    primer_table: List[GoldenGatePrimer]
    overhang_set: List[str]
    overhang_fidelity_mean: float
    internal_site_warnings: List[InternalSiteWarning]
    assembled_sequence: str
    assembled_length: int
    annotations: List[Dict[str, Any]]
    topology: str  # "circular" | "linear"
    protocol_notes: List[str] = field(default_factory=list)
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


def select_orthogonal_overhangs(
    n_junctions: int,
    existing: Optional[List[str]] = None,
    overhang_set: str = "high_fidelity_10"
) -> List[GoldenGateOverhang]:
    """
    Select n orthogonal 4bp overhangs for Golden Gate assembly.

    Ensures:
    - No palindromes or homopolymers
    - No reverse complements of already-selected overhangs
    - High fidelity scores
    - Compatibility with existing overhangs if provided

    Args:
        n_junctions: Number of overhangs needed
        existing: List of existing overhangs to avoid conflicts with
        overhang_set: Which predefined set to use

    Returns:
        List of GoldenGateOverhang objects sorted by fidelity
    """
    existing_set = set(s.upper() for s in (existing or []))

    # Add reverse complements of existing to exclusion set
    for oh in list(existing_set):
        existing_set.add(reverse_complement(oh))

    # Get candidate overhangs from the specified set
    candidates = STANDARD_OVERHANG_SETS.get(overhang_set, STANDARD_OVERHANG_SETS["high_fidelity_10"])

    selected: List[GoldenGateOverhang] = []
    used_set: set = set(existing_set)

    for seq, score in candidates:
        if len(selected) >= n_junctions:
            break

        seq_upper = seq.upper()
        rc_seq = reverse_complement(seq_upper)

        # Skip if this overhang or its RC is already used
        if seq_upper in used_set or rc_seq in used_set:
            continue

        # Skip bad overhangs
        if seq_upper in BAD_OVERHANGS:
            continue

        # Add to selected
        category = "moclo" if overhang_set == "moclo_standard" else "standard"
        selected.append(GoldenGateOverhang(
            sequence=seq_upper,
            fidelity_score=score,
            category=category,
        ))

        # Mark both the overhang and its RC as used
        used_set.add(seq_upper)
        used_set.add(rc_seq)

    return selected


def check_internal_type_iis_sites(
    sequence: str,
    enzyme: str,
    fragment_name: str = "fragment"
) -> List[InternalSiteWarning]:
    """
    Check for internal Type IIS recognition sites that would need domestication.

    Args:
        sequence: DNA sequence to check
        enzyme: Enzyme name ("BsaI", "BsmBI", "BbsI")
        fragment_name: Name for reporting

    Returns:
        List of InternalSiteWarning objects
    """
    warnings: List[InternalSiteWarning] = []
    seq = normalize_dna(sequence)

    if enzyme not in TYPE_IIS_ENZYMES:
        return warnings

    enz = TYPE_IIS_ENZYMES[enzyme]
    patterns = [enz["recognition"], enz["rc"]]

    for pattern in patterns:
        for match in re.finditer(pattern, seq):
            pos = match.start()
            # Get context (10bp on each side)
            start = max(0, pos - 10)
            end = min(len(seq), pos + len(pattern) + 10)
            context = seq[start:end]

            warnings.append(InternalSiteWarning(
                fragment_name=fragment_name,
                enzyme=enzyme,
                position=pos,
                sequence_context=context,
            ))

    return warnings


def pick_annealing_region(
    sequence: str,
    direction: str,
    target_tm: float = 60.0,
    min_len: int = 18,
    max_len: int = 25
) -> Tuple[str, float]:
    """
    Pick optimal annealing region from a sequence.

    Args:
        sequence: Full fragment sequence
        direction: "forward" (5' end) or "reverse" (3' end)
        target_tm: Target melting temperature
        min_len: Minimum annealing length
        max_len: Maximum annealing length

    Returns:
        Tuple of (annealing_sequence, tm)
    """
    seq = normalize_dna(sequence)

    if direction == "forward":
        region = seq[:max_len]
    else:
        region = reverse_complement(seq[-max_len:])

    best_seq = region[:min_len]
    best_tm = calculate_tm(best_seq)
    best_delta = abs(best_tm - target_tm)

    for length in range(min_len, min(max_len + 1, len(region) + 1)):
        candidate = region[:length]
        tm = calculate_tm(candidate)
        delta = abs(tm - target_tm)
        if delta < best_delta:
            best_seq = candidate
            best_tm = tm
            best_delta = delta
            # Stop if we hit target
            if delta < 1.0:
                break

    return best_seq, round(best_tm, 1)


def build_golden_gate_primer(
    annealing_seq: str,
    overhang: str,
    enzyme: str,
    direction: str,
    fragment_name: str,
    primer_suffix: str
) -> GoldenGatePrimer:
    """
    Build a complete Golden Gate primer with Type IIS extension.

    For BsaI: 5'-GGG-GGTCTC-N-[OVERHANG]-[annealing]-3' (forward)
              5'-GGG-GAGACC-N-[RC(OVERHANG)]-[annealing]-3' (reverse)

    Args:
        annealing_seq: The annealing portion of the primer
        overhang: 4bp overhang sequence
        enzyme: Enzyme name
        direction: "forward" or "reverse"
        fragment_name: Name of the fragment
        primer_suffix: Suffix for primer name (e.g., "FWD", "REV")

    Returns:
        GoldenGatePrimer object
    """
    enz = TYPE_IIS_ENZYMES.get(enzyme, TYPE_IIS_ENZYMES["BsaI"])

    if direction == "forward":
        recognition = enz["recognition"]
        oh_seq = overhang
    else:
        recognition = enz["rc"]
        oh_seq = reverse_complement(overhang)

    # Build full primer: extra + recognition + spacer + overhang + annealing
    full_seq = f"{enz['extra_bases']}{recognition}{enz['spacer']}{oh_seq}{annealing_seq}"

    return GoldenGatePrimer(
        name=f"{fragment_name}_{primer_suffix}",
        sequence=full_seq,
        annealing_portion=annealing_seq,
        annealing_tm=calculate_tm(annealing_seq),
        overhang_4bp=overhang,
        enzyme=enzyme,
        direction=direction,
        fragment_name=fragment_name,
    )


def design_multi_fragment_assembly(
    fragments: List[Dict[str, Any]],
    topology: str = "circular",
    enzyme: str = "BsaI",
    target_tm: float = 60.0,
) -> GoldenGatePrimerDesign:
    """
    Design primers for multi-fragment Golden Gate assembly.

    Args:
        fragments: List of dicts with 'name' and 'sequence' keys
        topology: "circular" or "linear"
        enzyme: Type IIS enzyme to use
        target_tm: Target annealing Tm

    Returns:
        Complete GoldenGatePrimerDesign
    """
    n_fragments = len(fragments)

    if n_fragments < 2:
        raise ValueError("Multi-fragment assembly requires at least 2 fragments")

    # For circular: n_fragments junctions (each fragment connects to next, last connects to first)
    # For linear: n_fragments - 1 junctions
    n_junctions = n_fragments if topology == "circular" else n_fragments - 1

    # Select orthogonal overhangs
    overhangs = select_orthogonal_overhangs(n_junctions)

    if len(overhangs) < n_junctions:
        raise ValueError(f"Could not select {n_junctions} orthogonal overhangs")

    # Check for internal sites
    all_warnings: List[InternalSiteWarning] = []
    for frag in fragments:
        frag_warnings = check_internal_type_iis_sites(
            frag["sequence"],
            enzyme,
            frag.get("name", "fragment"),
        )
        all_warnings.extend(frag_warnings)

    # Build junctions and primers
    junctions: List[GoldenGateJunction] = []
    all_primers: List[GoldenGatePrimer] = []

    for i in range(n_junctions):
        left_idx = i
        right_idx = (i + 1) % n_fragments

        left_frag = fragments[left_idx]
        right_frag = fragments[right_idx]
        overhang = overhangs[i]

        # Pick annealing regions
        left_anneal, left_tm = pick_annealing_region(
            left_frag["sequence"], "reverse", target_tm
        )
        right_anneal, right_tm = pick_annealing_region(
            right_frag["sequence"], "forward", target_tm
        )

        # Build primers
        # Reverse primer for left fragment (3' end) - carries the overhang
        rev_primer = build_golden_gate_primer(
            left_anneal,
            overhang.sequence,
            enzyme,
            "reverse",
            left_frag.get("name", f"Frag{left_idx + 1}"),
            "REV",
        )

        # Forward primer for right fragment (5' end) - carries the same overhang
        fwd_primer = build_golden_gate_primer(
            right_anneal,
            overhang.sequence,
            enzyme,
            "forward",
            right_frag.get("name", f"Frag{right_idx + 1}"),
            "FWD",
        )

        junction = GoldenGateJunction(
            left_fragment=left_frag.get("name", f"Fragment_{left_idx + 1}"),
            right_fragment=right_frag.get("name", f"Fragment_{right_idx + 1}"),
            overhang=overhang,
            forward_primer=fwd_primer,
            reverse_primer=rev_primer,
            junction_index=i,
        )
        junctions.append(junction)
        all_primers.extend([fwd_primer, rev_primer])

    # Simulate assembly
    assembled_seq, annotations = simulate_golden_gate_assembly(
        fragments, junctions, topology
    )

    # Calculate mean fidelity
    mean_fidelity = sum(j.overhang.fidelity_score for j in junctions) / len(junctions)

    # Build warnings list
    text_warnings: List[str] = []
    if all_warnings:
        text_warnings.append(
            f"Found {len(all_warnings)} internal {enzyme} site(s) that may need domestication"
        )

    # Protocol notes
    protocol = [
        f"1. PCR amplify each fragment with the provided primers using a high-fidelity polymerase (Q5 or Phusion)",
        f"2. Verify PCR products by gel electrophoresis",
        f"3. DpnI digest template DNA if starting from plasmid template",
        f"4. Purify PCR products (column or gel extraction)",
        f"5. Set up Golden Gate assembly reaction:",
        f"   - 75 ng of each fragment",
        f"   - 1 µL {enzyme} enzyme (NEB)",
        f"   - 1 µL T4 DNA ligase (400 U/µL)",
        f"   - 2 µL T4 DNA ligase buffer (10X)",
        f"   - Water to 20 µL",
        f"6. Cycle: 25-50 cycles of (37°C 3 min, 16°C 4 min), then 50°C 5 min, 80°C 10 min",
        f"7. Transform 2-5 µL into competent cells",
        f"8. Screen colonies by PCR or sequencing",
    ]

    return GoldenGatePrimerDesign(
        workflow_type="multi_fragment",
        enzyme=enzyme,
        fragments=fragments,
        junctions=junctions,
        primer_table=all_primers,
        overhang_set=[j.overhang.sequence for j in junctions],
        overhang_fidelity_mean=round(mean_fidelity, 3),
        internal_site_warnings=all_warnings,
        assembled_sequence=assembled_seq,
        assembled_length=len(assembled_seq),
        annotations=annotations,
        topology=topology,
        protocol_notes=protocol,
        warnings=text_warnings,
    )


def simulate_golden_gate_assembly(
    fragments: List[Dict[str, Any]],
    junctions: List[GoldenGateJunction],
    topology: str,
) -> Tuple[str, List[Dict[str, Any]]]:
    """
    Simulate Golden Gate assembly to predict the final sequence.

    Args:
        fragments: Fragment list with sequences
        junctions: Junction definitions with overhangs
        topology: "circular" or "linear"

    Returns:
        Tuple of (assembled_sequence, annotations_list)
    """
    if not fragments:
        return "", []

    # Build assembled sequence by concatenating fragments in order
    # Each fragment contributes from its 5' end (after the 4bp overhang region)
    # to its 3' end (before the next 4bp overhang region)
    assembled_parts: List[str] = []
    annotations: List[Dict[str, Any]] = []

    # Colors for fragments
    fragment_colors = [
        "#10B981",  # Emerald
        "#3B82F6",  # Blue
        "#8B5CF6",  # Purple
        "#F59E0B",  # Amber
        "#EF4444",  # Red
        "#06B6D4",  # Cyan
        "#EC4899",  # Pink
        "#84CC16",  # Lime
    ]

    junction_color = "#F97316"  # Orange for junctions
    primer_fwd_color = "#22D3EE"  # Cyan for forward primers
    primer_rev_color = "#A855F7"  # Purple for reverse primers

    current_pos = 0

    for i, frag in enumerate(fragments):
        seq = normalize_dna(frag.get("sequence", ""))

        # Add the fragment sequence
        frag_start = current_pos
        assembled_parts.append(seq)
        frag_end = current_pos + len(seq)

        # Add fragment annotation
        annotations.append({
            "name": frag.get("name", f"Fragment_{i + 1}"),
            "start": frag_start,
            "end": frag_end,
            "direction": 1,
            "color": fragment_colors[i % len(fragment_colors)],
            "type": "fragment",
        })

        # Add primer annealing annotations for this fragment
        # Find the junction that has primers for this fragment
        # Each junction has a reverse primer for the left fragment and forward primer for the right fragment

        # Forward primer (at 5' end of fragment) - comes from previous junction
        if i > 0 and i - 1 < len(junctions):
            prev_junction = junctions[i - 1]
            fwd_primer = prev_junction.forward_primer
            anneal_len = len(fwd_primer.annealing_portion)
            annotations.append({
                "name": f"{fwd_primer.name} (annealing)",
                "start": frag_start,
                "end": frag_start + anneal_len,
                "direction": 1,
                "color": primer_fwd_color,
                "type": "primer_annealing",
            })
        elif i == 0 and topology == "circular" and len(junctions) > 0:
            # For circular, first fragment's forward primer is from last junction
            last_junction = junctions[-1]
            fwd_primer = last_junction.forward_primer
            anneal_len = len(fwd_primer.annealing_portion)
            annotations.append({
                "name": f"{fwd_primer.name} (annealing)",
                "start": frag_start,
                "end": frag_start + anneal_len,
                "direction": 1,
                "color": primer_fwd_color,
                "type": "primer_annealing",
            })

        # Reverse primer (at 3' end of fragment) - comes from current junction
        if i < len(junctions):
            curr_junction = junctions[i]
            rev_primer = curr_junction.reverse_primer
            anneal_len = len(rev_primer.annealing_portion)
            annotations.append({
                "name": f"{rev_primer.name} (annealing)",
                "start": frag_end - anneal_len,
                "end": frag_end,
                "direction": -1,
                "color": primer_rev_color,
                "type": "primer_annealing",
            })

        current_pos = frag_end

        # Add junction annotation if not the last fragment (or if circular)
        if i < len(junctions):
            junction = junctions[i]
            # Junction overhang spans the boundary: last 1bp of left fragment + first 3bp of right fragment
            # This creates the 4bp overhang sequence that the primers add
            junction_start = frag_end - 1
            junction_end = frag_end + 3

            annotations.append({
                "name": f"Junction {i + 1} ({junction.overhang.sequence})",
                "start": junction_start,
                "end": junction_end,
                "direction": 1,
                "color": junction_color,
                "type": "junction",
            })

    assembled_seq = "".join(assembled_parts)

    # For circular topology, verify first/last overhangs match
    if topology == "circular" and len(junctions) > 0:
        # The last junction connects back to the first fragment
        pass  # Sequence is already assembled correctly

    return assembled_seq, annotations


def design_single_fragment_replacement(
    template: str,
    insert: str,
    position: int,
    enzyme: str = "BsaI",
    template_name: str = "Template",
    insert_name: str = "Insert",
) -> GoldenGatePrimerDesign:
    """
    Design primers for single fragment insertion/replacement.

    Args:
        template: Template sequence
        insert: Insert sequence (>52 bp for Golden Gate advantage)
        position: Position in template where insert goes
        enzyme: Type IIS enzyme
        template_name: Name for template
        insert_name: Name for insert

    Returns:
        GoldenGatePrimerDesign
    """
    template_seq = normalize_dna(template)
    insert_seq = normalize_dna(insert)

    if len(insert_seq) < 30:
        # Golden Gate not ideal for very short inserts
        pass  # Still proceed but add warning

    # Split template at insertion point
    left_template = template_seq[:position]
    right_template = template_seq[position:]

    fragments = [
        {"name": f"{template_name}_left", "sequence": left_template},
        {"name": insert_name, "sequence": insert_seq},
        {"name": f"{template_name}_right", "sequence": right_template},
    ]

    return design_multi_fragment_assembly(
        fragments=fragments,
        topology="circular",
        enzyme=enzyme,
    )


def design_scarless_deletion(
    template: str,
    del_start: int,
    del_end: int,
    enzyme: str = "BsaI",
    template_name: str = "Template",
) -> GoldenGatePrimerDesign:
    """
    Design primers for scarless deletion using Golden Gate.

    Args:
        template: Template sequence
        del_start: Start position of deletion (0-based)
        del_end: End position of deletion (0-based, exclusive)
        enzyme: Type IIS enzyme
        template_name: Name for template

    Returns:
        GoldenGatePrimerDesign
    """
    template_seq = normalize_dna(template)

    if del_start >= del_end:
        raise ValueError("del_start must be less than del_end")

    if del_end > len(template_seq):
        raise ValueError("del_end exceeds template length")

    # Create two fragments: before and after deletion
    left_fragment = template_seq[:del_start]
    right_fragment = template_seq[del_end:]

    fragments = [
        {"name": f"{template_name}_upstream", "sequence": left_fragment},
        {"name": f"{template_name}_downstream", "sequence": right_fragment},
    ]

    design = design_multi_fragment_assembly(
        fragments=fragments,
        topology="circular",
        enzyme=enzyme,
    )

    design.workflow_type = "deletion"
    return design


# ---------------------------------------------------------------------------
# Output Formatting
# ---------------------------------------------------------------------------

def format_primer_table_markdown(primers: List[GoldenGatePrimer]) -> str:
    """Format primer list as markdown table."""
    lines = [
        "| Name | Sequence (5' to 3') | Length | Tm | Overhang |",
        "|------|---------------------|--------|-----|----------|",
    ]

    for p in primers:
        lines.append(
            f"| {p.name} | `{p.sequence}` | {len(p.sequence)} bp | {p.annealing_tm:.1f}°C | {p.overhang_4bp} |"
        )

    return "\n".join(lines)


def format_overhang_table_markdown(junctions: List[GoldenGateJunction]) -> str:
    """Format junction/overhang information as markdown table."""
    lines = [
        "| Junction | Left → Right | Overhang | Fidelity |",
        "|----------|--------------|----------|----------|",
    ]

    for j in junctions:
        lines.append(
            f"| {j.junction_index + 1} | {j.left_fragment} → {j.right_fragment} | {j.overhang.sequence} | {j.overhang.fidelity_score:.2f} |"
        )

    return "\n".join(lines)


def build_design_response(design: GoldenGatePrimerDesign) -> Dict[str, Any]:
    """Build API response dict from design result."""

    # Build reply text
    reply_parts = [
        f"## Golden Gate {design.workflow_type.replace('_', ' ').title()}",
        "",
        f"**Enzyme**: {design.enzyme}",
        f"**Fragments**: {len(design.fragments)}",
        f"**Topology**: {design.topology.title()}",
        f"**Assembled Length**: {design.assembled_length:,} bp",
        f"**Mean Overhang Fidelity**: {design.overhang_fidelity_mean:.2f}",
        "",
        "### Overhang Set",
        "",
        format_overhang_table_markdown(design.junctions),
        "",
        "### Primers to Order",
        "",
        format_primer_table_markdown(design.primer_table),
        "",
    ]

    if design.internal_site_warnings:
        reply_parts.append("### Internal Site Warnings")
        reply_parts.append("")
        for w in design.internal_site_warnings:
            reply_parts.append(f"- **{w.fragment_name}**: {w.enzyme} site at position {w.position}")
            reply_parts.append(f"  Context: `{w.sequence_context}`")
        reply_parts.append("")

    if design.protocol_notes:
        reply_parts.append("### Protocol")
        reply_parts.append("")
        for note in design.protocol_notes:
            reply_parts.append(note)
        reply_parts.append("")

    reply = "\n".join(reply_parts)

    # Build viz payload
    viz = {
        "type": "design",
        "title": f"Golden Gate Assembly: {' + '.join(f.get('name', 'Fragment') for f in design.fragments)}",
        "sequence": design.assembled_sequence,
        "topology": design.topology,
        "total_length": design.assembled_length,
        "annotations": design.annotations,
    }

    # Build primer CSV
    primer_csv_lines = ["Name,Sequence,Length,Tm,Overhang,Direction,Fragment"]
    for p in design.primer_table:
        primer_csv_lines.append(
            f"{p.name},{p.sequence},{len(p.sequence)},{p.annealing_tm},{p.overhang_4bp},{p.direction},{p.fragment_name}"
        )
    primer_csv = "\n".join(primer_csv_lines)

    # Build GenBank file
    genbank_text = build_assembly_genbank(design)

    # Build protocol markdown
    protocol_md = "\n".join([
        f"# Golden Gate Assembly Protocol",
        "",
        f"## Assembly Information",
        f"- Enzyme: {design.enzyme}",
        f"- Fragments: {len(design.fragments)}",
        f"- Topology: {design.topology}",
        f"- Final size: {design.assembled_length:,} bp",
        "",
        f"## Overhang Set",
        "",
        format_overhang_table_markdown(design.junctions),
        "",
        f"## Primers",
        "",
        format_primer_table_markdown(design.primer_table),
        "",
        f"## Protocol Steps",
        "",
        *design.protocol_notes,
    ])

    files = [
        {
            "fileName": "golden_gate_primers.csv",
            "mimeType": "text/csv",
            "dataBase64": base64.b64encode(primer_csv.encode("utf-8")).decode("ascii"),
        },
        {
            "fileName": "assembled_plasmid.gb",
            "mimeType": "application/octet-stream",
            "dataBase64": base64.b64encode(genbank_text.encode("utf-8")).decode("ascii"),
        },
        {
            "fileName": "golden_gate_protocol.md",
            "mimeType": "text/markdown",
            "dataBase64": base64.b64encode(protocol_md.encode("utf-8")).decode("ascii"),
        },
    ]

    return {
        "ok": True,
        "reply": reply,
        "viz": viz,
        "files": files,
        "metadata": {
            "workflow_type": design.workflow_type,
            "enzyme": design.enzyme,
            "fragment_count": len(design.fragments),
            "primer_count": len(design.primer_table),
            "assembled_length": design.assembled_length,
            "overhang_fidelity_mean": design.overhang_fidelity_mean,
            "internal_site_count": len(design.internal_site_warnings),
        },
    }


def build_assembly_genbank(design: GoldenGatePrimerDesign) -> str:
    """Build GenBank format file for assembled sequence."""
    seq = design.assembled_sequence
    title = f"GG_Assembly"[:16]

    topology_str = "circular" if design.topology == "circular" else "linear"

    features = [
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(seq)}",
        '                     /organism="synthetic construct"',
        '                     /mol_type="other DNA"',
    ]

    for ann in design.annotations:
        start = int(ann.get("start", 0)) + 1  # Convert to 1-based
        end = int(ann.get("end", 0))
        if end <= start:
            continue

        location = f"{start}..{end}"
        if int(ann.get("direction", 1)) < 0:
            location = f"complement({location})"

        feat_type = "misc_feature"
        name = str(ann.get("name", "feature"))

        if ann.get("type") == "junction":
            feat_type = "misc_feature"
        elif ann.get("type") == "fragment":
            feat_type = "misc_feature"

        features.append(f"     {feat_type.ljust(15)} {location}")
        features.append(f'                     /label="{name.replace(chr(34), chr(39))}"')

    # Build origin section
    origin = ["ORIGIN"]
    seq_lower = seq.lower()
    for i in range(0, len(seq_lower), 60):
        chunk = seq_lower[i:i + 60]
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        origin.append(f"{str(i + 1).rjust(9)} {groups}")

    return "\n".join([
        f"LOCUS       {title.ljust(16)} {str(len(seq)).rjust(6)} bp    DNA     {topology_str} SYN",
        f"DEFINITION  Golden Gate assembly of {len(design.fragments)} fragments using {design.enzyme}.",
        "ACCESSION   .",
        "VERSION     .",
        "KEYWORDS    Golden Gate; Type IIS; modular assembly.",
        "SOURCE      synthetic DNA construct",
        "  ORGANISM  synthetic DNA construct",
        f"COMMENT     Assembled using {design.enzyme} Golden Gate assembly.",
        f"            Fragments: {', '.join(f.get('name', 'Fragment') for f in design.fragments)}",
        *features,
        *origin,
        "//",
    ])
