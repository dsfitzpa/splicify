"""
Junction — the atomic boundary between two adjacent modules.

All cloning operators act on junctions, not on whole plasmids.
Computed once from resolved modules; shared by all operator evaluators.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Default context window: bp of flanking sequence extracted per side.
# Must be >= (max_overlap_len + max_anneal_len) for primer design.
CONTEXT_WINDOW = 100

# Restriction enzymes to pre-scan at every junction.
# Relevant for: restriction operator feasibility, Golden Gate site detection.
_RE_SITES: Dict[str, str] = {
    # Common Type II (directional cloning)
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NheI": "GCTAGC",
    "XhoI": "CTCGAG",
    "NcoI": "CCATGG",
    "NotI": "GCGGCCGC",
    "XbaI": "TCTAGA",
    "SpeI": "ACTAGT",
    "PstI": "CTGCAG",
    "SalI": "GTCGAC",
    "KpnI": "GGTACC",
    "SmaI": "CCCGGG",
    "AgeI": "ACCGGT",
    "ClaI": "ATCGAT",
    "MluI": "ACGCGT",
    # Type IIs (Golden Gate)
    "BsmBI": "CGTCTC",
    "BsaI": "GGTCTC",
    "BbsI": "GAAGAC",
    "SapI": "GCTCTTC",
    "BtgZI": "GCGATG",
}

# Roles that carry open reading frames (reading-frame-sensitive junctions)
_CODING_ROLES = frozenset({
    "transgene", "reporter", "nuclease", "selection_marker",
    "backbone",  # sometimes contains ORFs (AmpR, etc.)
})
# Roles at regulatory boundaries
_PROMOTER_ROLES = frozenset({"promoter"})
_TERMINATOR_ROLES = frozenset({"polya", "terminator"})
_CDS_ROLES = frozenset({"transgene", "reporter", "nuclease", "selection_marker"})


@dataclass
class Junction:
    """
    Represents the boundary between two adjacent modules in a plasmid design.

    Computed once by build_junctions(); consumed by all cloning operators for
    feasibility analysis, primer design, and RE site detection.
    """
    junction_index: int
    left_module_index: int
    right_module_index: int

    # Human-readable names for display
    left_module_name: str
    right_module_name: str
    left_module_role: str
    right_module_role: str

    # Sequence context (up to CONTEXT_WINDOW bp each side)
    left_flank: str    # last N bp of left module sequence
    right_flank: str   # first N bp of right module sequence
    junction_context: str  # left_flank + right_flank (for site scanning)

    # Biological constraints
    reading_frame_continuation: bool = False  # ORF crosses this junction?
    scar_allowed: bool = True                 # can a scar/overhang live here?
    regulatory_boundary: bool = False         # promoter→CDS or CDS→polyA?

    # Pre-computed RE site scan (sites found within junction_context)
    internal_restriction_sites: Dict[str, List[int]] = field(default_factory=dict)

    # Sequence complexity flags
    has_tandem_repeat: bool = False
    repeat_note: str = ""
    left_gc_content: float = 0.0
    right_gc_content: float = 0.0
    left_has_homopolymer: bool = False
    right_has_homopolymer: bool = False


def build_junctions(
    modules: List[dict],
    topology: str = "circular",
    context_window: int = CONTEXT_WINDOW,
) -> List[Junction]:
    """
    Build Junction objects from a list of resolved module dicts.

    Args:
        modules: resolved module dicts (from plasmid_design_chat.py resolution step).
                 Each must have at minimum: 'sequence', 'role', and one of
                 'canonical_id' or 'description'.
        topology: 'circular' (n junctions) or 'linear' (n-1 junctions).
        context_window: bp of flanking sequence to extract per side.

    Returns:
        List of Junction objects in order around the construct.
    """
    n = len(modules)
    if n < 2:
        return []

    n_junctions = n if topology == "circular" else n - 1
    junctions: List[Junction] = []

    for i in range(n_junctions):
        left_idx = i
        right_idx = (i + 1) % n

        left_mod = modules[left_idx]
        right_mod = modules[right_idx]

        left_seq = left_mod.get("sequence") or ""
        right_seq = right_mod.get("sequence") or ""

        # Extract flanking context
        left_flank = left_seq[-context_window:] if len(left_seq) >= context_window else left_seq
        right_flank = right_seq[:context_window] if len(right_seq) >= context_window else right_seq
        context = left_flank + right_flank

        # Restriction site scan over junction context
        re_sites = _scan_restriction_sites(context)

        # Sequence complexity
        left_gc = _gc_content(left_flank)
        right_gc = _gc_content(right_flank)
        left_hp = _has_homopolymer(left_flank)
        right_hp = _has_homopolymer(right_flank)

        # Tandem repeat risk
        has_repeat, repeat_note = _detect_junction_repeat(left_seq, right_seq)

        # Role-based constraint inference
        left_role = left_mod.get("role", "")
        right_role = right_mod.get("role", "")
        reading_frame = _is_coding_junction(left_role, right_role, left_seq, right_seq)
        scar_allowed = not reading_frame
        regulatory = _is_regulatory_boundary(left_role, right_role)

        # Module display name: prefer canonical_id, fall back to description
        left_name = (
            left_mod.get("canonical_id")
            or left_mod.get("description")
            or f"module_{left_idx}"
        )
        right_name = (
            right_mod.get("canonical_id")
            or right_mod.get("description")
            or f"module_{right_idx}"
        )

        junctions.append(Junction(
            junction_index=i,
            left_module_index=left_idx,
            right_module_index=right_idx,
            left_module_name=left_name,
            right_module_name=right_name,
            left_module_role=left_role,
            right_module_role=right_role,
            left_flank=left_flank,
            right_flank=right_flank,
            junction_context=context,
            reading_frame_continuation=reading_frame,
            scar_allowed=scar_allowed,
            regulatory_boundary=regulatory,
            internal_restriction_sites=re_sites,
            has_tandem_repeat=has_repeat,
            repeat_note=repeat_note,
            left_gc_content=left_gc,
            right_gc_content=right_gc,
            left_has_homopolymer=left_hp,
            right_has_homopolymer=right_hp,
        ))

    return junctions


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _scan_restriction_sites(context: str) -> Dict[str, List[int]]:
    """Find all known RE recognition sites within the junction context."""
    ctx = context.upper()
    found: Dict[str, List[int]] = {}
    for name, site in _RE_SITES.items():
        positions = [m.start() for m in re.finditer(site, ctx)]
        rc = _reverse_complement(site)
        if rc != site:
            positions += [m.start() for m in re.finditer(rc, ctx)]
        if positions:
            found[name] = sorted(positions)
    return found


def _gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s)


def _has_homopolymer(seq: str, min_run: int = 5) -> bool:
    """True if any single-base run of >= min_run exists."""
    return bool(re.search(r"(A{5,}|T{5,}|G{5,}|C{5,})", seq.upper()))


def _detect_junction_repeat(
    left_seq: str, right_seq: str, check_len: int = 80, threshold: float = 0.85
) -> Tuple[bool, str]:
    """
    Check if the end of left_seq is highly similar to the start of right_seq,
    indicating tandem repeat misassembly risk (e.g., two 2A sequences, two CMV copies).
    """
    if not left_seq or not right_seq:
        return False, ""

    check_len = min(check_len, len(left_seq), len(right_seq))
    if check_len < 20:
        return False, ""

    left_end = left_seq[-check_len:].upper()
    right_start = right_seq[:check_len].upper()

    # Window identity check
    identity = sum(a == b for a, b in zip(left_end, right_start)) / check_len
    if identity >= threshold:
        return True, f"{identity:.0%} identity between end of left module and start of right module"
    return False, ""


def _is_coding_junction(
    left_role: str, right_role: str, left_seq: str, right_seq: str
) -> bool:
    """
    Returns True if the reading frame is expected to be continuous across this junction.
    Heuristic: both sides are coding roles AND left module has no stop codon at its end.
    """
    if left_role not in _CODING_ROLES or right_role not in _CODING_ROLES:
        return False
    # If left module's last codon is a stop codon → frame discontinuous
    if left_seq and len(left_seq) >= 3:
        last_codon = left_seq[-3:].upper()
        if last_codon in ("TAA", "TAG", "TGA"):
            return False
    # If right module starts with ATG → new ORF starts, not a fusion
    if right_seq and right_seq[:3].upper() == "ATG":
        return False
    # Both coding, no stop, no new ATG → likely in-frame fusion
    return left_role in _CODING_ROLES and right_role in _CODING_ROLES


def _is_regulatory_boundary(left_role: str, right_role: str) -> bool:
    """True for promoter→CDS or CDS→polyA junctions."""
    return (
        (left_role in _PROMOTER_ROLES and right_role in _CDS_ROLES)
        or (left_role in _CDS_ROLES and right_role in _TERMINATOR_ROLES)
    )


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ATGCatgc", "TACGtacg")
    return seq.translate(comp)[::-1]


def junction_summary(j: Junction) -> str:
    """One-line human-readable summary of a junction."""
    flags = []
    if j.reading_frame_continuation:
        flags.append("in-frame")
    if not j.scar_allowed:
        flags.append("no-scar")
    if j.regulatory_boundary:
        flags.append("reg-boundary")
    if j.has_tandem_repeat:
        flags.append("⚠ repeat")
    re_names = list(j.internal_restriction_sites.keys())
    if re_names:
        flags.append(f"RE:{','.join(re_names[:3])}")
    flag_str = f" [{', '.join(flags)}]" if flags else ""
    return f"J{j.junction_index}: {j.left_module_name} → {j.right_module_name}{flag_str}"
