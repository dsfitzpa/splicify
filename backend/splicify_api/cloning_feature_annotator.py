"""
Cloning feature annotator.

Fast, single-pass scanner that produces cloning-oriented annotations for the
plasmid viewer and backend designers. All emitted features carry
`layer="cloning_feature"` and are hidden by default in the frontend.

Feature families
----------------
- restriction_site_II   : Type II RE hits with strand-resolved top/bottom cut positions
- restriction_site_IIs  : Type IIs RE hits with downstream staggered cut window
- gateway_att           : Gateway att site hits (wraps scan_att_sites)
- primer_design_warning : Regions where PCR primer design is unfeasible

All fast; no secondary-structure ΔG. Inverted repeats use a direct palindromic
walk, not RNA folding.

Consumers
---------
- plannotate_router.py (Step 2.75) merges these into `hierarchical_annotations`
  and also surfaces them as a dedicated `cloning_features` payload with
  `cut_count_per_enzyme` + `non_cutters` for UI filtering.
- restriction_cloning_designer.py, golden_gate_primer_designer.py, and
  gateway_operator.py can call scan_cloning_features() directly to avoid
  re-scanning.
"""
from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .cloning.re_database import RE_DATABASE, RestrictionEnzyme
from .cloning.gateway_sites import scan_att_sites


# ---------------------------------------------------------------------------
# Enzyme tables
# ---------------------------------------------------------------------------

# Mirror of the Type IIs table in golden_gate_primer_designer.py. Keeping a
# local copy avoids a circular import and lets us add cut-offset resolution.
# cut_offset is (top_strand_offset, bottom_strand_offset) past the 3' end of
# the recognition site on the forward strand.
TYPE_IIS_ENZYMES: Dict[str, Dict] = {
    "BsaI":  {"recognition": "GGTCTC", "cut_offset": (1, 5), "overhang_len": 4},
    "BsmBI": {"recognition": "CGTCTC", "cut_offset": (1, 5), "overhang_len": 4},
    "BbsI":  {"recognition": "GAAGAC", "cut_offset": (2, 6), "overhang_len": 4},
    "SapI":  {"recognition": "GCTCTTC", "cut_offset": (1, 4), "overhang_len": 3},
    "BspQI": {"recognition": "GCTCTTC", "cut_offset": (1, 4), "overhang_len": 3},
    # # --- NEB-COMMON-TYPE-IIS-EXPANSION-2026-04-19 ---
    "AarI":   {"recognition": "CACCTGC", "cut_offset": (4, 8),   "overhang_len": 4},
    "Esp3I":  {"recognition": "CGTCTC",  "cut_offset": (1, 5),   "overhang_len": 4},
    "BtgZI":  {"recognition": "GCGATG",  "cut_offset": (10, 14), "overhang_len": 4},
    "BciVI":  {"recognition": "GTATCC",  "cut_offset": (12, 10), "overhang_len": -2},
    "AcuI":   {"recognition": "CTGAAG",  "cut_offset": (16, 14), "overhang_len": -2},
    "BpmI":   {"recognition": "CTGGAG",  "cut_offset": (16, 14), "overhang_len": -2},
    "BsgI":   {"recognition": "GTGCAG",  "cut_offset": (16, 14), "overhang_len": -2},
    "BpuEI":  {"recognition": "CTTGAG",  "cut_offset": (16, 14), "overhang_len": -2},
    "BseRI":  {"recognition": "GAGGAG",  "cut_offset": (10, 8),  "overhang_len": -2},
    "BsrDI":  {"recognition": "GCAATG",  "cut_offset": (2, 0),   "overhang_len": -2},
    "BtsI":   {"recognition": "GCAGTG",  "cut_offset": (2, 0),   "overhang_len": -2},
    "BtsCI":  {"recognition": "GGATG",   "cut_offset": (2, 0),   "overhang_len": -2},
    "BspMI":  {"recognition": "ACCTGC",  "cut_offset": (4, 8),   "overhang_len": 4},
    "EarI":   {"recognition": "CTCTTC",  "cut_offset": (1, 4),   "overhang_len": 3},
    "MlyI":   {"recognition": "GAGTC",   "cut_offset": (5, 5),   "overhang_len": 0},
}


# Default enzyme sets. Callers can override.
DEFAULT_RESTRICTION_II_SET: Tuple[str, ...] = tuple(RE_DATABASE.keys())
DEFAULT_RESTRICTION_IIS_SET: Tuple[str, ...] = ("BsaI", "BsmBI", "BbsI", "SapI")


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class CutProfile:
    """Absolute cut positions for a restriction hit.

    Positions are 0-indexed, between-base coordinates (i.e. cut_top=5 means
    the cut is between bases 4 and 5 of the forward strand).
    """
    cut_top: int
    cut_bottom: int
    overhang_seq: str
    overhang_type: str  # "5prime" | "3prime" | "blunt"
    overhang_len: int


@dataclass
class CloningFeature:
    feature_family: str      # restriction_site_II | restriction_site_IIs | gateway_att | primer_design_warning
    name: str                # display name (enzyme name, att site name, warning label)
    start: int               # recognition / region start (0-indexed, inclusive)
    end: int                 # recognition / region end   (0-indexed, exclusive)
    strand: int = 1          # 1 forward, -1 reverse, 0 n/a
    subtype: str = ""        # category-specific subtype
    cut_profile: Optional[CutProfile] = None
    metadata: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        d = asdict(self)
        if self.cut_profile is None:
            d["cut_profile"] = None
        return d


@dataclass
class ScanResult:
    features: List[CloningFeature]
    cut_count_per_enzyme: Dict[str, int]         # enzyme_name -> number of recognition-site hits (forward + reverse)
    non_cutters: List[str]                        # enzymes in the enabled set with zero hits
    enabled_sets: List[str]

    def to_dict(self) -> Dict:
        return {
            "features": [f.to_dict() for f in self.features],
            "cut_count_per_enzyme": self.cut_count_per_enzyme,
            "non_cutters": self.non_cutters,
            "enabled_sets": self.enabled_sets,
        }


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _iter_matches(pattern: str, sequence: str) -> Iterable[int]:
    if not pattern:
        return
    start = 0
    L = len(pattern)
    while True:
        idx = sequence.find(pattern, start)
        if idx < 0:
            return
        yield idx
        start = idx + 1


# ---------------------------------------------------------------------------
# Type II restriction scanner
# ---------------------------------------------------------------------------

def _scan_type_ii(
    sequence: str,
    enzymes: Sequence[str],
) -> Tuple[List[CloningFeature], Dict[str, int]]:
    seq = sequence.upper()
    L = len(seq)
    features: List[CloningFeature] = []
    counts: Dict[str, int] = defaultdict(int)

    for name in enzymes:
        enz: Optional[RestrictionEnzyme] = RE_DATABASE.get(name)
        if enz is None:
            continue
        recog = enz.recognition_seq.upper()
        rc_recog = _rc(recog)

        # Forward strand
        for pos in _iter_matches(recog, seq):
            cut_top = pos + enz.left_of_cut
            cut_bottom = pos + enz.right_of_cut
            features.append(CloningFeature(
                feature_family="restriction_site_II",
                name=enz.name,
                start=pos,
                end=pos + len(recog),
                strand=1,
                subtype=enz.overhang_type,
                cut_profile=CutProfile(
                    cut_top=cut_top,
                    cut_bottom=cut_bottom,
                    overhang_seq=enz.overhang_seq,
                    overhang_type=enz.overhang_type,
                    overhang_len=len(enz.overhang_seq) if enz.overhang_type != "blunt" else 0,
                ),
                metadata={
                    "recognition_seq": recog,
                    "left_of_cut": enz.left_of_cut,
                    "right_of_cut": enz.right_of_cut,
                    "dam_sensitive": enz.dam_sensitive,
                    "dcm_sensitive": enz.dcm_sensitive,
                    "buffer": enz.buffer,
                },
            ))
            counts[enz.name] += 1

        # Reverse strand (only if palindrome check would miss it)
        if rc_recog != recog:
            for pos in _iter_matches(rc_recog, seq):
                # Mirror the cut offsets to forward-strand coordinates.
                cut_top = pos + (len(recog) - enz.right_of_cut)
                cut_bottom = pos + (len(recog) - enz.left_of_cut)
                features.append(CloningFeature(
                    feature_family="restriction_site_II",
                    name=enz.name,
                    start=pos,
                    end=pos + len(recog),
                    strand=-1,
                    subtype=enz.overhang_type,
                    cut_profile=CutProfile(
                        cut_top=cut_top,
                        cut_bottom=cut_bottom,
                        overhang_seq=enz.overhang_seq,
                        overhang_type=enz.overhang_type,
                        overhang_len=len(enz.overhang_seq) if enz.overhang_type != "blunt" else 0,
                    ),
                    metadata={
                        "recognition_seq": recog,
                        "left_of_cut": enz.left_of_cut,
                        "right_of_cut": enz.right_of_cut,
                        "dam_sensitive": enz.dam_sensitive,
                        "dcm_sensitive": enz.dcm_sensitive,
                        "buffer": enz.buffer,
                    },
                ))
                counts[enz.name] += 1

    # Ensure every enabled enzyme is in the count dict (0 for non-cutters)
    for name in enzymes:
        if name in RE_DATABASE and name not in counts:
            counts[name] = 0

    return features, dict(counts)


# ---------------------------------------------------------------------------
# Type IIs restriction scanner
# ---------------------------------------------------------------------------

def _scan_type_iis(
    sequence: str,
    enzymes: Sequence[str],
) -> Tuple[List[CloningFeature], Dict[str, int]]:
    seq = sequence.upper()
    L = len(seq)
    features: List[CloningFeature] = []
    counts: Dict[str, int] = defaultdict(int)

    for name in enzymes:
        spec = TYPE_IIS_ENZYMES.get(name)
        if spec is None:
            continue
        recog = spec["recognition"].upper()
        rc_recog = _rc(recog)
        top_off, bottom_off = spec["cut_offset"]
        ohlen = spec["overhang_len"]

        # Forward hits: enzyme binds recognition, cuts top/bottom downstream.
        for pos in _iter_matches(recog, seq):
            rec_end = pos + len(recog)
            cut_top = rec_end + top_off
            cut_bottom = rec_end + bottom_off
            if cut_bottom > L:
                continue
            overhang = seq[cut_top:cut_top + ohlen]
            features.append(CloningFeature(
                feature_family="restriction_site_IIs",
                name=name,
                start=pos,
                end=rec_end,
                strand=1,
                subtype="5prime",
                cut_profile=CutProfile(
                    cut_top=cut_top,
                    cut_bottom=cut_bottom,
                    overhang_seq=overhang,
                    overhang_type="5prime",
                    overhang_len=ohlen,
                ),
                metadata={
                    "recognition_seq": recog,
                    "cut_offset": [top_off, bottom_off],
                    "spacer_len": top_off,
                },
            ))
            counts[name] += 1

        # Reverse hits: recognition sits on bottom strand, cuts upstream on forward.
        for pos in _iter_matches(rc_recog, seq):
            cut_bottom = pos - top_off
            cut_top = pos - bottom_off
            if cut_top < 0:
                continue
            overhang = seq[cut_top:cut_top + ohlen]
            features.append(CloningFeature(
                feature_family="restriction_site_IIs",
                name=name,
                start=pos,
                end=pos + len(recog),
                strand=-1,
                subtype="5prime",
                cut_profile=CutProfile(
                    cut_top=cut_top,
                    cut_bottom=cut_bottom,
                    overhang_seq=overhang,
                    overhang_type="5prime",
                    overhang_len=ohlen,
                ),
                metadata={
                    "recognition_seq": recog,
                    "cut_offset": [top_off, bottom_off],
                    "spacer_len": top_off,
                },
            ))
            counts[name] += 1

    for name in enzymes:
        if name in TYPE_IIS_ENZYMES and name not in counts:
            counts[name] = 0

    return features, dict(counts)


# ---------------------------------------------------------------------------
# Gateway att site wrapper
# ---------------------------------------------------------------------------

def _scan_gateway(sequence: str, fuzzy_threshold: int = 0) -> List[CloningFeature]:
    hits = scan_att_sites(sequence, fuzzy_threshold=fuzzy_threshold)
    features: List[CloningFeature] = []
    for h in hits:
        # Core sits inside the site; its exact position depends on site type.
        # attB: LEFT_SHORT(~7) + CORE(7) + RIGHT_SHORT(~7)
        # attP/L/R: long flanks (~80-100) around the core.
        # Find core inside the match by string search (cheap).
        subseq = h.sequence.upper() if h.strand == 1 else _rc(h.sequence.upper())
        core_offset = subseq.find(h.core_sequence)
        if core_offset < 0:
            core_offset = (h.end - h.start - len(h.core_sequence)) // 2
        if h.strand == 1:
            core_start = h.start + core_offset
        else:
            core_start = h.end - core_offset - len(h.core_sequence)
        core_end = core_start + len(h.core_sequence)

        features.append(CloningFeature(
            feature_family="gateway_att",
            name=h.site_type,
            start=h.start,
            end=h.end,
            strand=h.strand,
            subtype=h.site_type[:4],  # attB / attP / attL / attR
            cut_profile=CutProfile(
                # Recombination crossover falls within the 7bp core.
                # Use core midpoint as "cut" for display; both strands same.
                cut_top=(core_start + core_end) // 2,
                cut_bottom=(core_start + core_end) // 2,
                overhang_seq=h.core_sequence,
                overhang_type="blunt",  # recombination, not staggered cut
                overhang_len=len(h.core_sequence),
            ),
            metadata={
                "core_sequence": h.core_sequence,
                "core_start": core_start,
                "core_end": core_end,
                "match_quality": h.match_quality,
                "site_number": h.site_type[-1] if h.site_type[-1].isdigit() else "",
            },
        ))
    return features


# ---------------------------------------------------------------------------
# PCR-design feasibility scanner (fast only; no secondary-structure ΔG)
# ---------------------------------------------------------------------------

@dataclass
class PcrWarningThresholds:
    gc_window: int = 25
    gc_low: float = 0.30
    gc_high: float = 0.75
    homopolymer_min_len: int = 7        # >6
    microsat_unit_min: int = 1
    microsat_unit_max: int = 6
    microsat_min_tile: int = 13         # >12
    direct_repeat_min_len: int = 17     # ≥17
    inverted_repeat_min_stem: int = 6
    inverted_repeat_max_loop: int = 7   # <8
    palindrome_min_len: int = 10
    blast_seed_k: int = 11
    blast_min_hit_len: int = 20
    blast_min_identity: float = 0.80


def _scan_homopolymers(seq: str, min_len: int) -> List[Tuple[int, int, str]]:
    out = []
    if not seq:
        return out
    i = 0
    L = len(seq)
    while i < L:
        j = i + 1
        while j < L and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            out.append((i, j, seq[i]))
        i = j
    return out


def _scan_microsatellites(
    seq: str,
    unit_min: int,
    unit_max: int,
    min_tile: int,
) -> List[Tuple[int, int, str]]:
    """Find regions tiled by a repeat unit of length unit_min..unit_max
    covering ≥ min_tile bp. Greedy per unit-length; merges overlapping hits
    across unit lengths by keeping the longer span.
    """
    L = len(seq)
    spans: List[Tuple[int, int, str]] = []
    for k in range(unit_min, unit_max + 1):
        if L < 2 * k:
            continue
        i = 0
        while i + k <= L:
            unit = seq[i:i + k]
            if "N" in unit or len(set(unit)) == 1 and k > 1:
                # skip pure homopolymer at k>1 (already covered by _scan_homopolymers for k=1)
                i += 1
                continue
            j = i + k
            while j + k <= L and seq[j:j + k] == unit:
                j += k
            span_len = j - i
            if span_len >= min_tile and span_len >= 2 * k:
                spans.append((i, j, unit))
                i = j
            else:
                i += 1
    # Merge overlapping spans; keep longest overlap representative.
    if not spans:
        return []
    spans.sort(key=lambda s: (s[0], -(s[1] - s[0])))
    merged: List[Tuple[int, int, str]] = []
    for s in spans:
        if merged and s[0] < merged[-1][1]:
            # overlap: keep the longer
            prev = merged[-1]
            if (s[1] - s[0]) > (prev[1] - prev[0]):
                merged[-1] = s
        else:
            merged.append(s)
    return merged


def _scan_gc_extremes(
    seq: str,
    window: int,
    low: float,
    high: float,
) -> List[Tuple[int, int, float, str]]:
    """Sliding GC window. Returns (start, end, gc_fraction, tag) for maximal
    runs where GC is outside [low, high]."""
    L = len(seq)
    if L < window:
        return []
    # Rolling GC count.
    gc = sum(1 for c in seq[:window] if c in "GCgc")
    flags_low: List[int] = []
    flags_high: List[int] = []
    for i in range(L - window + 1):
        frac = gc / window
        if frac < low:
            flags_low.append(i)
        elif frac > high:
            flags_high.append(i)
        if i + window < L:
            out_c = seq[i]
            in_c = seq[i + window]
            if out_c in "GCgc":
                gc -= 1
            if in_c in "GCgc":
                gc += 1

    def _collapse(idxs: List[int], tag: str) -> List[Tuple[int, int, float, str]]:
        if not idxs:
            return []
        spans: List[Tuple[int, int, float, str]] = []
        start = idxs[0]
        prev = idxs[0]
        for i in idxs[1:]:
            if i == prev + 1:
                prev = i
            else:
                end = prev + window
                sub = seq[start:end]
                frac = sum(1 for c in sub if c in "GCgc") / max(1, len(sub))
                spans.append((start, end, frac, tag))
                start = prev = i
        end = prev + window
        sub = seq[start:end]
        frac = sum(1 for c in sub if c in "GCgc") / max(1, len(sub))
        spans.append((start, end, frac, tag))
        return spans

    return _collapse(flags_low, "low_gc") + _collapse(flags_high, "high_gc")


def _scan_direct_repeats(seq: str, min_len: int) -> List[Tuple[int, int, int, int]]:
    """Exact direct repeats of length ≥ min_len. Returns (a_start, a_end,
    b_start, b_end) with a_start < b_start. Uses k-mer index for speed; O(n)
    memory, O(n) expected time for typical plasmid sequences."""
    L = len(seq)
    k = min_len
    if L < 2 * k:
        return []
    index: Dict[str, List[int]] = defaultdict(list)
    for i in range(L - k + 1):
        kmer = seq[i:i + k]
        if "N" in kmer:
            continue
        index[kmer].append(i)
    hits: List[Tuple[int, int, int, int]] = []
    seen_pairs: set = set()
    for kmer, positions in index.items():
        if len(positions) < 2:
            continue
        for a in range(len(positions)):
            for b in range(a + 1, len(positions)):
                pa, pb = positions[a], positions[b]
                if pb - pa < k:
                    continue
                # Extend match rightward.
                ext = 0
                while (pa + k + ext < L and pb + k + ext < L
                       and seq[pa + k + ext] == seq[pb + k + ext]):
                    ext += 1
                span = k + ext
                key = (pa, pb, span)
                if key in seen_pairs:
                    continue
                # Skip sub-hits that will be found by their leftmost anchor.
                if pa > 0 and pb > 0 and seq[pa - 1] == seq[pb - 1]:
                    continue
                seen_pairs.add(key)
                hits.append((pa, pa + span, pb, pb + span))
    return hits


def _scan_blast_like_homology(
    seq: str,
    seed_k: int,
    min_hit_len: int,
    min_identity: float,
) -> List[Tuple[int, int, int, int, float]]:
    """BLAST-style seed-and-extend within a single sequence. Returns (a_start,
    a_end, b_start, b_end, identity) for ungapped matches ≥ min_hit_len with
    identity ≥ min_identity. The forward-vs-forward matches that are already
    exact direct repeats (identity=1.0) are subsumed; this pass looks for
    approximate matches."""
    L = len(seq)
    if L < 2 * min_hit_len:
        return []
    index: Dict[str, List[int]] = defaultdict(list)
    for i in range(L - seed_k + 1):
        kmer = seq[i:i + seed_k]
        if "N" in kmer:
            continue
        index[kmer].append(i)

    hits: List[Tuple[int, int, int, int, float]] = []
    emitted: set = set()
    for kmer, positions in index.items():
        if len(positions) < 2:
            continue
        for a in range(len(positions)):
            for b in range(a + 1, len(positions)):
                pa, pb = positions[a], positions[b]
                if pb - pa < min_hit_len:
                    continue
                # Extend left
                la = pa
                lb = pb
                left_mm = 0
                while la > 0 and lb > 0:
                    if seq[la - 1] == seq[lb - 1]:
                        la -= 1
                        lb -= 1
                    else:
                        # Allow one mismatch while we're still inside seed extension,
                        # but stop once we've accumulated enough to drop identity.
                        break
                # Extend right with a bounded mismatch budget.
                ra = pa + seed_k
                rb = pb + seed_k
                mm = 0
                total = seed_k + (pa - la)
                while ra < L and rb < L and (ra - la) < (L - pb):
                    if seq[ra] == seq[rb]:
                        ra += 1
                        rb += 1
                        total += 1
                    else:
                        # Allow a mismatch if running identity stays above min_identity
                        matches = total - mm
                        new_mm = mm + 1
                        new_total = total + 1
                        if (matches / new_total) >= min_identity:
                            mm = new_mm
                            total = new_total
                            ra += 1
                            rb += 1
                        else:
                            break
                length = ra - la
                if length < min_hit_len:
                    continue
                a_seg = seq[la:ra]
                b_seg = seq[lb:rb]
                match_count = sum(1 for x, y in zip(a_seg, b_seg) if x == y)
                identity = match_count / length
                if identity < min_identity:
                    continue
                # Skip exact full-length matches (handled by direct-repeat scan)
                if identity == 1.0:
                    continue
                key = (la, ra, lb, rb)
                if key in emitted:
                    continue
                emitted.add(key)
                hits.append((la, ra, lb, rb, identity))
    return hits


def _scan_palindromes_and_inverted_repeats(
    seq: str,
    palindrome_min: int,
    ir_min_stem: int,
    ir_max_loop: int,
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int, int, int]]]:
    """Walk every center (and every inter-base center) to find maximal
    palindromes (zero-loop) and inverted repeats (nonzero loop). Returns:
      palindromes: [(start, end)]
      inverted_repeats: [(stem_a_start, stem_a_end, stem_b_start, stem_b_end)]
    """
    L = len(seq)
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    palindromes: List[Tuple[int, int]] = []
    irs: List[Tuple[int, int, int, int]] = []

    # Zero-loop (even-length palindromes): center between bases i and i+1.
    for i in range(L - 1):
        a = i
        b = i + 1
        while a >= 0 and b < L and comp.get(seq[a], "") == seq[b]:
            a -= 1
            b += 1
        stem_len = (i - a)
        total_len = b - a - 1
        if total_len >= palindrome_min:
            palindromes.append((a + 1, b))

    # Odd-length palindromes centered on a base: include only if center base is
    # its own complement — impossible for standard ACGT, so skip.

    # Inverted repeats with loop > 0. Center loop between two positions
    # [loop_start, loop_end). For each possible loop size Lp in [1, ir_max_loop]:
    for loop in range(1, ir_max_loop + 1):
        for center in range(L - loop):
            # Loop region seq[center : center+loop]
            a = center - 1
            b = center + loop
            while a >= 0 and b < L and comp.get(seq[a], "") == seq[b]:
                a -= 1
                b += 1
            stem_len = (center - 1) - a
            if stem_len >= ir_min_stem:
                stem_a = (a + 1, center)
                stem_b = (center + loop, b)
                irs.append((stem_a[0], stem_a[1], stem_b[0], stem_b[1]))

    # Dedup palindromes by span; keep longest.
    palindromes = sorted(set(palindromes), key=lambda p: (p[0], -(p[1] - p[0])))
    deduped_p: List[Tuple[int, int]] = []
    for p in palindromes:
        if deduped_p and p[0] < deduped_p[-1][1]:
            if (p[1] - p[0]) > (deduped_p[-1][1] - deduped_p[-1][0]):
                deduped_p[-1] = p
        else:
            deduped_p.append(p)

    # Dedup IRs by (stem_a_start, stem_b_end) keeping the longest.
    irs_sorted = sorted(set(irs), key=lambda r: (r[0], -(r[3] - r[0])))
    deduped_irs: List[Tuple[int, int, int, int]] = []
    for r in irs_sorted:
        if deduped_irs and r[0] == deduped_irs[-1][0] and r[3] <= deduped_irs[-1][3]:
            continue
        deduped_irs.append(r)

    return deduped_p, deduped_irs


def _scan_pcr_warnings(
    sequence: str,
    thresholds: Optional[PcrWarningThresholds] = None,
) -> List[CloningFeature]:
    t = thresholds or PcrWarningThresholds()
    seq = sequence.upper()
    warnings: List[CloningFeature] = []

    # Homopolymers
    for (s, e, base) in _scan_homopolymers(seq, t.homopolymer_min_len):
        warnings.append(CloningFeature(
            feature_family="primer_design_warning",
            name=f"Homopolymer {base}×{e - s}",
            start=s, end=e, strand=0,
            subtype="homopolymer",
            metadata={"base": base, "length": e - s},
        ))

    # Microsatellite / short-unit repeat tiling
    for (s, e, unit) in _scan_microsatellites(
        seq, t.microsat_unit_min, t.microsat_unit_max, t.microsat_min_tile,
    ):
        # Skip k=1 (already reported as homopolymer).
        if len(unit) == 1:
            continue
        warnings.append(CloningFeature(
            feature_family="primer_design_warning",
            name=f"Repeat ({unit})×{(e - s) // len(unit)}",
            start=s, end=e, strand=0,
            subtype="microsatellite",
            metadata={"unit": unit, "unit_len": len(unit),
                      "copies": (e - s) // len(unit), "span": e - s},
        ))

    # GC extremes
    for (s, e, frac, tag) in _scan_gc_extremes(seq, t.gc_window, t.gc_low, t.gc_high):
        warnings.append(CloningFeature(
            feature_family="primer_design_warning",
            name=f"{'Low' if tag == 'low_gc' else 'High'} GC ({frac*100:.0f}%)",
            start=s, end=e, strand=0,
            subtype=tag,
            metadata={"gc_fraction": round(frac, 3), "window": t.gc_window},
        ))

    # Exact direct repeats ≥ min_len
    for (a_s, a_e, b_s, b_e) in _scan_direct_repeats(seq, t.direct_repeat_min_len):
        warnings.append(CloningFeature(
            feature_family="primer_design_warning",
            name=f"Direct repeat {a_e - a_s} bp",
            start=a_s, end=a_e, strand=0,
            subtype="direct_repeat",
            metadata={"partner_start": b_s, "partner_end": b_e,
                      "length": a_e - a_s, "identity": 1.0},
        ))
        warnings.append(CloningFeature(
            feature_family="primer_design_warning",
            name=f"Direct repeat {b_e - b_s} bp",
            start=b_s, end=b_e, strand=0,
            subtype="direct_repeat",
            metadata={"partner_start": a_s, "partner_end": a_e,
                      "length": b_e - b_s, "identity": 1.0},
        ))

    # BLAST-style approximate self-homology (80%+ over ≥20 bp)
    for (a_s, a_e, b_s, b_e, ident) in _scan_blast_like_homology(
        seq, t.blast_seed_k, t.blast_min_hit_len, t.blast_min_identity,
    ):
        warnings.append(CloningFeature(
            feature_family="primer_design_warning",
            name=f"Self-homology {a_e - a_s} bp ({ident*100:.0f}%)",
            start=a_s, end=a_e, strand=0,
            subtype="self_homology",
            metadata={"partner_start": b_s, "partner_end": b_e,
                      "length": a_e - a_s, "identity": round(ident, 3)},
        ))
        warnings.append(CloningFeature(
            feature_family="primer_design_warning",
            name=f"Self-homology {b_e - b_s} bp ({ident*100:.0f}%)",
            start=b_s, end=b_e, strand=0,
            subtype="self_homology",
            metadata={"partner_start": a_s, "partner_end": a_e,
                      "length": b_e - b_s, "identity": round(ident, 3)},
        ))

    # Palindromes and inverted repeats
    palindromes, irs = _scan_palindromes_and_inverted_repeats(
        seq, t.palindrome_min_len, t.inverted_repeat_min_stem, t.inverted_repeat_max_loop,
    )
    for (s, e) in palindromes:
        warnings.append(CloningFeature(
            feature_family="primer_design_warning",
            name=f"Palindrome {e - s} bp",
            start=s, end=e, strand=0,
            subtype="palindrome",
            metadata={"length": e - s},
        ))
    for (a_s, a_e, b_s, b_e) in irs:
        loop_len = b_s - a_e
        stem_len = a_e - a_s
        warnings.append(CloningFeature(
            feature_family="primer_design_warning",
            name=f"Inverted repeat stem {stem_len} / loop {loop_len}",
            start=a_s, end=b_e, strand=0,
            subtype="inverted_repeat",
            metadata={"stem_a": [a_s, a_e], "stem_b": [b_s, b_e],
                      "stem_len": stem_len, "loop_len": loop_len},
        ))

    return warnings


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def scan_cloning_features(
    sequence: str,
    *,
    enabled_sets: Optional[Sequence[str]] = None,
    type_ii_enzymes: Optional[Sequence[str]] = None,
    type_iis_enzymes: Optional[Sequence[str]] = None,
    gateway_fuzzy_threshold: int = 0,
    pcr_thresholds: Optional[PcrWarningThresholds] = None,
) -> ScanResult:
    """Run the unified cloning-feature scanner.

    Args:
        sequence: DNA sequence (linear view; caller handles circular wrap if needed).
        enabled_sets: subset of {"restriction_II", "restriction_IIs", "gateway", "pcr"}.
            Default = all four.
        type_ii_enzymes: restrict Type II scan to this subset.
        type_iis_enzymes: restrict Type IIs scan to this subset.
        gateway_fuzzy_threshold: passed through to scan_att_sites.
        pcr_thresholds: override for PCR feasibility thresholds.
    """
    enabled = list(enabled_sets) if enabled_sets is not None else [
        "restriction_II", "restriction_IIs", "gateway", "pcr",
    ]

    features: List[CloningFeature] = []
    cut_counts: Dict[str, int] = {}
    non_cutters: List[str] = []

    if "restriction_II" in enabled:
        enz = type_ii_enzymes or DEFAULT_RESTRICTION_II_SET
        f, c = _scan_type_ii(sequence, enz)
        features.extend(f)
        cut_counts.update(c)
        non_cutters.extend([n for n in enz if c.get(n, 0) == 0])

    if "restriction_IIs" in enabled:
        enz = type_iis_enzymes or DEFAULT_RESTRICTION_IIS_SET
        f, c = _scan_type_iis(sequence, enz)
        features.extend(f)
        cut_counts.update(c)
        non_cutters.extend([n for n in enz if c.get(n, 0) == 0])

    if "gateway" in enabled:
        features.extend(_scan_gateway(sequence, fuzzy_threshold=gateway_fuzzy_threshold))

    if "pcr" in enabled:
        features.extend(_scan_pcr_warnings(sequence, pcr_thresholds))

    return ScanResult(
        features=features,
        cut_count_per_enzyme=cut_counts,
        non_cutters=sorted(set(non_cutters)),
        enabled_sets=enabled,
    )


# ---------------------------------------------------------------------------
# Annotation conversion for hierarchical_annotations merge
# ---------------------------------------------------------------------------

_FAMILY_COLORS = {
    "restriction_site_II": "#C2185B",     # magenta (matches EcoRI bracket reference)
    "restriction_site_IIs": "#FF6F00",    # amber
    "gateway_att": "#6A1B9A",             # purple
    "primer_design_warning": "#D32F2F",   # red
}


def cloning_features_to_hierarchical(
    features: Sequence[CloningFeature],
) -> List[Dict]:
    """Convert CloningFeature list into hierarchical_annotation dicts ready
    to merge into the annotation response. Sets layer='cloning_feature' so
    the viewer can hide them by default."""
    out: List[Dict] = []
    for f in features:
        direction = f.strand if f.strand in (-1, 0, 1) else 0
        cut_profile = None
        if f.cut_profile is not None:
            cut_profile = asdict(f.cut_profile)
        out.append({
            "name": f.name,
            "start": f.start,
            "end": f.end,
            "direction": direction,
            "color": _FAMILY_COLORS.get(f.feature_family, "#607D8B"),
            "layer": "cloning_feature",
            "feature_family": f.feature_family,
            "subtype": f.subtype,
            "cut_profile": cut_profile,
            "metadata": f.metadata,
            "source": "cloning_feature_annotator",
        })
    return out
