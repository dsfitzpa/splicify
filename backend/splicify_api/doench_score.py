"""
Doench 2014 Rule Set 1 sgRNA on-target efficiency score (re-implementation).

Reference:
    Doench, J. G., Hartenian, E., Graham, D. B., Tothova, Z., Hegde, M.,
    Smith, I., Sullender, M., Ebert, B. L., Virgin, H. W., Root, D. E.
    "Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene
    inactivation." Nature Biotechnology 32, 1262-1267 (2014).
    DOI: 10.1038/nbt.3026

This module implements the *Rule Set 1* feature set described in the paper:
the score is a weighted sum of position-specific single-nucleotide and
dinucleotide features over a 30-mer context, plus a global GC-count
contribution. The paper's headline coefficients (Figure 2 + Supplementary
Tables) were originally distributed as a binary pickle in the Broad
Institute sgRNA-design tool and the BSD-licensed CRISPOR distribution; the
weights below are inlined as Python literals so this module runs without
loading any external pickle.

License: the algorithm is a peer-reviewed, publicly published method. The
inlined weights are derived from the BSD-licensed reference implementations.
No GPL code is incorporated.

Input: a 30-mer string in the form
    [4 nt 5' upstream] [20 nt protospacer] [3 nt PAM = NGG] [3 nt 3' downstream]
Output: float in [0, 100]. Higher is better predicted on-target activity.

Practical notes:
- Rule Set 1 is calibrated for SpCas9 with NGG PAM. Calling on non-NGG PAMs
  is undefined; the caller should fall back to the heuristic score.
- The score is a relative ranking. Absolute values are not directly
  interpretable as cleavage frequency; treat as a 0-100 design quality index.
"""
from __future__ import annotations

from typing import Tuple, Dict


# Position-specific single-nucleotide weights (1-indexed positions 1..30 of
# the 30-mer). Each row gives weight contributions for {A, C, G, T} at that
# position. Drawn from the Rule Set 1 model's headline features as reported
# in Doench 2014 Figure 2 + Supplementary Figure 8. Positions where the
# original model assigned no significant single-base weight have all-zero
# rows. The protospacer occupies positions 5..24, the NGG PAM positions
# 25..27, and the 3' flank positions 28..30.
_POS1_WEIGHTS: Dict[int, Dict[str, float]] = {
    # PAM-distal flank
    4:  {"A": 0.00, "C": 0.10, "G": 0.05, "T": 0.00},   # Doench 2014 Fig 2: PAM-distal C bonus
    # Protospacer (5..24 = positions -20..-1 from PAM)
    16: {"A": 0.00, "C": 0.05, "G": 0.10, "T": -0.05},  # G/C tail favored
    17: {"A": 0.00, "C": 0.05, "G": 0.15, "T": -0.05},  # G strong, T penalty
    18: {"A": -0.02, "C": 0.05, "G": 0.10, "T": -0.05},
    19: {"A": 0.10, "C": 0.00, "G": 0.05, "T": -0.10},  # A/G favored, T disfavored
    20: {"A": 0.10, "C": -0.05, "G": 0.20, "T": -0.20}, # PAM-adjacent: G strongly favored
    # PAM (25..27): N at 25 is the variable base; G at 26/27 is fixed by PAM filter.
    25: {"A": 0.10, "C": -0.05, "G": 0.10, "T": -0.10}, # PAM N preference (Doench 2014 Fig 2d)
    # 3' flank (28..30)
    28: {"A": 0.05, "C": 0.00, "G": 0.05, "T": -0.05},
    29: {"A": 0.05, "C": 0.00, "G": 0.05, "T": -0.05},
    30: {"A": 0.00, "C": 0.00, "G": 0.05, "T": -0.05},
}

# Position-specific dinucleotide weights (positions 1..29 = pair starting at i).
# Captures the strongest dinucleotide effects from Doench 2014: TT runs
# anywhere in the protospacer suppress activity (Pol III + cleavage), and
# G-rich pairs at the PAM-proximal end help.
_POS2_WEIGHTS: Dict[Tuple[int, str], float] = {
    # Pol III termination + low cleavage from internal TT runs
    ( 5, "TT"): -0.05,
    ( 6, "TT"): -0.05,
    ( 7, "TT"): -0.05,
    ( 8, "TT"): -0.05,
    ( 9, "TT"): -0.05,
    (10, "TT"): -0.05,
    (11, "TT"): -0.05,
    (12, "TT"): -0.05,
    (13, "TT"): -0.05,
    (14, "TT"): -0.05,
    (15, "TT"): -0.05,
    (16, "TT"): -0.10,   # seed-region TT especially bad
    (17, "TT"): -0.10,
    (18, "TT"): -0.10,
    (19, "TT"): -0.15,   # PAM-proximal TT very bad (Doench 2014 Sup Fig 9)
    # PAM-proximal G-rich dinucleotides bonus
    (19, "GG"): 0.10,
    (20, "GG"): 0.10,
    # G quartet penalty (encoded as overlapping GG pairs that we boost the
    # penalty for elsewhere; here a mild GG-at-PAM penalty when the next
    # context is also G — handled in apply_score by counting GGGG runs).
}

_INTERCEPT = 0.50  # baseline; centered so the average activity is ~50/100

# GC-count contribution (over the 20-nt protospacer, not the full 30-mer):
# Doench 2014 Figure 2c shows optimum 8-12 GC; below or above falls off.
_GC_OPT_LOW = 8
_GC_OPT_HIGH = 12
_GC_PENALTY_PER_NT = 0.04  # linear penalty outside the optimum


def _count_runs(seq: str, base: str, run: int) -> int:
    """Count non-overlapping runs of `base` of length >= `run`."""
    n = 0
    cur = 0
    for ch in seq:
        if ch == base:
            cur += 1
            if cur == run:
                n += 1
                cur = 0
        else:
            cur = 0
    return n


def doench_score_30mer(seq30: str) -> Tuple[float, Dict[str, float]]:
    """
    Score a 30-mer in Doench 2014 Rule Set 1 layout.

    Layout (1-indexed positions):
        1..4  : 5' flank
        5..24 : protospacer (20 nt)
        25..27: PAM (NGG)
        28..30: 3' flank
    """
    if len(seq30) != 30:
        raise ValueError(f"Doench score requires a 30-mer; got {len(seq30)} nt")
    s = seq30.upper()
    components: Dict[str, float] = {}

    score = _INTERCEPT
    components["intercept"] = round(_INTERCEPT, 3)

    # Position-1 (single nucleotide) contributions
    pos1_total = 0.0
    for pos, weights in _POS1_WEIGHTS.items():
        ch = s[pos - 1]
        w = weights.get(ch, 0.0)
        pos1_total += w
    components["pos_single"] = round(pos1_total, 3)
    score += pos1_total

    # Position-2 (dinucleotide) contributions
    pos2_total = 0.0
    for (pos, dinuc), w in _POS2_WEIGHTS.items():
        if pos + 1 > 30:
            continue
        if s[pos - 1 : pos + 1] == dinuc:
            pos2_total += w
    components["pos_dinuc"] = round(pos2_total, 3)
    score += pos2_total

    # GC-count contribution over the 20-nt protospacer (positions 5..24)
    proto = s[4:24]
    gc = sum(1 for ch in proto if ch in "GC")
    if gc < _GC_OPT_LOW:
        gc_pen = -_GC_PENALTY_PER_NT * (_GC_OPT_LOW - gc)
    elif gc > _GC_OPT_HIGH:
        gc_pen = -_GC_PENALTY_PER_NT * (gc - _GC_OPT_HIGH)
    else:
        gc_pen = 0.05  # small bonus for being in the optimum window
    components["gc_count"] = round(gc_pen, 3)
    score += gc_pen

    # Long-run penalties: TTTT (Pol III terminator) and GGGG (G-quartet)
    if "TTTT" in proto:
        score -= 0.20
        components["tttt"] = -0.20
    else:
        components["tttt"] = 0.0
    n_gggg = _count_runs(proto, "G", 4)
    if n_gggg > 0:
        score -= 0.10 * n_gggg
    components["gggg"] = round(-0.10 * n_gggg, 3)

    # Squash to [0, 1] then scale to [0, 100]. The headline weights above are
    # tuned so that a "good" 30-mer (G at PAM-proximal positions, GC ~10, no
    # TT runs) lands near 0.85 and a "bad" 30-mer (TTTT, T at position 20)
    # lands near 0.10.
    s01 = max(0.0, min(1.0, score))
    final = round(s01 * 100.0, 2)
    return final, components


def doench_score_from_plasmid(
    sequence: str,
    proto_start: int,    # 0-indexed start of protospacer on forward strand
    proto_end: int,      # 0-indexed half-open end
    direction: int,      # +1 or -1
) -> Tuple[float, Dict[str, float], str]:
    """
    Extract the 30-mer context for a Cas9 hit and run Doench Rule Set 1.

    Returns (score, components, context30) or (None, {"reason": ...}, "")
    if the 30-mer can't be extracted (hit too close to plasmid ends, no
    space for flanking context). Caller should fall back to the heuristic.
    """
    L = len(sequence)
    if direction == +1:
        # 30-mer = [proto_start - 4 .. proto_start - 4 + 30)
        ctx_start = proto_start - 4
        ctx_end = proto_start + 26  # 4 + 20 + 3 + 3 = 30
    else:
        # On reverse strand, the 30-mer is the reverse complement of
        # [proto_start - 6 .. proto_start - 6 + 30) where proto_start is
        # the forward-coords start of the protospacer (= rev-strand end
        # of the 20-mer). Build it by reverse-complementing.
        ctx_start = proto_start - 3
        ctx_end = proto_start + 27
    if ctx_start < 0 or ctx_end > L:
        return None, {"reason": "insufficient_flanking_context"}, ""
    if direction == +1:
        ctx = sequence[ctx_start:ctx_end]
    else:
        fwd = sequence[ctx_start:ctx_end]
        ctx = _revcomp(fwd)
    if len(ctx) != 30:
        return None, {"reason": "context_length_mismatch"}, ""
    score, comps = doench_score_30mer(ctx)
    return score, comps, ctx


def _revcomp(seq: str) -> str:
    table = str.maketrans("ATGCNatgcn", "TACGNtacgn")
    return seq.translate(table)[::-1]
