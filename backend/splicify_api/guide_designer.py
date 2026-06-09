"""
sgRNA designer for the plasmid viewer.

This is a pure-Python guide designer. Earlier drafts shelled out to R for
crisprVerse scoring; we instead reimplement the parts of crisprVerse that
matter most for *plasmid-scale* guide design (template ≤ 50 kb, no genome
off-target search) so there's nothing to install.

Design contract — for a target region of a plasmid sequence:
1. Scan both strands for a user-specified PAM (IUPAC, e.g. NGG, NAG, TTTV).
2. Extract the protospacer (Cas9 = 20 nt 5' of NGG, Cas12a = 23 nt 3' of TTTV).
3. Score each candidate on a 0-100 composite:
    - GC content: 30 pts when 40-60 %, linear falloff outside
    - Pol III termination: 20 pts iff no internal "TTTT" run
    - Homopolymer: 15 pts (no 5+ runs), 7.5 (4-run), 0 (5+)
    - Seed (3' last 5 nt) GC: up to 20 pts; favours G/C-rich seed
    - Off-targets on plasmid: 15 pts (unique), 7 (1 hit elsewhere), 0 (>1)
4. Sort descending, cap at `max_guides`, return SeqViz-shaped data the
   frontend visualises directly.

Off-target counting compares the spacer + PAM against both strands of the
*same* plasmid only — appropriate for cloning-scale design. Genome-wide
off-target requires bowtie/BWA against a genome index, which the
crisprVerse R stack provides; that is intentionally out of scope here.

Score components are designed to be cheap, deterministic, and reasonably
correlated with crisprVerse's CRISPRko Cas9 ranking. They are NOT a
substitute for Rule Set 2 / CFD when those matter (e.g. genome-wide
panels). The score field reports the integer 0-100 used for ranking, and
`score_components` exposes each contribution so the frontend can show a
breakdown in the gene-card popup.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple

from .doench_score import doench_score_from_plasmid

# IUPAC nucleotide -> set of bases it can match.
_IUPAC: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _pam_to_regex(pam: str) -> re.Pattern:
    parts: List[str] = []
    for ch in pam.upper():
        bases = _IUPAC.get(ch)
        if not bases:
            raise ValueError(f"Invalid PAM character: {ch!r}")
        parts.append(f"[{bases}]" if len(bases) > 1 else bases)
    return re.compile("".join(parts))


def _pam_matches(pam_re: re.Pattern, seq: str, i: int, pam_len: int) -> bool:
    if i < 0 or i + pam_len > len(seq):
        return False
    return bool(pam_re.fullmatch(seq[i : i + pam_len]))


def _revcomp(seq: str) -> str:
    table = str.maketrans("ATGCNatgcn", "TACGNtacgn")
    return seq.translate(table)[::-1]


def _gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(1 for ch in seq if ch in "GCgc") / len(seq)


def _max_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    best, run = 1, 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
            if run > best:
                best = run
        else:
            run = 1
    return best


@dataclass
class GuideHit:
    name: str
    spacer: str
    pam: str
    start: int       # 0-indexed, plasmid coords (forward strand)
    end: int         # half-open
    direction: int   # +1 forward, -1 reverse
    score: float
    score_components: Dict[str, float]
    score_method: str   # "doench2014" or "heuristic"
    context_30mer: str
    gc_fraction: float
    max_homopolymer: int
    n_offtargets: int


def _score_guide(spacer: str, n_offtargets: int) -> Tuple[float, Dict[str, float]]:
    components: Dict[str, float] = {}
    score = 0.0

    gc = _gc_fraction(spacer)
    if 0.40 <= gc <= 0.60:
        gc_pts = 30.0
    else:
        d = max(0.40 - gc, gc - 0.60, 0.0)
        gc_pts = max(0.0, 30.0 - 30.0 * (d / 0.20))
    components["gc"] = round(gc_pts, 2)
    score += gc_pts

    if "TTTT" in spacer.upper():
        tcomp = 0.0
    else:
        tcomp = 20.0
    components["pol3_term"] = tcomp
    score += tcomp

    homo_run = _max_homopolymer(spacer)
    if homo_run >= 5:
        hpoly_pts = 0.0
    elif homo_run == 4:
        hpoly_pts = 7.5
    else:
        hpoly_pts = 15.0
    components["homopolymer"] = hpoly_pts
    score += hpoly_pts

    seed = spacer[-5:] if len(spacer) >= 5 else spacer
    seed_gc = _gc_fraction(seed)
    if seed_gc >= 0.40:
        seed_pts = 20.0 * min(1.0, seed_gc / 0.60)
    else:
        seed_pts = 8.0 * (seed_gc / 0.40) if seed_gc > 0 else 0.0
    components["seed_gc"] = round(seed_pts, 2)
    score += seed_pts

    if n_offtargets == 0:
        ot_pts = 15.0
    elif n_offtargets == 1:
        ot_pts = 7.0
    else:
        ot_pts = 0.0
    components["off_target"] = ot_pts
    score += ot_pts

    return round(score, 2), components


def _scan_pam_hits(
    sequence: str,
    region_start: int,
    region_end: int,
    pam: str,
    guide_length: int,
    pam_position: str,
) -> List[Tuple[int, int, int, str]]:
    """List of (proto_start, proto_end_excl, direction, pam_observed_bases)
    where the protospacer falls fully inside [region_start, region_end)."""
    pam_re = _pam_to_regex(pam)
    pam_len = len(pam)
    seq = sequence
    rc = _revcomp(sequence)
    L = len(sequence)
    hits: List[Tuple[int, int, int, str]] = []

    if pam_position == "3prime":
        # Cas9-style: protospacer 5' of PAM. Forward-strand scan.
        for i in range(region_start, region_end - guide_length + 1):
            j = i + guide_length
            if j + pam_len > L:
                break
            if _pam_matches(pam_re, seq, j, pam_len):
                hits.append((i, j, +1, seq[j : j + pam_len]))
        # Reverse-strand scan: search the reverse-complement, then map back.
        for i in range(L - region_end, L - region_start - guide_length + 1):
            j = i + guide_length
            if j + pam_len > L:
                break
            if _pam_matches(pam_re, rc, j, pam_len):
                fwd_start = L - (j + pam_len) + pam_len
                fwd_end = fwd_start + guide_length
                if fwd_start < region_start or fwd_end > region_end:
                    continue
                hits.append((fwd_start, fwd_end, -1, rc[j : j + pam_len]))
    else:  # 5prime — Cas12a style
        for i in range(region_start - pam_len, region_end - guide_length - pam_len + 1):
            if i < 0 or i + pam_len + guide_length > L:
                continue
            if _pam_matches(pam_re, seq, i, pam_len):
                ps_start = i + pam_len
                if ps_start < region_start or ps_start + guide_length > region_end:
                    continue
                hits.append((ps_start, ps_start + guide_length, +1, seq[i : i + pam_len]))
        for i in range(L - region_end - pam_len, L - region_start - guide_length - pam_len + 1):
            if i < 0 or i + pam_len + guide_length > L:
                continue
            if _pam_matches(pam_re, rc, i, pam_len):
                ps_start_rc = i + pam_len
                fwd_end = L - ps_start_rc
                fwd_start = fwd_end - guide_length
                if fwd_start < region_start or fwd_end > region_end:
                    continue
                hits.append((fwd_start, fwd_end, -1, rc[i : i + pam_len]))

    seen: set = set()
    unique: List[Tuple[int, int, int, str]] = []
    for h in hits:
        key = (h[0], h[1], h[2])
        if key in seen:
            continue
        seen.add(key)
        unique.append(h)
    return unique


def _count_offtargets(sequence: str, spacer: str, pam: str, pam_position: str) -> int:
    pam_re = _pam_to_regex(pam)
    pam_len = len(pam)
    rc_seq = _revcomp(sequence)
    spacer_up = spacer.upper()
    total = 0

    def _count_on(strand: str, query: str) -> int:
        n = 0
        i = 0
        q_len = len(query)
        while True:
            j = strand.find(query, i)
            if j < 0:
                break
            if pam_position == "3prime":
                if _pam_matches(pam_re, strand, j + q_len, pam_len):
                    n += 1
            else:
                if j - pam_len >= 0 and _pam_matches(pam_re, strand, j - pam_len, pam_len):
                    n += 1
            i = j + 1
        return n

    total += _count_on(sequence, spacer_up)
    total += _count_on(rc_seq, spacer_up)
    return max(0, total - 1)


def design_guides(
    sequence: str,
    region_start: int,             # 1-indexed inclusive
    region_end: int,               # 1-indexed inclusive
    pam: str = "NGG",
    guide_length: int = 20,
    pam_position: str = "3prime",  # "3prime" (Cas9) | "5prime" (Cas12a)
    max_guides: int = 50,
    min_score: float = 0.0,
    score_method: str = "doench2014",  # "doench2014" | "heuristic"
) -> Dict:
    seq = re.sub(r"[^ATGCNatgcn]", "", sequence or "").upper()
    if not seq:
        return {"ok": False, "error": "Empty plasmid sequence", "guides": [], "summary": {}}

    L = len(seq)
    rs = max(0, int(region_start) - 1)
    re_ = max(rs + 1, min(L, int(region_end)))
    if re_ - rs < guide_length:
        return {"ok": False,
                "error": f"Region too short: {re_ - rs} bp < guide length {guide_length}",
                "guides": [], "summary": {}}

    raw_hits = _scan_pam_hits(seq, rs, re_, pam, guide_length, pam_position)

    # Doench 2014 Rule Set 1 was calibrated for SpCas9 (NGG PAM, 20-nt
    # protospacer, PAM 3' of protospacer). Use it as the default scorer
    # when those preconditions hold; otherwise fall back to the heuristic.
    cas9_compatible = (
        score_method == "doench2014"
        and pam.upper() == "NGG"
        and guide_length == 20
        and pam_position == "3prime"
    )

    guides: List[GuideHit] = []
    for i, (ps, pe, direction, pam_seq) in enumerate(raw_hits):
        spacer = seq[ps:pe] if direction == +1 else _revcomp(seq[ps:pe])
        n_ot = _count_offtargets(seq, spacer, pam, pam_position)
        homo_run = _max_homopolymer(spacer)
        gc_frac = _gc_fraction(spacer)

        used_method = "heuristic"
        components: Dict[str, float] = {}
        score_val: float
        ctx30 = ""
        if cas9_compatible:
            d_score, d_components, ctx30 = doench_score_from_plasmid(seq, ps, pe, direction)
            if d_score is not None:
                score_val = d_score
                components = dict(d_components)
                # Off-target contribution layered on top of the Doench score:
                # multiplicatively penalise non-unique guides, which the
                # plasmid-scoped Doench score does not see.
                if n_ot >= 2:
                    score_val *= 0.5
                    components["off_target_factor"] = 0.5
                elif n_ot == 1:
                    score_val *= 0.8
                    components["off_target_factor"] = 0.8
                else:
                    components["off_target_factor"] = 1.0
                used_method = "doench2014"

        if used_method == "heuristic":
            score_val, components = _score_guide(spacer, n_ot)
            used_method = "heuristic"

        if score_val < min_score:
            continue

        guides.append(GuideHit(
            name=f"guide_{i + 1:03d}_{spacer[:6]}",
            spacer=spacer,
            pam=pam_seq,
            start=ps,
            end=pe,
            direction=direction,
            score=round(score_val, 2),
            score_components=components,
            score_method=used_method,
            context_30mer=ctx30,
            gc_fraction=round(gc_frac, 3),
            max_homopolymer=homo_run,
            n_offtargets=n_ot,
        ))

    guides.sort(key=lambda g: (-g.score, g.start))
    if max_guides and len(guides) > max_guides:
        guides = guides[:max_guides]

    methods_used = sorted({g.score_method for g in guides})
    return {
        "ok": True,
        "guides": [asdict(g) for g in guides],
        "summary": {
            "n_candidates": len(raw_hits),
            "n_returned": len(guides),
            "region_1based": f"{rs + 1}..{re_}",
            "pam": pam,
            "guide_length": guide_length,
            "pam_position": pam_position,
            "score_method": (
                methods_used[0] if len(methods_used) == 1
                else "+".join(methods_used) if methods_used
                else score_method
            ),
            "score_method_requested": score_method,
            "doench2014_eligible": cas9_compatible,
        },
    }
