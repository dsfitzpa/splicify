"""Score primers for Sanger sequencing quality.

The scoring model is built around the failure modes documented in the
GENEWIZ "Solutions Guide for Sanger Sequencing — Poor Quality" handout,
adapted to the things that are actually inferable from a primer + template
+ target position alone (i.e. design-time, not wet-lab):

  * Tm window — primer Tm should be 50–60 °C (GENEWIZ uses a 50 °C
    annealing temperature for Sanger; primers below that anneal poorly
    and are the dominant cause of "poor priming" results).
  * Self-dimer ΔG — primer3 calcHomodimer; problems arise at ΔG more
    negative than −10 kcal/mol.
  * Hairpin ΔG — primer3 calcHairpin; same threshold.
  * Length / GC content — outside the 18–25 nt and 40–60 % GC bands the
    primer is more likely to mis-prime or anneal off-target.
  * 3' GC clamp — 1–2 G/C in the last 5 bases is ideal; 3+ encourages
    mispriming, 0 destabilises the 3' end.
  * Homopolymer runs — runs ≥ 5 of the same base anywhere in the primer
    cause smearing in the read.
  * Read-window distance — Sanger traces are unreliable for ~50 bp
    immediately after the primer (the "lead-in" Type-2 region) and
    typically fade past ~700 bp. The target feature should land in the
    50–500 bp sweet spot for a clean read.
  * Mispriming — extra exact matches of the primer in the template are
    Type-1 "dominant trace near background" failures.
  * Template secondary structure between the primer and the target — a
    GC-rich stretch or a palindromic hairpin in the read window is the
    classic Type-2 "signal diminishes before pos 500" cause.

`score_sanger_primer(template_id, primers, target_feature_position)`
returns a list of `SangerPrimerScore` dicts, one per input primer, with
`overall_score` in [0, 100] and a `breakdown` of every axis with its
weight, raw value, score, and a human-readable note.

`template_id` may be either a raw template sequence (str) or a dict of
the form `{"sequence": "...", "id": "..."}` so the function plays nicely
with whatever upstream caller wants to pass in.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, Sequence, Union

import primer3

from .utils import normalize_dna, reverse_complement


# ────────────────────────────────────────────────────────────────────────────
# Tunables — chosen to mirror the GENEWIZ doc thresholds. Tweak in one place.
# ────────────────────────────────────────────────────────────────────────────
TM_OPT_RANGE = (50.0, 60.0)        # ideal primer Tm
TM_HARD_LOW = 45.0                  # below this is severe (Sanger anneals at 50)
TM_HARD_HIGH = 65.0                 # above this dimer / mispriming risk grows fast

DIMER_DG_HARD_THRESHOLD = -10_000.0  # cal/mol; <= this is a hard fail
DIMER_DG_SOFT_THRESHOLD = -8_000.0   # cal/mol; linear penalty between soft↔hard

LENGTH_OPT_RANGE = (18, 25)
LENGTH_HARD_RANGE = (15, 35)

GC_OPT_RANGE = (0.40, 0.60)
GC_HARD_RANGE = (0.20, 0.80)

GC_CLAMP_OPT = (1, 2)               # GC count in last 5 bases
HOMOPOLYMER_HARD = 5                # run length >= this is bad

# Read-window scoring (distance from primer 3' end to target_feature_position).
# Anything in the 50–500 bp window is full marks; we taper smoothly outside.
READ_WINDOW_OPT = (50, 500)
READ_WINDOW_USABLE_END = 700
READ_WINDOW_LIMIT = 900             # past this the trace is ~always background

# Template-structure scan (between primer 3' end and read-window end).
STRUCT_GC_WINDOW = 30               # window size for the GC-content scan
STRUCT_GC_HARD = 0.80               # >= this in any window is a hard penalty
STRUCT_PAL_MIN_STEM = 8             # palindrome stem length to flag as a hairpin

# Per-axis weights. Sum normalised at score-aggregation time.
WEIGHTS = {
    "tm": 1.4,
    "self_dimer": 1.2,
    "hairpin": 1.2,
    "length": 0.6,
    "gc_content": 0.6,
    "gc_clamp": 0.5,
    "homopolymer": 0.5,
    "read_window": 1.6,
    "mispriming": 1.0,
    "template_structure": 1.0,
}


# ────────────────────────────────────────────────────────────────────────────
# Public types
# ────────────────────────────────────────────────────────────────────────────
@dataclass
class AxisScore:
    name: str
    score: float            # 0..100
    weight: float
    value: Any              # the raw measurement (Tm, ΔG, count, …)
    note: str               # one-line human-readable rationale
    severity: str           # "ok" | "warn" | "fail"


@dataclass
class SangerPrimerScore:
    primer_index: int
    primer_sequence: str
    direction: str          # "forward" | "reverse"
    binding_start: Optional[int]   # 0-indexed bp on template (5' end of primer as it sits)
    binding_end: Optional[int]     # 0-indexed bp on template (3' end of primer as it sits)
    target_position: int
    distance_to_target: Optional[int]  # bp from 3' end to target_feature_position (signed; − if upstream)
    overall_score: float    # 0..100
    rating: str             # "excellent" | "good" | "marginal" | "poor"
    warnings: List[str]
    breakdown: Dict[str, Dict[str, Any]]


# ────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ────────────────────────────────────────────────────────────────────────────
def _coerce_template(template_id: Union[str, Dict[str, Any]]) -> str:
    """Accept either a raw sequence or a {sequence,id} dict."""
    if isinstance(template_id, str):
        return normalize_dna(template_id)
    if isinstance(template_id, dict):
        seq = template_id.get("sequence") or template_id.get("template") or template_id.get("seq")
        if not seq:
            raise ValueError("template_id dict must include a 'sequence' field")
        return normalize_dna(seq)
    raise TypeError(f"template_id must be str or dict, got {type(template_id).__name__}")


def _coerce_primer(p: Union[str, Dict[str, Any]], idx: int) -> Dict[str, Any]:
    """Normalise a primer-spec into {sequence, direction, binding_start, binding_end, name}."""
    if isinstance(p, str):
        return {"sequence": normalize_dna(p), "direction": "forward",
                "binding_start": None, "binding_end": None, "name": f"primer_{idx}"}
    if not isinstance(p, dict):
        raise TypeError(f"primer at index {idx} must be str or dict")
    seq = p.get("sequence") or p.get("seq") or p.get("primer")
    if not seq:
        raise ValueError(f"primer at index {idx} missing 'sequence'")
    direction = (p.get("direction") or p.get("strand") or "forward").lower()
    if direction in ("rev", "reverse", "-1", "-"):
        direction = "reverse"
    else:
        direction = "forward"
    return {
        "sequence": normalize_dna(seq),
        "direction": direction,
        "binding_start": p.get("binding_start", p.get("start")),
        "binding_end": p.get("binding_end", p.get("end")),
        "name": p.get("name") or f"primer_{idx}",
    }


def _gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GCgc")
    return gc / len(seq)


def _max_homopolymer(seq: str) -> int:
    best = 0
    cur = 0
    last = ""
    for b in seq:
        if b == last:
            cur += 1
        else:
            cur = 1
            last = b
        if cur > best:
            best = cur
    return best


def _gc_clamp_count(seq: str, window: int = 5) -> int:
    tail = seq[-window:]
    return sum(1 for b in tail if b in "GCgc")


def _find_primer_binding(template: str, primer: str, direction: str) -> Optional[tuple[int, int]]:
    """Return (start, end_inclusive) 0-indexed bp on the FORWARD strand where the
    primer anneals. For reverse primers we look for the reverse complement on
    the forward strand. Returns None if no exact match is found.
    """
    if direction == "forward":
        i = template.upper().find(primer.upper())
    else:
        i = template.upper().find(reverse_complement(primer).upper())
    if i < 0:
        return None
    return (i, i + len(primer) - 1)


def _all_exact_hits(template: str, primer: str) -> List[int]:
    """Indices of all exact occurrences of primer (forward) and its reverse
    complement (reverse) on the template."""
    hits: List[int] = []
    t = template.upper()
    for needle in (primer.upper(), reverse_complement(primer).upper()):
        start = 0
        while True:
            idx = t.find(needle, start)
            if idx < 0:
                break
            hits.append(idx)
            start = idx + 1
    return sorted(set(hits))


def _piecewise_score(value: float, *, opt_range: tuple[float, float],
                     hard_range: tuple[float, float]) -> float:
    """Score 100 inside opt_range, linearly decaying to 0 at hard_range edges."""
    lo_opt, hi_opt = opt_range
    lo_hard, hi_hard = hard_range
    if lo_opt <= value <= hi_opt:
        return 100.0
    if value < lo_opt:
        if value <= lo_hard:
            return 0.0
        return 100.0 * (value - lo_hard) / (lo_opt - lo_hard)
    # value > hi_opt
    if value >= hi_hard:
        return 0.0
    return 100.0 * (hi_hard - value) / (hi_hard - hi_opt)


def _dimer_dg(seq: str, kind: str) -> Optional[float]:
    """Return ΔG (cal/mol) for the requested primer3 thermodynamic check."""
    try:
        if kind == "hairpin":
            r = primer3.bindings.calcHairpin(seq)
        else:
            r = primer3.bindings.calcHomodimer(seq)
    except Exception:
        return None
    dg = getattr(r, "dg", None)
    if dg is None or not getattr(r, "structure_found", False):
        return 0.0  # primer3 reports no structure → effectively neutral
    try:
        return float(dg)
    except Exception:
        return None


# ────────────────────────────────────────────────────────────────────────────
# Per-axis scorers
# ────────────────────────────────────────────────────────────────────────────
def _score_tm(seq: str) -> AxisScore:
    try:
        tm = float(primer3.bindings.calcTm(seq))
    except Exception:
        tm = float("nan")
    if math.isnan(tm):
        return AxisScore("tm", 0.0, WEIGHTS["tm"], None,
                         "Tm could not be computed", "fail")
    score = _piecewise_score(tm, opt_range=TM_OPT_RANGE,
                             hard_range=(TM_HARD_LOW, TM_HARD_HIGH))
    severity = "ok" if score >= 80 else "warn" if score >= 40 else "fail"
    note = (f"Tm {tm:.1f} °C — ideal 50–60 °C for Sanger annealing"
            + (" (too low: weak priming)" if tm < TM_OPT_RANGE[0] else
               " (too high: mispriming risk)" if tm > TM_OPT_RANGE[1] else ""))
    return AxisScore("tm", score, WEIGHTS["tm"], round(tm, 2), note, severity)


def _score_dimer(seq: str, axis: str) -> AxisScore:
    dg = _dimer_dg(seq, "homodimer" if axis == "self_dimer" else "hairpin")
    label = "self-dimer" if axis == "self_dimer" else "hairpin"
    if dg is None:
        return AxisScore(axis, 0.0, WEIGHTS[axis], None,
                         f"{label} ΔG could not be computed", "fail")
    if dg >= DIMER_DG_SOFT_THRESHOLD:
        score = 100.0
        severity = "ok"
    elif dg <= DIMER_DG_HARD_THRESHOLD:
        score = 0.0
        severity = "fail"
    else:
        # linear between soft (-8 kcal/mol → 100) and hard (-10 kcal/mol → 0)
        score = 100.0 * (dg - DIMER_DG_HARD_THRESHOLD) / (
            DIMER_DG_SOFT_THRESHOLD - DIMER_DG_HARD_THRESHOLD)
        severity = "warn"
    note = f"{label} ΔG {dg / 1000:.2f} kcal/mol (problem at ≤ −10 kcal/mol)"
    return AxisScore(axis, score, WEIGHTS[axis], round(dg / 1000, 2), note, severity)


def _score_length(seq: str) -> AxisScore:
    L = len(seq)
    score = _piecewise_score(L, opt_range=LENGTH_OPT_RANGE,
                             hard_range=LENGTH_HARD_RANGE)
    severity = "ok" if score >= 80 else "warn" if score >= 40 else "fail"
    note = f"Length {L} nt (ideal 18–25)"
    return AxisScore("length", score, WEIGHTS["length"], L, note, severity)


def _score_gc(seq: str) -> AxisScore:
    gc = _gc_fraction(seq)
    score = _piecewise_score(gc, opt_range=GC_OPT_RANGE, hard_range=GC_HARD_RANGE)
    severity = "ok" if score >= 80 else "warn" if score >= 40 else "fail"
    note = f"GC {gc * 100:.0f}% (ideal 40–60%)"
    return AxisScore("gc_content", score, WEIGHTS["gc_content"],
                     round(gc, 3), note, severity)


def _score_gc_clamp(seq: str) -> AxisScore:
    n = _gc_clamp_count(seq)
    if GC_CLAMP_OPT[0] <= n <= GC_CLAMP_OPT[1]:
        score, severity, note = 100.0, "ok", f"3' clamp = {n} G/C in last 5 nt (ideal 1–2)"
    elif n == 0:
        score, severity, note = 40.0, "warn", "No G/C in last 5 nt — 3' end may be unstable"
    elif n >= 4:
        score, severity, note = 30.0, "warn", f"{n} G/C in last 5 nt — mispriming risk"
    else:  # n == 3
        score, severity, note = 70.0, "warn", "3 G/C in last 5 nt — borderline clamp"
    return AxisScore("gc_clamp", score, WEIGHTS["gc_clamp"], n, note, severity)


def _score_homopolymer(seq: str) -> AxisScore:
    run = _max_homopolymer(seq)
    if run < HOMOPOLYMER_HARD:
        score, severity, note = 100.0, "ok", f"Longest run {run} nt (max safe: 4)"
    elif run == HOMOPOLYMER_HARD:
        score, severity, note = 50.0, "warn", "5-nt run — borderline; may cause stutter"
    else:
        score = max(0.0, 100.0 - (run - HOMOPOLYMER_HARD) * 25.0)
        severity = "fail" if run >= 7 else "warn"
        note = f"{run}-nt homopolymer run — Sanger trace will stutter"
    return AxisScore("homopolymer", score, WEIGHTS["homopolymer"], run, note, severity)


def _score_read_window(distance: Optional[int]) -> AxisScore:
    if distance is None:
        return AxisScore("read_window", 0.0, WEIGHTS["read_window"], None,
                         "Primer does not bind upstream of the target",
                         "fail")
    if distance < 0:
        return AxisScore(
            "read_window", 0.0, WEIGHTS["read_window"], distance,
            f"Target sits {abs(distance)} bp upstream of primer 3' end "
            "— sequencing reads forward, will never reach target",
            "fail")
    if READ_WINDOW_OPT[0] <= distance <= READ_WINDOW_OPT[1]:
        return AxisScore("read_window", 100.0, WEIGHTS["read_window"], distance,
                         f"Target lands {distance} bp downstream — in the "
                         "50–500 bp Sanger sweet spot", "ok")
    if distance < READ_WINDOW_OPT[0]:
        # too close — target sits in the unreliable lead-in region
        score = 100.0 * (distance / READ_WINDOW_OPT[0])
        return AxisScore("read_window", score, WEIGHTS["read_window"], distance,
                         f"Target only {distance} bp from primer 3' end — "
                         "in the unreliable lead-in region (Type-2 risk)",
                         "warn" if score >= 40 else "fail")
    if distance <= READ_WINDOW_USABLE_END:
        # 500..700 — graded decay
        score = 100.0 - 50.0 * (distance - READ_WINDOW_OPT[1]) / (
            READ_WINDOW_USABLE_END - READ_WINDOW_OPT[1])
        return AxisScore("read_window", score, WEIGHTS["read_window"], distance,
                         f"Target {distance} bp downstream — past the "
                         "ideal window, signal will be degrading", "warn")
    if distance <= READ_WINDOW_LIMIT:
        # 700..900 — fast decay to 0
        score = 50.0 * (READ_WINDOW_LIMIT - distance) / (
            READ_WINDOW_LIMIT - READ_WINDOW_USABLE_END)
        return AxisScore("read_window", score, WEIGHTS["read_window"], distance,
                         f"Target {distance} bp downstream — likely past the "
                         "Sanger read length", "fail")
    return AxisScore("read_window", 0.0, WEIGHTS["read_window"], distance,
                     f"Target {distance} bp downstream — beyond Sanger "
                     "read length", "fail")


def _score_mispriming(template: str, primer_seq: str) -> AxisScore:
    hits = _all_exact_hits(template, primer_seq)
    n = len(hits)
    if n <= 1:
        return AxisScore("mispriming", 100.0, WEIGHTS["mispriming"], n,
                         "Single binding site on template", "ok")
    if n == 2:
        return AxisScore("mispriming", 40.0, WEIGHTS["mispriming"], n,
                         "Primer matches template at 2 sites — risk of "
                         "mixed traces", "warn")
    return AxisScore("mispriming", 0.0, WEIGHTS["mispriming"], n,
                     f"Primer matches template at {n} sites — Sanger trace "
                     "will be ambiguous", "fail")


def _score_template_structure(template: str, primer_3p_idx: Optional[int],
                              direction: str, target_position: int) -> AxisScore:
    """Scan the read window for high-GC stretches and palindromic stems."""
    if primer_3p_idx is None:
        return AxisScore("template_structure", 50.0, WEIGHTS["template_structure"],
                         None, "Primer binding not located — structure not "
                         "scanned", "warn")
    # Read window = from primer 3' end forward to ~READ_WINDOW_USABLE_END bp
    # (or up to the target position, whichever is further). For a reverse
    # primer the read travels back toward index 0.
    if direction == "forward":
        win_start = primer_3p_idx + 1
        win_end = min(len(template), max(target_position + 50,
                                         primer_3p_idx + READ_WINDOW_USABLE_END))
        window = template[win_start:win_end]
    else:
        win_end = primer_3p_idx
        win_start = max(0, min(target_position - 50,
                               primer_3p_idx - READ_WINDOW_USABLE_END))
        # Reverse-complement so the "read direction" runs left-to-right.
        window = reverse_complement(template[win_start:win_end])
    if len(window) < STRUCT_GC_WINDOW:
        return AxisScore("template_structure", 100.0,
                         WEIGHTS["template_structure"], None,
                         "Read window too short to scan for secondary structure",
                         "ok")

    # 1) GC-rich windows — anything >= STRUCT_GC_HARD (default 80%) is a flag.
    worst_gc = 0.0
    for i in range(0, len(window) - STRUCT_GC_WINDOW + 1, 5):
        gc = _gc_fraction(window[i:i + STRUCT_GC_WINDOW])
        if gc > worst_gc:
            worst_gc = gc

    # 2) Palindromic stems (a quick proxy for hairpins): scan for any
    #    substring of length >= STRUCT_PAL_MIN_STEM whose reverse complement
    #    appears within ~80 nt of it.
    pal_found = False
    pal_seq = ""
    upper = window.upper()
    for stem_len in (12, 10, STRUCT_PAL_MIN_STEM):
        for i in range(len(upper) - stem_len):
            stem = upper[i:i + stem_len]
            if "N" in stem:
                continue
            rc = reverse_complement(stem)
            j = upper.find(rc, i + stem_len, i + stem_len + 80)
            if j >= 0:
                pal_found = True
                pal_seq = stem
                break
        if pal_found:
            break

    # Aggregate.
    score = 100.0
    sev = "ok"
    notes: List[str] = []
    if worst_gc >= STRUCT_GC_HARD:
        score -= 50.0
        sev = "warn"
        notes.append(f"GC-rich stretch {worst_gc * 100:.0f}% in read window — "
                     "polymerase may stall (Type-2 signal loss)")
    elif worst_gc >= 0.70:
        score -= 25.0
        sev = "warn"
        notes.append(f"Moderately GC-rich stretch {worst_gc * 100:.0f}% in "
                     "read window")
    if pal_found:
        score -= 35.0
        sev = "fail" if score < 40 else "warn"
        notes.append(f"Palindromic stem '{pal_seq}' in read window — likely "
                     "hairpin / Type-2 signal loss")
    if not notes:
        notes.append("Read window has no flagged secondary-structure features")
        sev = "ok"
    score = max(0.0, score)
    return AxisScore("template_structure", score,
                     WEIGHTS["template_structure"],
                     {"max_gc_window": round(worst_gc, 3),
                      "palindrome_stem": pal_seq or None},
                     "; ".join(notes), sev)


# ────────────────────────────────────────────────────────────────────────────
# Public entry point
# ────────────────────────────────────────────────────────────────────────────
def _rate(score: float) -> str:
    if score >= 85:
        return "excellent"
    if score >= 70:
        return "good"
    if score >= 50:
        return "marginal"
    return "poor"


def score_sanger_primer(
    template_id: Union[str, Dict[str, Any]],
    primers: Sequence[Union[str, Dict[str, Any]]],
    target_feature_position: int,
) -> List[Dict[str, Any]]:
    """Score each primer for likelihood of yielding a clean Sanger trace
    over `target_feature_position`.

    Args:
        template_id: the template sequence (str) or a dict with a
            "sequence" key. Anything else with a "sequence"/"seq" field
            also works.
        primers: a list of primer specs. Each entry may be a raw primer
            sequence string OR a dict of the form:
                {"sequence": "...", "direction": "forward"|"reverse",
                 "binding_start": int|None, "binding_end": int|None,
                 "name": "..."}
            If binding positions are omitted the function locates the
            primer by exact match (forward strand for "forward" primers,
            reverse complement for "reverse" primers).
        target_feature_position: 0-indexed bp on the template that should
            land inside the readable Sanger window.

    Returns:
        A list of `SangerPrimerScore`-shaped dicts, one per input primer,
        each with `overall_score` ∈ [0, 100], a 4-tier `rating`, a list of
        `warnings`, and a per-axis `breakdown`.
    """
    template = _coerce_template(template_id)
    if not template:
        raise ValueError("Empty template sequence")
    if not isinstance(target_feature_position, int):
        raise TypeError("target_feature_position must be int (0-indexed bp)")
    if target_feature_position < 0 or target_feature_position >= len(template):
        # Soft-clamp for circular plasmids (caller should pre-mod, but be lenient).
        target_feature_position %= len(template)

    out: List[Dict[str, Any]] = []
    for idx, p in enumerate(primers):
        spec = _coerce_primer(p, idx)
        seq = spec["sequence"]
        direction = spec["direction"]

        # Resolve binding positions on the forward strand if not supplied.
        bs, be = spec.get("binding_start"), spec.get("binding_end")
        if bs is None or be is None:
            located = _find_primer_binding(template, seq, direction)
            if located is not None:
                bs, be = located

        # Distance from primer 3' end → target. Forward primers sequence
        # downstream (target index > 3' end); reverse primers sequence
        # upstream (target index < 3' end).
        if be is not None:
            three_prime = be if direction == "forward" else bs
            if direction == "forward":
                distance = target_feature_position - three_prime
            else:
                distance = three_prime - target_feature_position
        else:
            distance = None

        axes: List[AxisScore] = [
            _score_tm(seq),
            _score_dimer(seq, "self_dimer"),
            _score_dimer(seq, "hairpin"),
            _score_length(seq),
            _score_gc(seq),
            _score_gc_clamp(seq),
            _score_homopolymer(seq),
            _score_read_window(distance),
            _score_mispriming(template, seq),
            _score_template_structure(template,
                                      be if direction == "forward" else bs,
                                      direction, target_feature_position),
        ]

        total_w = sum(a.weight for a in axes)
        overall = sum(a.score * a.weight for a in axes) / total_w if total_w else 0.0
        warnings = [a.note for a in axes if a.severity != "ok"]

        breakdown = {
            a.name: {
                "score": round(a.score, 1),
                "weight": a.weight,
                "value": a.value,
                "note": a.note,
                "severity": a.severity,
            }
            for a in axes
        }

        result = SangerPrimerScore(
            primer_index=idx,
            primer_sequence=seq,
            direction=direction,
            binding_start=bs,
            binding_end=be,
            target_position=target_feature_position,
            distance_to_target=distance,
            overall_score=round(overall, 1),
            rating=_rate(overall),
            warnings=warnings,
            breakdown=breakdown,
        )
        out.append(asdict(result))
    return out
