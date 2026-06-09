from __future__ import annotations

from typing import Optional, Dict, Any, List, Tuple
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field, ConfigDict
import primer3

from .utils import normalize_dna, ensure_session_id, is_valid_dna, reverse_complement
from .sanger_scoring import score_sanger_primer

router = APIRouter(tags=["pcr"])


# Application presets ──────────────────────────────────────────────────────
# Adapters are prepended to the 5' end of the corresponding designed primer.
# Tm is intentionally calculated for the annealing sequence only (primer3
# already designs against the template, so the Tm it reports IS the
# annealing Tm — no recomputation needed when adapters are added).
APPLICATIONS = ("fragment", "sanger", "illumina")
ILLUMINA_FWD_ADAPTER = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
ILLUMINA_REV_ADAPTER = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
ILLUMINA_MAX_AMPLICON_BP = 600


class PrimerRequest(BaseModel):
    """
    PCR primer design request.

    Notes:
    - fragments_in: your template sequence (string)
    - excluded region: use excluded_start/excluded_length OR legacy target_start/target_length
      to force primers to flank that excluded interval.
    """
    model_config = ConfigDict(populate_by_name=True, extra="ignore")

    fragments_in: str

    # Excluded region (preferred names)
    excluded_start: Optional[int] = None
    excluded_length: Optional[int] = None
    excluded_end: Optional[int] = None  # NEW (accept end-style)


    # Legacy names (keep compatibility with existing n8n/router)
    target_start: Optional[int] = Field(default=None, description="Legacy alias for excluded_start")
    target_length: Optional[int] = Field(default=None, description="Legacy alias for excluded_length")

    # Amplicon constraints
    product_size_min: Optional[int] = 100
    product_size_max: Optional[int] = 300

    # Tm tuning (allow nulls -> FastAPI defaults applied here)
    primer_min_tm: Optional[float] = None
    primer_opt_tm: Optional[float] = None
    primer_max_tm: Optional[float] = None

    # Primer length constraints (optional; if null, primer3 defaults apply)
    primer_min_size: Optional[int] = None
    primer_opt_size: Optional[int] = None
    primer_max_size: Optional[int] = None

    # How many candidate pairs to retrieve from primer3 before selecting best one
    num_return: int = 5

    # Downstream application — drives adapter injection / size guard.
    application: Optional[str] = "fragment"

    session_id: Optional[str] = Field(default=None, alias="sessionId")


def _coalesce_excluded(req: PrimerRequest) -> Optional[Tuple[int, int]]:
    """Return (start, length) excluded region if provided, else None."""
    start = req.excluded_start if req.excluded_start is not None else req.target_start
    length = req.excluded_length if req.excluded_length is not None else req.target_length
    if length is None and req.excluded_end is not None and start is not None:
        length = int(req.excluded_end) - int(start)

    
    if start is None or length is None:
        return None
    try:
        s = int(start)
        L = int(length)
    except Exception:
        return None
    if s < 0 or L <= 0:
        return None
    return (s, L)


def _primer_left_pos(result: Dict[str, Any], i: int) -> Optional[Tuple[int, int]]:
    # PRIMER_LEFT_i is [start, length]
    v = result.get(f"PRIMER_LEFT_{i}")
    if isinstance(v, (list, tuple)) and len(v) >= 2:
        try:
            return (int(v[0]), int(v[1]))
        except Exception:
            return None
    return None


def _primer_right_pos(result: Dict[str, Any], i: int) -> Optional[Tuple[int, int]]:
    # PRIMER_RIGHT_i is [start, length] where start is the 3'-most base index
    v = result.get(f"PRIMER_RIGHT_{i}")
    if isinstance(v, (list, tuple)) and len(v) >= 2:
        try:
            return (int(v[0]), int(v[1]))
        except Exception:
            return None
    return None


def _pair_flanks_excluded(
    left_pos: Tuple[int, int],
    right_pos: Tuple[int, int],
    excluded: Tuple[int, int],
) -> bool:
    """
    Ensure:
      left primer binds entirely to the LEFT of excluded region
      right primer binds entirely to the RIGHT of excluded region

    excluded is (start, length), interval = [start, start+length)
    left primer spans [Lstart, Lstart+Llen)
    right primer spans [Rstart-Rlen+1, Rstart+1) because primer3 uses 3'-end index
    """
    ex_start, ex_len = excluded
    ex_end = ex_start + ex_len

    Lstart, Llen = left_pos
    Lend = Lstart + Llen

    R3, Rlen = right_pos
    Rstart = (R3 - Rlen + 1)
    Rend = R3 + 1

    # left entirely before excluded start
    if not (Lend <= ex_start):
        return False

    # right entirely after excluded end
    if not (Rstart >= ex_end):
        return False

    # sanity: left before right
    if not (Lstart < Rstart):
        return False

    return True


def _mispriming_sites(template: str, primer: str) -> List[int]:
    """
    Very simple mispriming indicator:
    returns all start indices where primer OR its reverse-complement appears in template.
    (Exact match only. You can upgrade this later to allow mismatches.)
    """
    t = template
    p = primer
    rc = reverse_complement(primer)

    hits: List[int] = []
    for q in (p, rc):
        start = 0
        while True:
            idx = t.find(q, start)
            if idx == -1:
                break
            hits.append(idx)
            start = idx + 1
    hits = sorted(set(hits))
    return hits


def _thermo_scores(primer: str) -> Dict[str, Any]:
    """
    Provide hairpin / self-dimer / end-self-dimer approximations.
    primer3-py returns objects; we normalize to floats where possible.
    """
    out: Dict[str, Any] = {}
    try:
        hd = primer3.bindings.calcHairpin(primer)
        out["hairpin_th"] = float(getattr(hd, "tm", None)) if getattr(hd, "tm", None) is not None else None
    except Exception:
        out["hairpin_th"] = None

    try:
        sd = primer3.bindings.calcHomodimer(primer)
        out["any_th"] = float(getattr(sd, "tm", None)) if getattr(sd, "tm", None) is not None else None
    except Exception:
        out["any_th"] = None

    try:
        esd = primer3.bindings.calcEndStability(primer, primer)
        out["3p_th"] = float(getattr(esd, "tm", None)) if getattr(esd, "tm", None) is not None else None
    except Exception:
        out["3p_th"] = None

    return out


def _sanger_scores_for_pair(
    template: str,
    left_seq: str,
    right_seq: str,
    left_pos,
    right_pos,
    excluded,
):
    """Run score_sanger_primer against the AMPLICON the primer pair
    produces — that is what the downstream Sanger reaction actually uses
    as template after PCR. We slice the amplicon out of the plasmid,
    translate the primer binding positions and the target into amplicon
    coordinates, then call the scorer.

    For Sanger sequencing of a PCR product, mispriming and template
    structure that lie outside the amplicon are irrelevant — only the
    amplicon matters. Scoring the full plasmid would unfairly penalise
    primers whose plasmid hits don't carry through to the amplicon.

    Forward primer reads downstream toward the target; reverse primer
    reads upstream toward it. The target defaults to the centre of the
    excluded region (the user-selected feature ± padding) and falls
    back to the centre of the amplicon when no excluded region is set.
    """
    if not left_pos or not right_pos or not left_seq or not right_seq:
        return []

    left_start = int(left_pos[0])
    left_len = int(left_pos[1])
    rp_3p = int(right_pos[0])
    rp_len = int(right_pos[1])
    right_end_excl = rp_3p + 1  # half-open

    # Linear amplicon slice. Origin-wrapping circular amplicons are rare
    # in the Sanger flow (the design region is constrained); if we ever
    # see one we fall back to the full template so the score is still
    # informative.
    if left_start < right_end_excl <= len(template):
        amplicon = template[left_start:right_end_excl]
        amp_offset = left_start
    else:
        amplicon = template
        amp_offset = 0

    # Target in template coords.
    if excluded:
        target_template = excluded[0] + excluded[1] // 2
    else:
        target_template = (left_start + right_end_excl) // 2
    target_amp = target_template - amp_offset
    if not (0 <= target_amp < len(amplicon)):
        target_amp = len(amplicon) // 2

    primers = [
        {
            "name": "forward",
            "sequence": left_seq,
            "direction": "forward",
            "binding_start": 0,
            "binding_end": left_len - 1,
        },
        {
            "name": "reverse",
            "sequence": right_seq,
            "direction": "reverse",
            "binding_start": len(amplicon) - rp_len,
            "binding_end": len(amplicon) - 1,
        },
    ]
    try:
        return score_sanger_primer(amplicon, primers, int(target_amp))
    except Exception as e:
        return [{"error": f"sanger scoring failed: {e}"}]


class SangerScoreRequest(BaseModel):
    """Direct entry point for the score_sanger_primer scorer. Accepts a
    template + a list of primers + a target feature position (1-indexed
    or 0-indexed via `target_is_one_indexed`)."""

    model_config = ConfigDict(populate_by_name=True, extra="ignore")

    template: str
    primers: List[Dict[str, Any]]
    target_feature_position: int
    target_is_one_indexed: bool = True


@router.post("/score-sanger-primers")
def score_sanger_primers_endpoint(req: SangerScoreRequest):
    template = normalize_dna(req.template)
    if not is_valid_dna(template):
        raise HTTPException(status_code=400, detail="Invalid DNA template")
    pos = int(req.target_feature_position)
    if req.target_is_one_indexed:
        pos -= 1
    if pos < 0 or pos >= len(template):
        pos %= len(template)
    try:
        scores = score_sanger_primer(template, req.primers, pos)
    except (ValueError, TypeError) as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {"scores": scores, "target_position_zero_indexed": pos,
            "template_length": len(template)}


@router.post("/design-primers")
def design_primers(req: PrimerRequest):
    template = normalize_dna(req.fragments_in)
    if not is_valid_dna(template):
        raise HTTPException(status_code=400, detail="Invalid DNA fragments_in")

    excluded = _coalesce_excluded(req)
    if excluded:
        ex_start, ex_len = excluded
        if ex_start + ex_len > len(template):
            raise HTTPException(status_code=400, detail="Excluded region exceeds template length")

    # If excluded region exists, force product_size_min >= excluded length (helps prevent primers inside)
    product_min = int(req.product_size_min) if req.product_size_min is not None else 100
    product_max = int(req.product_size_max) if req.product_size_max is not None else 300
    if excluded:
        product_min = max(product_min, excluded[1])

    if product_min <= 0 or product_max <= 0 or product_min > product_max:
        raise HTTPException(status_code=400, detail="Invalid product size range")

    application = (req.application or "fragment").lower()
    if application not in APPLICATIONS:
        raise HTTPException(status_code=400,
                            detail=f"Unknown application '{application}'. "
                                   f"Choose one of {APPLICATIONS}.")
    if application == "illumina" and product_max > ILLUMINA_MAX_AMPLICON_BP:
        raise HTTPException(
            status_code=400,
            detail=f"Illumina primer design is limited to amplicons \u2264 "
                   f"{ILLUMINA_MAX_AMPLICON_BP} bp; requested max is {product_max} bp.",
        )

    seq_args: Dict[str, Any] = {"SEQUENCE_TEMPLATE": template}
    if excluded:
        # For applications that intend to AMPLIFY the feature (Sanger reads
        # across the feature; Illumina sequences across it), the excluded
        # region IS the target region — primer3 must place primers OUTSIDE
        # it AND make the amplicon span it. SEQUENCE_TARGET enforces both.
        # For fragment cloning the original semantics (primers just avoid
        # the region) are preserved.
        if application in ("sanger", "illumina"):
            seq_args["SEQUENCE_TARGET"] = [[int(excluded[0]), int(excluded[1])]]
        else:
            seq_args["SEQUENCE_EXCLUDED_REGION"] = [[int(excluded[0]), int(excluded[1])]]

    # Sanger uses banded design: split the [min, max] range into 3
    # equal product-size bands and run primer3 once per band asking for
    # PRIMER_NUM_RETURN=3, giving ~9 candidates spread across the
    # length range instead of all clustered at primer3's optimum. Other
    # applications keep the single-pass design.
    SANGER_BANDS = 3
    SANGER_PER_BAND = 3
    sanger_size_bands: List[Tuple[int, int]] = []
    if application == "sanger":
        rng = product_max - product_min
        if rng >= SANGER_BANDS * 2:
            cuts = [product_min + int(round(k * rng / SANGER_BANDS))
                    for k in range(SANGER_BANDS + 1)]
            cuts[0] = product_min
            cuts[-1] = product_max
            sanger_size_bands = [(cuts[k], cuts[k + 1]) for k in range(SANGER_BANDS)]
        else:
            sanger_size_bands = [(product_min, product_max)]

    effective_num_return = int(max(1, req.num_return))

    global_args: Dict[str, Any] = {
        "PRIMER_TASK": "generic",
        "PRIMER_NUM_RETURN": effective_num_return,
        "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
    }

    # Only set TM args if user provided; otherwise let primer3 defaults apply.
    if req.primer_min_tm is not None:
        global_args["PRIMER_MIN_TM"] = float(req.primer_min_tm)
    if req.primer_opt_tm is not None:
        global_args["PRIMER_OPT_TM"] = float(req.primer_opt_tm)
    if req.primer_max_tm is not None:
        global_args["PRIMER_MAX_TM"] = float(req.primer_max_tm)

    # Optional size constraints
    if req.primer_min_size is not None:
        global_args["PRIMER_MIN_SIZE"] = int(req.primer_min_size)
    if req.primer_opt_size is not None:
        global_args["PRIMER_OPT_SIZE"] = int(req.primer_opt_size)
    if req.primer_max_size is not None:
        global_args["PRIMER_MAX_SIZE"] = int(req.primer_max_size)

    # Run primer3 — once per Sanger band, otherwise a single pass.
    designs: List[Dict[str, Any]] = []  # {result, num, band_min, band_max, band_index}
    try:
        if application == "sanger" and sanger_size_bands:
            for bi, (bmin, bmax) in enumerate(sanger_size_bands):
                ga = dict(global_args)
                ga["PRIMER_PRODUCT_SIZE_RANGE"] = [[bmin, bmax]]
                ga["PRIMER_NUM_RETURN"] = SANGER_PER_BAND
                try:
                    r = primer3.bindings.designPrimers(seq_args, ga)
                except Exception:
                    continue
                designs.append({"result": r, "num": SANGER_PER_BAND,
                                "band_min": bmin, "band_max": bmax,
                                "band_index": bi})
        else:
            r = primer3.bindings.designPrimers(seq_args, global_args)
            designs.append({"result": r, "num": effective_num_return,
                            "band_min": product_min, "band_max": product_max,
                            "band_index": 0})
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"primer3 failed: {e}")

    # Build candidate list across ALL bands. When the excluded region is
    # present, prefer candidates that flank it but keep the others as a
    # fallback pool.
    candidates: List[Dict[str, Any]] = []
    fallback_pool: List[Dict[str, Any]] = []
    for d in designs:
        result = d["result"]
        for i in range(d["num"]):
            pair_pen = result.get(f"PRIMER_PAIR_{i}_PENALTY")
            if pair_pen is None:
                continue
            try:
                pen = float(pair_pen)
            except Exception:
                continue
            l_pos = _primer_left_pos(result, i)
            r_pos = _primer_right_pos(result, i)
            if not l_pos or not r_pos:
                continue
            l_seq = result.get(f"PRIMER_LEFT_{i}_SEQUENCE")
            r_seq = result.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
            if not l_seq or not r_seq:
                continue
            rec = {
                "index": i,
                "primer3_index": i,
                "band_index": d["band_index"],
                "band_min": d["band_min"],
                "band_max": d["band_max"],
                "result": result,
                "penalty": pen,
                "left_pos": l_pos,
                "right_pos": r_pos,
                "left_seq": l_seq,
                "right_seq": r_seq,
            }
            if excluded and not _pair_flanks_excluded(l_pos, r_pos, excluded):
                fallback_pool.append(rec)
            else:
                candidates.append(rec)

    if not candidates:
        candidates = fallback_pool

    if not candidates:
        raise HTTPException(status_code=500, detail="primer3 returned no primer pairs")

    # The downstream code refers to a single `result` dict to fetch
    # things like PRIMER_PAIR_{i}_PRODUCT_SIZE / PRIMER_LEFT_{i}_TM /
    # PRIMER_RIGHT_{i}_TM. Once we've picked the winner (`chosen`)
    # below, we rebind `result` to the candidate's source primer3
    # result so those lookups work in the banded case.

    candidate_scores: List[Dict[str, Any]] = []
    selection_method = "primer3_lowest_penalty"
    selection_rationale = ""

    if application == "sanger":
        # Score every candidate's primer pair with score_sanger_primer and
        # combine with primer3's penalty. This is what makes the Sanger
        # design considerations actually drive selection rather than just
        # appear as a post-hoc badge.
        for c in candidates:
            scores = _sanger_scores_for_pair(
                template, c["left_seq"], c["right_seq"],
                c["left_pos"], c["right_pos"], excluded,
            )
            avg = (sum(s.get("overall_score", 0.0) for s in scores) / len(scores)
                   if scores else 0.0)
            # Map primer3 penalty into a 0..100 quality score (lower penalty
            # → higher quality; clamp at 10 since beyond that is poor).
            primer3_q = max(0.0, 100.0 - 10.0 * min(c["penalty"], 10.0))
            # 0.5 / 0.5 weighting: Sanger axes (read window, template
            # structure on the AMPLICON, mispriming) get equal say with
            # primer3's Tm / dimer / GC compatibility penalty.
            combined = 0.5 * avg + 0.5 * primer3_q
            c["sanger_avg"] = avg
            c["primer3_q"] = primer3_q
            c["combined"] = combined
            candidate_scores.append({
                "index": c["index"],
                "band_index": c["band_index"],
                "band_min": c["band_min"],
                "band_max": c["band_max"],
                "primer3_penalty": round(c["penalty"], 3),
                "sanger_avg": round(avg, 1),
                "combined_score": round(combined, 1),
                "left_score": round(scores[0]["overall_score"], 1) if len(scores) > 0 else None,
                "right_score": round(scores[1]["overall_score"], 1) if len(scores) > 1 else None,
            })
        # Highest combined wins; primer3 penalty is the tiebreaker.
        candidates.sort(key=lambda c: (-c["combined"], c["penalty"]))
        chosen = candidates[0]
        selection_method = "sanger_aware"
        selection_rationale = (
            f"Picked candidate #{chosen['index']} (band {chosen['band_min']}–"
            f"{chosen['band_max']} bp) from {len(candidates)} primer3 pairs "
            f"across {len(designs)} size band(s) by combined Sanger "
            f"{chosen['sanger_avg']:.0f}/100 + primer3 penalty "
            f"{chosen['penalty']:.2f}."
        )
    else:
        candidates.sort(key=lambda c: c["penalty"])
        chosen = candidates[0]
        selection_rationale = (
            f"Picked candidate #{chosen['index']} (primer3 penalty "
            f"{chosen['penalty']:.2f}) from {len(candidates)} pairs."
        )

    best_i = chosen["index"]
    best_pen = chosen["penalty"]
    # Banded design uses one primer3 result dict per band; switch to the
    # winner's dict so downstream `result.get(f"PRIMER_…_{i}")` lookups
    # resolve to the right candidate.
    result = chosen["result"]

    i = best_i
    left_seq = result.get(f"PRIMER_LEFT_{i}_SEQUENCE")
    right_seq = result.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
    if not left_seq or not right_seq:
        raise HTTPException(status_code=500, detail="primer3 missing primer sequences for selected pair")

    left_pos = _primer_left_pos(result, i)
    right_pos = _primer_right_pos(result, i)

    # Mispriming: report all exact occurrences (including the “true” site)
    left_misprime = _mispriming_sites(template, left_seq)
    right_misprime = _mispriming_sites(template, right_seq)

    # Annealing portions = exactly what primer3 designed against the template.
    left_anneal = left_seq
    right_anneal = right_seq

    # For Illumina, primers ship as adapter+anneal. The Tm primer3 returned
    # is for the annealing portion only (which is what we want — the adapter
    # is non-templated overhang and is not part of the binding event).
    if application == "illumina":
        full_left = ILLUMINA_FWD_ADAPTER + left_anneal
        full_right = ILLUMINA_REV_ADAPTER + right_anneal
    else:
        full_left = left_anneal
        full_right = right_anneal

    sid = ensure_session_id(req.session_id)

    return {
        "sessionId": sid,
        "fragments_in": template,
        "application": application,

        "pair_index": i,
        "pair_penalty": best_pen,
        "product_size": result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
        "selection_method": selection_method,
        "selection_rationale": selection_rationale,
        "candidate_scores": candidate_scores,
        "num_candidates_considered": len(candidates),

        "excluded_region": {"start": excluded[0], "length": excluded[1]} if excluded else None,

        # Full ordered primer (adapter+anneal for Illumina, anneal-only otherwise).
        "left_primer": full_left,
        "right_primer": full_right,

        # Annealing portion only — what binds the template, what gets annotated
        # on the plasmid, and what the reported Tm refers to.
        "left_annealing": left_anneal,
        "right_annealing": right_anneal,

        # Adapter overhangs (empty unless Illumina).
        "left_adapter": ILLUMINA_FWD_ADAPTER if application == "illumina" else "",
        "right_adapter": ILLUMINA_REV_ADAPTER if application == "illumina" else "",

        # Tm refers to the annealing sequence (primer3's value, unchanged).
        "left_tm": result.get(f"PRIMER_LEFT_{i}_TM"),
        "right_tm": result.get(f"PRIMER_RIGHT_{i}_TM"),

        # Binding positions on the template — annealing portion only.
        "left_pos": {"start": left_pos[0], "len": left_pos[1]} if left_pos else None,
        "right_pos": {"start_3prime": right_pos[0], "len": right_pos[1]} if right_pos else None,

        # “Mispriming” locations (exact hits) — annealing portion only.
        "left_mispriming_sites": left_misprime,
        "right_mispriming_sites": right_misprime,

        # Thermo-ish scores (computed on annealing portion).
        "left_scores": _thermo_scores(left_anneal),
        "right_scores": _thermo_scores(right_anneal),

        # Sanger-sequencing-quality scores (per primer). The target the
        # primer is meant to read across is the midpoint of the excluded
        # region when one is given, otherwise the centre of the amplicon.
        "sanger_scores": _sanger_scores_for_pair(
            template, left_anneal, right_anneal,
            left_pos, right_pos, excluded,
        ),
    }
