from __future__ import annotations

from typing import Optional, Dict, Any, List, Tuple
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field, ConfigDict
import primer3

from .utils import normalize_dna, ensure_session_id, is_valid_dna, reverse_complement

router = APIRouter(tags=["pcr"])


class BatchPrimerRequest(BaseModel):
    """
    Batch PCR primer design request.

    fragments_in: list of template sequences (strings)
    Other params match single-template endpoint.
    """
    model_config = ConfigDict(populate_by_name=True, extra="ignore")

    fragments_in: List[str]

    excluded_start: Optional[int] = None
    excluded_length: Optional[int] = None
    excluded_end: Optional[int] = None

    target_start: Optional[int] = Field(default=None, description="Legacy alias for excluded_start")
    target_length: Optional[int] = Field(default=None, description="Legacy alias for excluded_length")

    product_size_min: Optional[int] = 100
    product_size_max: Optional[int] = 300

    primer_min_tm: Optional[float] = None
    primer_opt_tm: Optional[float] = None
    primer_max_tm: Optional[float] = None

    primer_min_size: Optional[int] = None
    primer_opt_size: Optional[int] = None
    primer_max_size: Optional[int] = None

    num_return: int = 5

    session_id: Optional[str] = Field(default=None, alias="sessionId")
    include_ai_explanation: Optional[bool] = Field(default=None)

def _coalesce_excluded(req: BatchPrimerRequest) -> Optional[Tuple[int, int]]:
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
    v = result.get(f"PRIMER_LEFT_{i}")
    if isinstance(v, (list, tuple)) and len(v) >= 2:
        try:
            return (int(v[0]), int(v[1]))
        except Exception:
            return None
    return None


def _primer_right_pos(result: Dict[str, Any], i: int) -> Optional[Tuple[int, int]]:
    v = result.get(f"PRIMER_RIGHT_{i}")
    if isinstance(v, (list, tuple)) and len(v) >= 2:
        try:
            return (int(v[0]), int(v[1]))
        except Exception:
            return None
    return None


def _pair_flanks_excluded(left_pos: Tuple[int, int], right_pos: Tuple[int, int], excluded: Tuple[int, int]) -> bool:
    ex_start, ex_len = excluded
    ex_end = ex_start + ex_len

    Lstart, Llen = left_pos
    Lend = Lstart + Llen

    R3, Rlen = right_pos
    Rstart = (R3 - Rlen + 1)
    Rend = R3 + 1

    if not (Lend <= ex_start):
        return False
    if not (Rstart >= ex_end):
        return False
    if not (Lstart < Rstart):
        return False

    return True


def _mispriming_sites(template: str, primer: str) -> List[int]:
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
    return sorted(set(hits))


def _thermo_scores(primer: str) -> Dict[str, Any]:
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


def _design_one(template_in: str, req: BatchPrimerRequest, sid: str) -> Dict[str, Any]:
    template = normalize_dna(template_in)
    if not is_valid_dna(template):
        raise HTTPException(status_code=400, detail="Invalid DNA in fragments_in[]")

    excluded = _coalesce_excluded(req)
    if excluded:
        ex_start, ex_len = excluded
        if ex_start + ex_len > len(template):
            raise HTTPException(status_code=400, detail="Excluded region exceeds template length")

    product_min = int(req.product_size_min) if req.product_size_min is not None else 100
    product_max = int(req.product_size_max) if req.product_size_max is not None else 300
    if excluded:
        product_min = max(product_min, excluded[1])

    if product_min <= 0 or product_max <= 0 or product_min > product_max:
        raise HTTPException(status_code=400, detail="Invalid product size range")

    seq_args: Dict[str, Any] = {"SEQUENCE_TEMPLATE": template}
    if excluded:
        seq_args["SEQUENCE_EXCLUDED_REGION"] = [[int(excluded[0]), int(excluded[1])]]

    global_args: Dict[str, Any] = {
        "PRIMER_TASK": "generic",
        "PRIMER_NUM_RETURN": int(max(1, req.num_return)),
        "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
    }

    if req.primer_min_tm is not None:
        global_args["PRIMER_MIN_TM"] = float(req.primer_min_tm)
    if req.primer_opt_tm is not None:
        global_args["PRIMER_OPT_TM"] = float(req.primer_opt_tm)
    if req.primer_max_tm is not None:
        global_args["PRIMER_MAX_TM"] = float(req.primer_max_tm)

    if req.primer_min_size is not None:
        global_args["PRIMER_MIN_SIZE"] = int(req.primer_min_size)
    if req.primer_opt_size is not None:
        global_args["PRIMER_OPT_SIZE"] = int(req.primer_opt_size)
    if req.primer_max_size is not None:
        global_args["PRIMER_MAX_SIZE"] = int(req.primer_max_size)

    try:
        result = primer3.bindings.designPrimers(seq_args, global_args)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"primer3 failed: {e}")

    best_i: Optional[int] = None
    best_pen: float = float("inf")

    n = int(global_args.get("PRIMER_NUM_RETURN", 1))
    for i in range(n):
        pair_pen = result.get(f"PRIMER_PAIR_{i}_PENALTY")
        if pair_pen is None:
            continue
        try:
            pen = float(pair_pen)
        except Exception:
            continue

        left_pos = _primer_left_pos(result, i)
        right_pos = _primer_right_pos(result, i)
        if not left_pos or not right_pos:
            continue

        if excluded and not _pair_flanks_excluded(left_pos, right_pos, excluded):
            continue

        if pen < best_pen:
            best_pen = pen
            best_i = i

    if best_i is None:
        for i in range(n):
            pair_pen = result.get(f"PRIMER_PAIR_{i}_PENALTY")
            if pair_pen is None:
                continue
            try:
                pen = float(pair_pen)
            except Exception:
                continue
            if pen < best_pen:
                best_pen = pen
                best_i = i

    if best_i is None:
        raise HTTPException(status_code=500, detail="primer3 returned no primer pairs")

    i = best_i
    left_seq = result.get(f"PRIMER_LEFT_{i}_SEQUENCE")
    right_seq = result.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
    if not left_seq or not right_seq:
        raise HTTPException(status_code=500, detail="primer3 missing primer sequences for selected pair")

    left_pos = _primer_left_pos(result, i)
    right_pos = _primer_right_pos(result, i)

    left_misprime = _mispriming_sites(template, left_seq)
    right_misprime = _mispriming_sites(template, right_seq)

    return {
        "sessionId": sid,
        "fragments_in": template,

        "pair_index": i,
        "pair_penalty": best_pen,
        "product_size": result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),

        "excluded_region": {"start": excluded[0], "length": excluded[1]} if excluded else None,

        "left_primer": left_seq,
        "right_primer": right_seq,
        "left_tm": result.get(f"PRIMER_LEFT_{i}_TM"),
        "right_tm": result.get(f"PRIMER_RIGHT_{i}_TM"),

        "left_pos": {"start": left_pos[0], "len": left_pos[1]} if left_pos else None,
        "right_pos": {"start_3prime": right_pos[0], "len": right_pos[1]} if right_pos else None,

        "left_mispriming_sites": left_misprime,
        "right_mispriming_sites": right_misprime,

        "left_scores": _thermo_scores(left_seq),
        "right_scores": _thermo_scores(right_seq),
    }


@router.post("/batch-design-primers")
def batch_design_primers(req: BatchPrimerRequest):
    sid = ensure_session_id(req.session_id)

    if not req.fragments_in or not isinstance(req.fragments_in, list):
        raise HTTPException(status_code=400, detail="fragments_in must be a non-empty list of DNA templates")

    results: List[Dict[str, Any]] = []
    for idx, tpl in enumerate(req.fragments_in):
        one = _design_one(tpl, req, sid)
        one["template_index"] = idx
        one["template_name"] = f"Template_{idx+1}"
        results.append(one)

    return {
        "sessionId": sid,
        "include_ai_explanation": req.include_ai_explanation,
        "count": len(results),
        "results": results,
    }
