from __future__ import annotations

from typing import List, Union, Optional, Literal, Dict, Any, Tuple
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field, ConfigDict, TypeAdapter

from .utils import normalize_dna, reverse_complement, safe_tm, ensure_session_id, is_valid_dna

router = APIRouter(tags=["gibson_plan"])

# Simple in-memory session cache for planning re-use (optional)
SESSIONS: Dict[str, Dict[str, Any]] = {}


class GibsonFragment(BaseModel):
    name: str
    sequence: str


class GibsonPlanRequest(BaseModel):
    model_config = ConfigDict(populate_by_name=True, extra="ignore")

    fragments: List[Union[GibsonFragment, str]]
    overlap_len: int = Field(25, ge=10, le=60)

    # Planning defaults to circular now (per your spec)
    assembly: Literal["linear", "circular"] = "circular"

    # Existing-homology detection window and threshold
    homology_window: int = Field(40, ge=10, le=200)
    homology_min: int = Field(15, ge=8, le=60)

    # ✅ IMPORTANT: these are referenced later (apply_order_and_flip)
    order: Optional[List[str]] = None
    flip: Optional[List[str]] = None

    session_id: Optional[str] = Field(default=None, alias="sessionId")


# ✅ Robust fragment coercion: accepts GibsonFragment, string, or dict {name, sequence}
FragmentAdapter = TypeAdapter(Union[GibsonFragment, str])


def normalize_fragments(frags: List[Union[GibsonFragment, str]]) -> List[GibsonFragment]:
    out: List[GibsonFragment] = []
    for i, f in enumerate(frags):
        try:
            f2 = FragmentAdapter.validate_python(f)
        except Exception:
            raise HTTPException(status_code=400, detail="Invalid fragments format")

        if isinstance(f2, GibsonFragment):
            name = f2.name or f"Fragment_{i+1}"
            seq = normalize_dna(f2.sequence)
        else:
            name = f"Fragment_{i+1}"
            seq = normalize_dna(f2)

        if not is_valid_dna(seq):
            raise HTTPException(status_code=400, detail=f"Invalid DNA in fragment {name}")

        out.append(GibsonFragment(name=name, sequence=seq))
    return out


def find_existing_overlap(
    left_seq: str, right_seq: str, window: int = 40, min_overlap: int = 15
) -> Optional[Dict[str, Any]]:
    """
    Finds longest exact suffix/prefix overlap between:
      left 3' end (last `window`) and right 5' end (first `window`)
    Returns {overlap_length, overlap_sequence} if len >= min_overlap else None.
    """
    left = normalize_dna(left_seq)
    right = normalize_dna(right_seq)
    if not left or not right:
        return None

    w = max(1, min(int(window), len(left), len(right)))
    left_tail = left[-w:]
    right_head = right[:w]

    best_k = 0
    for k in range(1, w + 1):
        if left_tail[-k:] == right_head[:k]:
            best_k = k

    if best_k >= int(min_overlap):
        return {"overlap_length": best_k, "overlap_sequence": right_head[:best_k]}
    return None


def apply_order_and_flip(
    frags: List[GibsonFragment], order: Optional[List[str]], flip: Optional[List[str]]
) -> List[Dict[str, Any]]:
    """
    Returns list of dict fragments with fields:
      {id, name, sequence, orientation}
    Orientation is "fwd" or "rev".
    """
    enriched = [
        {"id": f"frag{i+1}", "name": f.name, "sequence": f.sequence, "orientation": "fwd"}
        for i, f in enumerate(frags)
    ]

    flip_set = set([str(x) for x in (flip or [])])
    for f in enriched:
        if f["name"] in flip_set:
            f["sequence"] = reverse_complement(f["sequence"])
            f["orientation"] = "rev"

    if order:
        name_to_frag = {f["name"]: f for f in enriched}
        reordered = []
        for nm in order:
            if nm not in name_to_frag:
                raise HTTPException(status_code=400, detail=f"Order references unknown fragment name: {nm}")
            reordered.append(name_to_frag[nm])
        remaining = [f for f in enriched if f["name"] not in set(order)]
        enriched = reordered + remaining

    return enriched


def stitch_construct(
    fragments: List[Dict[str, Any]],
    assembly: Literal["linear", "circular"],
    window: int,
    min_ol: int,
) -> Tuple[str, List[Dict[str, Any]], List[Dict[str, Any]]]:
    if not fragments:
        raise HTTPException(status_code=400, detail="No fragments provided")

    annotations: List[Dict[str, Any]] = []
    junctions: List[Dict[str, Any]] = []

    seq = ""
    frag_spans: List[Dict[str, Any]] = []

    for i, f in enumerate(fragments):
        name = f["name"]
        frag_seq = f["sequence"]
        if i == 0:
            seq = frag_seq
            frag_spans.append(
                {
                    "name": name,
                    "start": 0,
                    "end": len(seq),
                    "direction": 1 if f["orientation"] == "fwd" else -1,
                }
            )
            continue

        left = fragments[i - 1]
        overlap = find_existing_overlap(left["sequence"], frag_seq, window=window, min_overlap=min_ol)

        if overlap:
            k = overlap["overlap_length"]
            seq += frag_seq[k:]

            start = len(seq) - (len(frag_seq) - k)
            end = len(seq)
            frag_spans.append(
                {
                    "name": name,
                    "start": start,
                    "end": end,
                    "direction": 1 if f["orientation"] == "fwd" else -1,
                }
            )

            left_span = frag_spans[i - 1]
            ov_start = left_span["end"] - k
            ov_end = left_span["end"]
            annotations.append(
                {
                    "name": f"Overlap {left_span['name']}→{name} (existing {k}bp)",
                    "start": ov_start,
                    "end": ov_end,
                    "direction": 1,
                }
            )
            junctions.append(
                {
                    "from": left_span["name"],
                    "to": name,
                    "overlap_sequence": overlap["overlap_sequence"],
                    "overlap_length": k,
                    "overlap_tm": safe_tm(overlap["overlap_sequence"]),
                    "source": "existing",
                }
            )
        else:
            seq += frag_seq
            start = len(seq) - len(frag_seq)
            end = len(seq)
            frag_spans.append(
                {
                    "name": name,
                    "start": start,
                    "end": end,
                    "direction": 1 if f["orientation"] == "fwd" else -1,
                }
            )
            junctions.append(
                {
                    "from": fragments[i - 1]["name"],
                    "to": name,
                    "overlap_sequence": None,
                    "overlap_length": 0,
                    "overlap_tm": None,
                    "source": "none",
                }
            )

    annotations = (
        [{"name": sp["name"], "start": sp["start"], "end": sp["end"], "direction": sp["direction"]} for sp in frag_spans]
        + annotations
    )

    if assembly == "circular" and len(fragments) >= 2:
        first = fragments[0]
        last = fragments[-1]

        overlap = find_existing_overlap(last["sequence"], first["sequence"], window=window, min_overlap=min_ol)
        if overlap:
            k = overlap["overlap_length"]
            if k < len(seq):
                # trim redundant overlap at end
                seq = seq[:-k]

                # adjust last fragment end
                last_frag = annotations[len(frag_spans) - 1]
                last_frag["end"] = max(last_frag["start"], last_frag["end"] - k)

                junctions.append(
                    {
                        "from": last["name"],
                        "to": first["name"],
                        "overlap_sequence": overlap["overlap_sequence"],
                        "overlap_length": k,
                        "overlap_tm": safe_tm(overlap["overlap_sequence"]),
                        "source": "existing",
                    }
                )

                # NOTE: For circular overlap across boundary, we annotate on 0..k
                annotations.append(
                    {
                        "name": f"Overlap {last['name']}→{first['name']} (existing {k}bp)",
                        "start": 0,
                        "end": k,
                        "direction": 1,
                    }
                )
        else:
            junctions.append(
                {
                    "from": last["name"],
                    "to": first["name"],
                    "overlap_sequence": None,
                    "overlap_length": 0,
                    "overlap_tm": None,
                    "source": "none",
                }
            )

    return seq, annotations, junctions


@router.post("/plan-gibson-assembly")
def plan_gibson_assembly(req: GibsonPlanRequest):
    frags = normalize_fragments(req.fragments)
    sid = ensure_session_id(req.session_id)

    fragments = apply_order_and_flip(frags, req.order, req.flip)

    seq, annotations, junctions = stitch_construct(
        fragments=fragments,
        assembly=req.assembly,
        window=req.homology_window,
        min_ol=req.homology_min,
    )

    construct = {
        "assembly": req.assembly,
        "sequence": seq,
        "annotations": annotations,
        "type": "gibson",
    }

    import time

    SESSIONS[sid] = {
        "saved_at": time.time(),
        "assembly": req.assembly,
        "overlap_len": req.overlap_len,
        "homology_window": req.homology_window,
        "homology_min": req.homology_min,
        "fragments": fragments,
        "construct": construct,
    }

    return {
        "sessionId": sid,
        "assembly": req.assembly,
        "overlap_len": req.overlap_len,
        "junctions": junctions,
        "construct": construct,
        "viz": {
            "type": "gibson",
            "sequence": seq,
            "annotations": annotations,
        },
    }

