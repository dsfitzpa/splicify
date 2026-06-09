from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

router = APIRouter()

# --- Tunables (env overridable) ---
MAX_INVENTORY_FILES = int(os.environ.get("INV_GIB_MAX_INVENTORY_FILES", "10"))
SEED_K = int(os.environ.get("INV_GIB_SEED_K", "20"))  # k-mer seed size
MIN_MATCH_BP = int(os.environ.get("INV_GIB_MIN_MATCH_BP", "200"))  # min contiguous exact match
MAX_KMER_POS_PER_SEED = int(os.environ.get("INV_GIB_MAX_KMER_POS_PER_SEED", "200"))  # cap to avoid repeats
MAX_TOTAL_SEED_HITS = int(os.environ.get("INV_GIB_MAX_TOTAL_SEED_HITS", "4000"))  # per inventory seq
MIN_NEW_COVERAGE_BP = int(os.environ.get("INV_GIB_MIN_NEW_COVERAGE_BP", "50"))  # accept a hit if it adds this much
SYNTH_GAP_MIN_BP = int(os.environ.get("INV_GIB_SYNTH_GAP_MIN_BP", "20"))  # fill gaps > this with synth

RE_DNA_RUN = re.compile(r"[ACGTNacgtn]{20,}")

def _parse_bool(x: Any, default: Optional[bool] = None) -> Optional[bool]:
    if x is None:
        return default
    s = str(x).strip().lower()
    if s in ("1", "true", "t", "yes", "y", "on"):
        return True
    if s in ("0", "false", "f", "no", "n", "off"):
        return False
    if s == "":
        return default
    return default

def _sanitize_name(name: str) -> str:
    name = (name or "").strip()
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    return name[:120] if name else "seq"


def _guess_basename(upload_filename: Optional[str]) -> str:
    if not upload_filename:
        return "plasmid"
    return _sanitize_name(Path(upload_filename).stem or "plasmid")


def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def _read_seq_from_upload(raw: bytes, filename: Optional[str]) -> Tuple[str, str]:
    """
    Returns: (record_name, SEQ upper)
    Accepts GenBank (.gb/.gbk/.genbank) or FASTA (.fa/.fasta/.fna).
    """
    from io import StringIO
    from Bio import SeqIO

    if not raw:
        raise HTTPException(status_code=400, detail="Empty upload")

    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")

    fmt = "genbank"
    if filename:
        ext = Path(filename).suffix.lower()
        if ext in [".fa", ".fasta", ".fna"]:
            fmt = "fasta"

    rec = SeqIO.read(StringIO(text), fmt)
    seq = str(rec.seq).upper()
    if not seq:
        raise HTTPException(status_code=400, detail=f"No sequence found in {filename or 'upload'}")

    record_name = _sanitize_name(getattr(rec, "name", "") or getattr(rec, "id", "") or _guess_basename(filename))
    return record_name, seq


@dataclass(frozen=True)
class Hit:
    inv_name: str
    inv_len: int
    inv_orientation: str  # "+" or "-"
    inv_start: int
    inv_end: int
    tgt_start: int  # normalized to [0, L)
    tgt_end: int    # normalized to (tgt_start, tgt_start+L], may be > L to represent wrap
    match_len: int

    def wraps(self, L: int) -> bool:
        return self.tgt_end > L


def _build_kmer_index(target2: str, k: int, max_pos_per_seed: int) -> Dict[str, List[int]]:
    idx: Dict[str, List[int]] = {}
    n = len(target2)
    if n < k:
        return idx

    for i in range(0, n - k + 1):
        kmer = target2[i : i + k]
        if "N" in kmer:
            continue
        arr = idx.get(kmer)
        if arr is None:
            idx[kmer] = [i]
        else:
            if len(arr) < max_pos_per_seed:
                arr.append(i)
    return idx


def _extend_exact(
    target2: str,
    inv2: str,
    tgt_pos: int,
    inv_pos: int,
    k: int,
) -> Tuple[int, int, int, int, int]:
    """
    Extend exact match from a k-mer seed between target2 and inv2.
    """
    t0 = tgt_pos
    i0 = inv_pos
    while t0 > 0 and i0 > 0 and target2[t0 - 1] == inv2[i0 - 1]:
        t0 -= 1
        i0 -= 1

    t1 = tgt_pos + k
    i1 = inv_pos + k
    tN = len(target2)
    iN = len(inv2)
    while t1 < tN and i1 < iN and target2[t1] == inv2[i1]:
        t1 += 1
        i1 += 1

    return t0, t1, i0, i1, (t1 - t0)


def _normalize_target2_interval_to_circular_first_copy(L: int, t0: int, t1: int) -> Optional[Tuple[int, int]]:
    """
    Return (t0n, t1n) where:
      - t0n in [0, L)
      - t1n in (t0n, t0n+L]
      - t1n > L indicates wrap across origin
    """
    if t1 <= t0:
        return None

    if t0 >= L:
        t0 -= L
        t1 -= L

    if not (0 <= t0 < L):
        return None

    if t1 > t0 + L:
        t1 = t0 + L

    return t0, t1


def _best_exact_hit_for_inventory(
    target: str,
    inv_seq: str,
    inv_name: str,
    kmer_index: Dict[str, List[int]],
    k: int,
    min_match_bp: int,
    max_total_seed_hits: int,
) -> Optional[Hit]:
    """
    Best exact contiguous match between circular target and circular inventory (both strands).
    """
    L = len(target)
    target2 = target + target

    best: Optional[Hit] = None

    def scan_one(orient_seq: str, orientation: str) -> None:
        nonlocal best

        invL = len(orient_seq)
        if invL < k or len(target2) < k:
            return

        inv2 = orient_seq + orient_seq

        seed_hits_seen = 0

        # seeds only in first copy of inventory to avoid duplicate starts
        for inv_pos in range(0, invL - k + 1):
            kmer = orient_seq[inv_pos : inv_pos + k]
            if "N" in kmer:
                continue

            tgt_positions = kmer_index.get(kmer)
            if not tgt_positions:
                continue

            seed_hits_seen += len(tgt_positions)
            if seed_hits_seen > max_total_seed_hits:
                break

            for tgt_pos in tgt_positions:
                t0, t1, i0, i1, mlen = _extend_exact(target2, inv2, tgt_pos, inv_pos, k)
                if mlen < min_match_bp:
                    continue

                # dedupe inventory starts from second copy
                if i0 >= invL:
                    continue

                norm = _normalize_target2_interval_to_circular_first_copy(L, t0, t1)
                if norm is None:
                    continue
                t0n, t1n = norm

                if best is None or mlen > best.match_len or (mlen == best.match_len and t0n < best.tgt_start):
                    best = Hit(
                        inv_name=inv_name,
                        inv_len=invL,
                        inv_orientation=orientation,
                        inv_start=i0,
                        inv_end=i1,
                        tgt_start=t0n,
                        tgt_end=t1n,
                        match_len=mlen,
                    )

    scan_one(inv_seq, "+")
    scan_one(_revcomp(inv_seq), "-")

    return best


def _mark_covered(covered: bytearray, start: int, end: int) -> int:
    L = len(covered)
    start = max(0, min(L, start))
    end = max(0, min(L, end))
    if end <= start:
        return 0
    new = 0
    for i in range(start, end):
        if covered[i] == 0:
            covered[i] = 1
            new += 1
    return new


def _new_coverage_for_hit(L: int, covered: bytearray, t0: int, t1: int) -> int:
    if t1 <= L:
        return sum(1 for i in range(t0, t1) if covered[i] == 0)

    a_new = sum(1 for i in range(t0, L) if covered[i] == 0)
    b_end = min(L, t1 - L)
    b_new = sum(1 for i in range(0, b_end) if covered[i] == 0)
    return a_new + b_new


def _apply_hit_coverage(L: int, covered: bytearray, t0: int, t1: int) -> int:
    if t1 <= L:
        return _mark_covered(covered, t0, t1)

    new = _mark_covered(covered, t0, L)
    b_end = min(L, t1 - L)
    new += _mark_covered(covered, 0, b_end)
    return new


def _slice_circular(target: str, t0: int, t1: int) -> str:
    """
    Slice target in circular coordinates where:
      - t0 in [0,L)
      - t1 in (t0, t0+L]
      - t1 > L indicates wrap
    """
    L = len(target)
    if t1 <= L:
        return target[t0:t1]
    return target[t0:L] + target[0 : min(L, t1 - L)]


def _trim_hit_circular_edges_to_uncovered(L: int, covered: bytearray, t0: int, t1: int) -> Optional[Tuple[int, int]]:
    """
    Trim already-covered bases from the START and END of a circular interval [t0,t1)
    without splitting wrap-hits into two fragments.

    Works in "target2 coordinates":
      t0 in [0,L), t1 in (t0, t0+L]
    """
    if t1 <= t0:
        return None

    a = t0
    b = t1

    # Trim left edge
    while a < b and covered[a % L] == 1:
        a += 1

    # Trim right edge
    while b > a and covered[(b - 1) % L] == 1:
        b -= 1

    if b <= a:
        return None

    return a, b


def _compute_gaps(target: str, covered: bytearray, min_gap: int) -> List[Dict[str, Any]]:
    L = len(target)
    gaps: List[Tuple[int, int]] = []

    i = 0
    while i < L:
        if covered[i] == 1:
            i += 1
            continue
        j = i
        while j < L and covered[j] == 0:
            j += 1
        gaps.append((i, j))
        i = j

    if not gaps:
        return []

    # merge wrapping gap
    if gaps[0][0] == 0 and gaps[-1][1] == L and len(gaps) >= 2:
        merged = (gaps[-1][0], gaps[0][1])
        gaps = [merged] + gaps[1:-1]

    synth: List[Dict[str, Any]] = []
    sidx = 1
    for (a, b) in gaps:
        glen = b - a
        if glen <= min_gap:
            continue
        synth.append(
            {
                "name": f"Synthesis_{sidx}",
                "start": a,
                "end": b,
                "length_bp": glen,
                "sequence": target[a:b],
            }
        )
        sidx += 1
    return synth


@router.post("/inv-gib")
async def inv_gib_endpoint(
    file: UploadFile = File(...),
    inventory_files: List[UploadFile] = File(default=[]),
    session_id_form: str = Form("", alias="session_id"),
    sessionId_form: str = Form("", alias="sessionId"),
    include_ai_explanation_form: str = Form("", alias="include_ai_explanation"),
    message: str = Form(""),
):
    raw = await file.read()
    session_id = (session_id_form or sessionId_form or "").strip()
    include_ai_explanation = _parse_bool(include_ai_explanation_form, default=None)
    basename = _guess_basename(file.filename)

    inv_list: List[UploadFile] = list(inventory_files or [])
    if not inv_list:
        raise HTTPException(status_code=400, detail="inv_gib requires inventory_files (one or more plasmids)")
    if len(inv_list) > MAX_INVENTORY_FILES:
        raise HTTPException(status_code=400, detail=f"Too many inventory files: {len(inv_list)} (max {MAX_INVENTORY_FILES})")

    _, target_seq = _read_seq_from_upload(raw, file.filename)
    L = len(target_seq)
    if L < SEED_K:
        raise HTTPException(status_code=400, detail=f"Target too short for SEED_K={SEED_K}")

    target2 = target_seq + target_seq
    kmer_index = _build_kmer_index(target2, SEED_K, MAX_KMER_POS_PER_SEED)

    inv_display_names: Dict[str, str] = {}

    hits: List[Hit] = []
    inv_debug: List[Dict[str, Any]] = []

    for inv in inv_list:
        inv_raw = await inv.read()
        inv_bn = _guess_basename(inv.filename)
        inv_rec_name, inv_seq = _read_seq_from_upload(inv_raw, inv.filename)

        inv_key = inv_bn or inv_rec_name or f"inventory_{len(inv_display_names)+1}"
        inv_display_names[inv_key] = inv_bn or inv_rec_name

        best = _best_exact_hit_for_inventory(
            target=target_seq,
            inv_seq=inv_seq,
            inv_name=inv_key,
            kmer_index=kmer_index,
            k=SEED_K,
            min_match_bp=MIN_MATCH_BP,
            max_total_seed_hits=MAX_TOTAL_SEED_HITS,
        )

        if best:
            hits.append(best)
            inv_debug.append(
                {
                    "inventory": inv_display_names.get(inv_key, inv_key),
                    "best_match_len": best.match_len,
                    "orientation": best.inv_orientation,
                    "tgt_start": best.tgt_start,
                    "tgt_end": best.tgt_end,
                    "wraps": bool(best.tgt_end > L),
                }
            )
        else:
            inv_debug.append({"inventory": inv_display_names.get(inv_key, inv_key), "best_match_len": 0, "note": "no hit >= MIN_MATCH_BP"})

    # Choose hits by added coverage (selection coverage)
    hits_sorted = sorted(hits, key=lambda h: h.match_len, reverse=True)
    covered_select = bytearray(L)
    chosen: List[Hit] = []
    for h in hits_sorted:
        new_cov = _new_coverage_for_hit(L, covered_select, h.tgt_start, h.tgt_end)
        if new_cov >= MIN_NEW_COVERAGE_BP:
            _apply_hit_coverage(L, covered_select, h.tgt_start, h.tgt_end)
            chosen.append(h)

    # Emit fragments WITHOUT splitting wrap hits:
    #   - trim only at the circular edges against already-emitted coverage
    #   - keep a single fragment with tgt_end possibly > L
    fragments_out: List[Dict[str, Any]] = []
    covered_emit = bytearray(L)

    for idx, h in enumerate(chosen, start=1):
        trimmed = _trim_hit_circular_edges_to_uncovered(L, covered_emit, h.tgt_start, h.tgt_end)
        if trimmed is None:
            continue
        t0, t1 = trimmed

        # Apply coverage for the trimmed hit
        _apply_hit_coverage(L, covered_emit, t0 % L, t1 if t1 <= L else t1)

        seg_seq = _slice_circular(target_seq, t0 % L, t1)

        fragments_out.append(
            {
                "name": f"InvFrag_{idx}_{_sanitize_name(inv_display_names.get(h.inv_name, h.inv_name))}",
                "source_inventory": inv_display_names.get(h.inv_name, h.inv_name),
                "source_orientation": h.inv_orientation,
                "target_start": int(t0 % L),
                "target_end": int(t1),  # may be > L (wrap) — stays one fragment
                "wraps": bool(t1 > L),
                "length_bp": int(len(seg_seq)),
                "sequence": seg_seq,
            }
        )

    # Compute synth gaps relative to emitted coverage (not selection coverage)
    synth_out = _compute_gaps(target_seq, covered_emit, SYNTH_GAP_MIN_BP)

    covered_bp = int(sum(covered_emit))
    uncovered_bp = L - covered_bp
    reply_lines = [
        f"**inv_gib fragment search complete** for **{basename}** (inventory-only, circular-aware).",
        f"- Target length: {L} bp",
        f"- Inventory plasmids: {len(inv_list)}",
        f"- Matched fragments emitted: {len(fragments_out)}",
        f"- Covered bp: {covered_bp} ({(covered_bp / L * 100.0):.1f}%)",
        f"- Gaps > {SYNTH_GAP_MIN_BP} bp synthesized: {len(synth_out)}",
    ]

    return {
        "ok": True,
        "sessionId": session_id,
        "include_ai_explanation": include_ai_explanation,
        "reply": "\n".join(reply_lines),
        "fragments_in": fragments_out + synth_out,
        "target_sequence": target_seq,
        "anneal_overrides": [],
        "inv_gib_summary": {
            "target_len": L,
            "chosen_inventory_hits": len(chosen),
            "emitted_inventory_fragments": len(fragments_out),
            "covered_bp": covered_bp,
            "uncovered_bp": uncovered_bp,
            "synth_gap_count": len(synth_out),
            "seed_k": SEED_K,
            "min_match_bp": MIN_MATCH_BP,
        },
        "meta": {
            "record_name": basename,
            "input_filename": file.filename,
            "uploaded_mime_type": file.content_type,
            "message": message,
            "inventory_files_count": len(inv_list),
            "notes": (
                "Wrap-hit fix: trimming now happens on circular hit edges without splitting wrap matches into two linear pieces. "
                "Junction duplication fix retained by trimming already-covered bases at hit start/end."
            ),
            "inventory_debug": inv_debug,
        },
    }
