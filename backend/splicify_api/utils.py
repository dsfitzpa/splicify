from __future__ import annotations

from typing import Optional
import re
import uuid
import primer3

_DNA = set("ACGTN")


def normalize_dna(seq: str) -> str:
    s = (seq or "").replace("\r", "").replace("\t", "").strip()
    # Remove FASTA headers if present
    if s.startswith(">"):
        s = "\n".join([ln.strip() for ln in s.splitlines() if ln and not ln.startswith(">")])
    s = s.replace(" ", "").replace("\n", "").upper()
    return s


def reverse_complement(seq: str) -> str:
    s = normalize_dna(seq)
    comp = str.maketrans("ACGTN", "TGCAN")
    return s.translate(comp)[::-1]


def safe_tm(seq: str) -> Optional[float]:
    try:
        s = normalize_dna(seq)
        if not s:
            return None
        if any(c not in _DNA for c in s):
            return None
        return float(primer3.bindings.calcTm(s))
    except Exception:
        return None


def ensure_session_id(maybe: Optional[str]) -> str:
    """
    Accepts:
      - None -> generates new UUID
      - UUID-like string (possibly with whitespace/newlines) -> extracts UUID
      - any other string -> returns stripped string (still stable)
    """
    if not maybe:
        return str(uuid.uuid4())
    raw = str(maybe).strip()
    m = re.search(
        r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
        raw,
        re.I,
    )
    return m.group(0) if m else raw


def is_valid_dna(seq: str) -> bool:
    return bool(seq) and all(c in _DNA for c in seq)
