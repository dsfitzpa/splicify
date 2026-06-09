"""Output emitter helpers shared across guides_csv / guides_gb / parts_order /
protocol / workflow_trace / assembled_gb."""
from __future__ import annotations

import re
from typing import Any, Iterable, Optional


_RESIDUE_RX = re.compile(r"([A-Z]\d+[A-Z*]|[A-Z]\d+ins|[A-Z]\d+del)")


def _safe_token(s: str, max_len: int = 40) -> str:
    """Sanitise a string for use as a filename prefix: keep [A-Za-z0-9_],
    collapse runs, trim length."""
    if not s:
        return ""
    s = re.sub(r"[^A-Za-z0-9_]+", "_", s.strip())
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:max_len]


def derive_descriptor(args: dict[str, Any]) -> Optional[str]:
    """Best-effort prefix derivation from emit-args. Used only when the LLM
    didn't pass an explicit `descriptor`. Combines distinct residue codes
    found in pegRNA / guide names (e.g. 'D431S_K479E')."""
    explicit = args.get("descriptor")
    if isinstance(explicit, str) and explicit.strip():
        return _safe_token(explicit)

    residues: list[str] = []
    seen: set[str] = set()
    for key in ("pegrnas", "guides"):
        for item in (args.get(key) or []):
            if not isinstance(item, dict):
                continue
            name = item.get("name") or ""
            for hit in _RESIDUE_RX.findall(name):
                if hit not in seen:
                    seen.add(hit)
                    residues.append(hit)
            if len(residues) >= 4:
                break
    if residues:
        return _safe_token("_".join(residues))
    return None


def prefixed_filename(base: str, descriptor: Optional[str]) -> str:
    """`guides.gb` + descriptor='KEAP1_R15C' -> 'KEAP1_R15C_guides.gb'."""
    d = _safe_token(descriptor or "")
    if not d:
        return base
    return f"{d}_{base}"
