"""
Sequence and parameter extraction from user messages.
"""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

# Common English words to skip when capturing fragment names
_SKIP_WORDS = frozenset({
    "the", "for", "and", "from", "with", "this", "that", "these",
    "have", "has", "are", "its", "use", "can", "our", "all",
    "design", "primer", "primers", "fragment", "sequence", "template",
    "assembly", "gibson", "pcr",
})


def extract_sequences(message: str) -> Dict[str, Any]:
    """
    Extract DNA sequences and named fragments from message text.

    Handles multi-line sequences where a fragment's DNA is split across
    lines or contains embedded whitespace (common when pasting from docs).

    Returns:
        sequences: ordered list of unique uppercase sequences
        fragments: dict mapping name -> sequence (preserves order)
        count: number of sequences found
    """
    sequences: List[str] = []
    fragments: Dict[str, str] = {}

    # ── Pattern 1: Numbered labels  (Frag1/Fragment 2/Template_3) ─────────────
    # For these we extract ALL DNA characters between consecutive labels so that
    # multi-line sequences (where each line is a separate DNA chunk) are joined
    # into one complete sequence per fragment.
    numbered_re = re.compile(
        r'(?:frag(?:ment)?[_\s]*(\d+)|template[_\s]*(\d+))\s*[:\s=]+',
        re.IGNORECASE,
    )
    numbered_matches = list(numbered_re.finditer(message))

    if numbered_matches:
        for idx, m in enumerate(numbered_matches):
            num_frag = m.group(1)
            num_tmpl = m.group(2)
            name = f"Fragment_{num_frag}" if num_frag else f"Template_{num_tmpl}"

            # Text span: from after this label to the start of the next label
            end_pos = numbered_matches[idx + 1].start() if idx + 1 < len(numbered_matches) else len(message)
            segment = message[m.end():end_pos]

            # Strip everything that isn't an ATGC base to join multi-line DNA
            seq = re.sub(r'[^ATGCatgc]', '', segment).upper()

            if len(seq) >= 20 and seq not in sequences:
                sequences.append(seq)
                fragments[name] = seq

    else:
        # ── Pattern 2: Generic word labels  "MyFrag: ATGC..." ─────────────────
        # Only the first continuous DNA block is captured (no multi-line join)
        # because generic labels can't reliably delimit fragment boundaries.
        for m in re.finditer(
            r'([A-Za-z]\w{0,19})[\s:=]+([ATGCatgc]{20,})',
            message,
            re.IGNORECASE,
        ):
            label = m.group(1)
            seq = m.group(2).upper()
            if label.lower() not in _SKIP_WORDS and seq not in sequences:
                sequences.append(seq)
                fragments[label] = seq

        # ── Pattern 3: Standalone sequences (≥50 bp, not already captured) ───
        for seq in re.findall(r'(?<![A-Za-z])([ATGCatgc]{50,})(?![A-Za-z])', message):
            seq_upper = seq.upper()
            if seq_upper not in sequences:
                name = f"Fragment_{len(fragments) + 1}"
                sequences.append(seq_upper)
                fragments[name] = seq_upper

    return {
        "sequences": sequences,
        "fragments": fragments,
        "count": len(sequences),
    }


def extract_parameters(message: str) -> Dict[str, Any]:
    """
    Extract design parameters from message text.

    Recognised parameters:
        target_tm        – melting temperature (°C)
        overlap_length   – Gibson overlap length (bp)
        excluded_start   – excluded region start position
        excluded_length  – excluded region length
    """
    params: Dict[str, Any] = {}
    msg = message.lower()

    # Tm: "Tm 60", "Tm of 65", "Tm=58.5"
    tm_match = re.search(r'\btm\s*(?:of\s*|=\s*)?(\d+(?:\.\d+)?)', msg)
    if tm_match:
        params["target_tm"] = float(tm_match.group(1))

    # Overlap: "overlap 30", "overlap of 25 bp"
    overlap_match = re.search(r'\boverlap\s*(?:of\s*)?(\d+)', msg)
    if overlap_match:
        params["overlap_length"] = int(overlap_match.group(1))

    # Excluded region: "exclude 100 to 200", "excluding positions 50-150"
    excl_match = re.search(
        r'\bexclu(?:de|ding)\b.*?(\d+)\s*(?:to|through|-)\s*(\d+)', msg
    )
    if excl_match:
        start = int(excl_match.group(1))
        end = int(excl_match.group(2))
        if end > start:
            params["excluded_start"] = start
            params["excluded_length"] = end - start

    return params


def redact_sequences(message: str) -> str:
    """
    Replace DNA sequences in the message with short placeholders.
    Used so the LLM receives a clean, compact message for intent detection.
    """
    counter = [0]

    def replacer(m: re.Match) -> str:
        counter[0] += 1
        seq = m.group(0)
        return f"[DNA_SEQ_{counter[0]}:{len(seq)}bp]"

    # Redact labeled sequences (e.g. "Frag1: ATGC...") keeping the label
    labeled_redacted = re.sub(
        r'(?<=[:\s=])([ATGCatgc]{20,})',
        replacer,
        message,
    )
    # Redact any remaining standalone sequences
    return re.sub(r'(?<![A-Za-z])([ATGCatgc]{50,})(?![A-Za-z])', replacer, labeled_redacted)


def build_fragment_objects(seq_data: Dict[str, Any]) -> List[Any]:
    """
    Convert extracted sequence data into GibsonFragment objects when names are available,
    otherwise return plain strings.  Import happens lazily to avoid circular deps.
    """
    from .gibson_primers import GibsonFragment

    if seq_data["fragments"]:
        return [
            GibsonFragment(name=name, sequence=seq)
            for name, seq in seq_data["fragments"].items()
        ]
    return list(seq_data["sequences"])
