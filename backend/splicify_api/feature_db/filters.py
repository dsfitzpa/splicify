"""Provenance gate — classify each GenoLIB-shaped part into one of three
tiers per the LLM_ANNOTATION_WORKFLOW.md spec:

    "main"  — searched in the main BLAST/MMseqs tier (feature_reference
              for nt / feature_protein for translated CDS).
    "motif" — searched at lower priority via blastn-short with relaxed
              thresholds (preserves short operators, RBS, tag motifs).
    "drop"  — filtered out entirely.

Rules (from the spec, Provenance gate section):
    1. Licence check  — must be one of {CC-BY-4.0, CC0, public-domain}.
    2. Length floors:
         main tier:   >=20 bp non-CDS  OR  >=30 aa (>=90 bp) CDS
         motif tier:  6-19 bp
         drop:        <6 bp
    3. Primer filter  — drop primer_bind type or any name/description
       containing "primer".
    4. Low-complexity filter — drop if a single base occupies >80 % of
       the sequence.

Returns the classification + a reason string for audit/debug.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from typing import Optional

_ALLOWED_LICENSES = {"CC-BY-4.0", "CC0", "public-domain"}
_MIN_NONCDS_MAIN = 20      # bp
_MIN_CDS_MAIN = 90         # bp (= 30 aa)
_MIN_MOTIF = 6             # bp
_LOW_COMPLEXITY_FRAC = 0.80
_PRIMER_RE = re.compile(r"primer", re.IGNORECASE)


@dataclass
class PartClassification:
    tier: str        # "main" / "motif" / "drop"
    reason: str      # short audit string


def classify_part(
    *,
    sequence: str,
    feature_type: str,
    name: str = "",
    description: str = "",
    license: str = "CC-BY-4.0",
) -> PartClassification:
    """Classify a single candidate part."""
    if license not in _ALLOWED_LICENSES:
        return PartClassification("drop", f"licence:{license}")

    seq = (sequence or "").strip()
    if not seq:
        return PartClassification("drop", "no_sequence")

    length = len(seq)
    ft = (feature_type or "").lower()

    # Primer filter
    if ft == "primer_bind":
        return PartClassification("drop", "primer_type")
    if _PRIMER_RE.search(name) or _PRIMER_RE.search(description):
        return PartClassification("drop", "primer_text")

    # Low-complexity filter (most common base > 80%)
    if length >= 6:
        most_common_frac = Counter(seq.upper()).most_common(1)[0][1] / length
        if most_common_frac > _LOW_COMPLEXITY_FRAC:
            return PartClassification("drop", "low_complexity")

    # Length floors
    if length < _MIN_MOTIF:
        return PartClassification("drop", "too_short")
    if length < _MIN_NONCDS_MAIN:
        return PartClassification("motif", "len<20bp")

    # CDS gets a stricter floor (the protein search needs at least 30 aa
    # to give meaningful hits).
    if ft == "cds" and length < _MIN_CDS_MAIN:
        return PartClassification("motif", "cds<90bp")

    return PartClassification("main", "ok")
