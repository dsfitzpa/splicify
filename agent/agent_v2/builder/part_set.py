"""Part + PartSet: the unit of work the interpreter ships to the
builder. Every Part carries 50 bp of upstream + downstream junction
sequence from its source plasmid so the builder can decide which
construction method (Gibson, Golden Gate, restriction) is feasible
without adding new sequence."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any, Optional


JUNCTION_BP = 50


def _rc(seq: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed((seq or "").upper()))


@dataclass
class Part:
    """One annotation-derived part sliced from a source plasmid.

    `body_sequence` is the part itself; `upstream_junction` and
    `downstream_junction` are JUNCTION_BP (50) bp slices flanking
    the body in the SOURCE plasmid (NOT the eventual construct).
    The builder uses the junctions to detect Gibson overlap, Type IIs
    site coverage for Golden Gate, and shared restriction sites for
    restriction cloning.
    """
    name: str
    role: str                       # "promoter", "cds", "polya", "selection_marker", "origin", "ltr", "scaffold", etc.
    source_plasmid_id: str
    source_start: int               # 0-based, in source plasmid coords (after junction extension is applied)
    source_end: int                 # 0-based exclusive
    source_strand: int              # 1 or -1
    body_sequence: str              # the part itself (no junctions)
    upstream_junction: str          # JUNCTION_BP bp from the source plasmid, 5' of source_start
    downstream_junction: str        # JUNCTION_BP bp from the source plasmid, 3' of source_end
    length_bp: int = 0
    feature_class: Optional[str] = None
    kb_data: Optional[dict[str, Any]] = None

    def __post_init__(self):
        self.length_bp = self.length_bp or len(self.body_sequence)

    @property
    def part_id(self) -> str:
        """Stable identifier for this part — used in BuildJournal entries."""
        h = hashlib.sha1(
            (self.source_plasmid_id + str(self.source_start) +
             str(self.source_end) + str(self.source_strand)).encode()
        ).hexdigest()[:8]
        return f"part_{h}"

    def sequence_with_junctions(self, orientation: int = 1) -> str:
        """Body + flanking junctions in the requested orientation. The
        builder uses these in pair-wise compatibility checks."""
        seq = self.upstream_junction + self.body_sequence + self.downstream_junction
        if orientation == -1:
            seq = _rc(seq)
        return seq

    def body(self, orientation: int = 1) -> str:
        return self.body_sequence if orientation != -1 else _rc(self.body_sequence)

    def to_dict(self) -> dict[str, Any]:
        return {
            "part_id": self.part_id,
            "name": self.name,
            "role": self.role,
            "source_plasmid_id": self.source_plasmid_id,
            "source_start": self.source_start,
            "source_end": self.source_end,
            "source_strand": self.source_strand,
            "length_bp": self.length_bp,
            "feature_class": self.feature_class,
        }


@dataclass
class PartSet:
    """The collection of parts the interpreter resolved for a build
    request, plus coverage metadata describing which intent
    requirements are satisfied and which still have gaps."""
    parts: list[Part] = field(default_factory=list)
    # coverage: required_module/interaction → "satisfied" | "partial" | "missing"
    coverage: dict[str, str] = field(default_factory=dict)
    # parts the interpreter considered but rejected (kept around so the
    # builder can ask for them back when a swap is needed)
    candidates: list[Part] = field(default_factory=list)
    # free-text gap descriptions when something is missing
    gaps: list[str] = field(default_factory=list)

    def by_role(self, role: str) -> list[Part]:
        return [p for p in self.parts if p.role == role]

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_parts": len(self.parts),
            "parts": [p.to_dict() for p in self.parts],
            "coverage": dict(self.coverage),
            "candidates_n": len(self.candidates),
            "gaps": list(self.gaps),
        }


def extract_part(
    plasmid_id: str,
    plasmid_sequence: str,
    annotation: dict[str, Any],
    *,
    role: str,
    junction_bp: int = JUNCTION_BP,
    circular: bool = True,
) -> Optional[Part]:
    """Slice one Part out of a source plasmid's annotation, including
    junction_bp of flanking sequence on each side. Wraps around the
    sequence boundary when the plasmid is circular."""
    try:
        start = int(annotation.get("start"))
        end = int(annotation.get("end"))
    except (TypeError, ValueError):
        return None
    if start is None or end is None or end <= start:
        return None

    L = len(plasmid_sequence or "")
    if L == 0:
        return None
    strand = int(annotation.get("direction") or annotation.get("strand") or 1)
    if strand == 0:
        strand = 1

    def _slice(a: int, b: int) -> str:
        # circular-safe slice
        if a < 0:
            if circular:
                return plasmid_sequence[a % L:] + plasmid_sequence[:b]
            return plasmid_sequence[max(0, a):b]
        if b > L:
            if circular:
                return plasmid_sequence[a:] + plasmid_sequence[:b % L]
            return plasmid_sequence[a:L]
        return plasmid_sequence[a:b]

    body = _slice(start, end)
    up = _slice(start - junction_bp, start)
    dn = _slice(end, end + junction_bp)
    if not body:
        return None

    return Part(
        name=annotation.get("name") or "?",
        role=role,
        source_plasmid_id=plasmid_id,
        source_start=start,
        source_end=end,
        source_strand=strand,
        body_sequence=body.upper(),
        upstream_junction=up.upper(),
        downstream_junction=dn.upper(),
        feature_class=(annotation.get("kb_data") or {}).get("feature_class"),
        kb_data=annotation.get("kb_data"),
    )
