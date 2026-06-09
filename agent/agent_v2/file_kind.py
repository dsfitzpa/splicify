"""Classify an uploaded GenBank file as plasmid vs genomic.

The same GenBank file format is used for engineered plasmids (synthetic
constructs, circular, KB-driven annotation) and for genomic gene records
(natural organism, linear, NCBI/RefSeq-style multi-exon CDS). The agent's
behaviour must diverge between the two:

- **plasmid**: the existing KB-driven annotator runs over the sequence to
  classify modules / cloning features / interactions. Topology = circular.
- **genomic**: trust the GenBank file's native feature table. Just dedup,
  reclassify ambiguous misc_features as CDS when the /db_xref or /note
  qualifies, and walk joined exon coordinates for codon -> residue
  resolution. Topology = linear.

This module is a pure classifier — no IO, no annotation. Routers call it
once at upload time and stash the result; subagents read the cached kind
to decide which downstream pipeline to follow.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class FileKind:
    kind: str  # "plasmid" | "genomic" | "unknown"
    topology: str  # "circular" | "linear" | "unknown"
    organism: Optional[str] = None
    division: Optional[str] = None  # GenBank division: SYN | CON | PRI | ROD | ...
    accession: Optional[str] = None
    n_genes: int = 0
    n_mrna: int = 0
    n_cds: int = 0
    n_exon: int = 0
    n_misc_feature: int = 0
    has_join_coords: bool = False
    has_refseq_keyword: bool = False
    confidence: float = 0.0  # 0..1; higher = more confident in `kind`
    signals: list[str] = field(default_factory=list)  # human-readable rationale


# -- header parsers ---------------------------------------------------------
_LOCUS_RE = re.compile(
    r"^LOCUS\s+\S+\s+\d+\s+bp\s+\S+\s+(circular|linear)\s+(\S+)\s+",
    re.MULTILINE | re.IGNORECASE,
)
_ACCESSION_RE = re.compile(r"^ACCESSION\s+(\S+)", re.MULTILINE)
_ORGANISM_RE = re.compile(r"^\s+ORGANISM\s+(.+?)$", re.MULTILINE)
_SOURCE_RE = re.compile(r"^SOURCE\s+(.+?)$", re.MULTILINE)
_KEYWORDS_RE = re.compile(r"^KEYWORDS\s+(.+?)$", re.MULTILINE)
_FEATURE_TYPE_RE = re.compile(
    r"^     (gene|mRNA|CDS|exon|misc_feature|tRNA|rRNA|ncRNA|source)\s",
    re.MULTILINE,
)
_JOIN_RE = re.compile(r"\bjoin\s*\(", re.IGNORECASE)


# Accession prefixes that strongly imply genomic / RefSeq data.
_GENOMIC_ACCESSION_PREFIXES = (
    "NC_", "NG_", "NW_", "NT_", "NM_", "NR_", "NP_",
    "AC_", "AP_", "WP_", "YP_", "XM_", "XR_",
)


def _organism_signal(organism: str) -> tuple[str, str]:
    """Return (vote, signal_msg). vote in {plasmid, genomic, unknown}."""
    if not organism:
        return ("unknown", "no ORGANISM line")
    o = organism.lower()
    synthetic_markers = ("synthetic", "construct", "vector", "unidentified")
    if any(m in o for m in synthetic_markers):
        return ("plasmid", f"ORGANISM='{organism}' (synthetic)")
    # Anything with a real binomial Latin name (two words, both capitalised
    # in the original) is treated as genomic by default. We do a light
    # heuristic: at least two words, first word starts with a capital.
    parts = organism.split()
    if len(parts) >= 2 and parts[0][:1].isupper():
        return ("genomic", f"ORGANISM='{organism}' (natural organism)")
    return ("unknown", f"ORGANISM='{organism}' (unclassified)")


def _division_signal(division: str) -> tuple[str, str]:
    if not division:
        return ("unknown", "no GenBank division")
    d = division.upper()
    if d == "SYN":
        return ("plasmid", "GenBank division=SYN")
    if d in {"CON", "PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "PHG", "ENV", "UNA"}:
        return ("genomic", f"GenBank division={d}")
    return ("unknown", f"GenBank division={d}")


def _accession_signal(accession: str) -> tuple[str, str]:
    if not accession or accession == ".":
        return ("unknown", "no ACCESSION")
    # Strip the REGION suffix sometimes seen (e.g. "NC_000019 REGION:...").
    head = accession.split()[0]
    for pre in _GENOMIC_ACCESSION_PREFIXES:
        if head.startswith(pre):
            return ("genomic", f"ACCESSION='{head}' (RefSeq-style)")
    return ("unknown", f"ACCESSION='{head}'")


def classify_genbank(gb_text: str) -> FileKind:
    """Run the heuristic classifier. Returns FileKind with rationale."""
    if not gb_text:
        return FileKind(kind="unknown", topology="unknown", confidence=0.0,
                         signals=["empty input"])

    # Header parses
    locus = _LOCUS_RE.search(gb_text)
    topology = locus.group(1).lower() if locus else "unknown"
    division = locus.group(2).upper() if locus else None
    accession_match = _ACCESSION_RE.search(gb_text)
    accession = accession_match.group(1).strip() if accession_match else None
    org_match = _ORGANISM_RE.search(gb_text)
    organism = org_match.group(1).strip() if org_match else None
    keywords_match = _KEYWORDS_RE.search(gb_text)
    keywords = keywords_match.group(1).strip() if keywords_match else ""
    has_refseq_keyword = "refseq" in keywords.lower()

    # Feature-table counts (only counts at the start of feature lines so
    # /note="...gene..." etc. don't false-positive).
    feature_counts: dict[str, int] = {}
    for ft in _FEATURE_TYPE_RE.findall(gb_text):
        feature_counts[ft] = feature_counts.get(ft, 0) + 1

    # join(...) coords are an extremely strong genomic signal — multi-exon
    # CDS or mRNA records use them; synthetic plasmid annotators almost
    # never emit join() because every feature is a contiguous interval.
    has_join = bool(_JOIN_RE.search(gb_text))

    # Tally votes
    votes: dict[str, float] = {"plasmid": 0.0, "genomic": 0.0}
    signals: list[str] = []

    org_vote, org_msg = _organism_signal(organism or "")
    if org_vote != "unknown":
        votes[org_vote] += 3.0
    signals.append(org_msg)

    div_vote, div_msg = _division_signal(division or "")
    if div_vote != "unknown":
        votes[div_vote] += 2.0
    signals.append(div_msg)

    acc_vote, acc_msg = _accession_signal(accession or "")
    if acc_vote != "unknown":
        votes[acc_vote] += 2.0
    signals.append(acc_msg)

    if topology == "circular":
        votes["plasmid"] += 1.5
        signals.append("topology=circular")
    elif topology == "linear":
        votes["genomic"] += 1.0
        signals.append("topology=linear")

    if has_refseq_keyword:
        votes["genomic"] += 2.0
        signals.append("KEYWORDS includes RefSeq")

    if has_join:
        votes["genomic"] += 2.0
        signals.append("feature coords use join() — multi-exon")

    # Feature-table shape: a record with multiple mRNA + CDS + gene features
    # is genomic; a record dominated by misc_feature (typical of SnapGene
    # plasmid exports) is plasmid.
    n_cds = feature_counts.get("CDS", 0)
    n_mrna = feature_counts.get("mRNA", 0)
    n_gene = feature_counts.get("gene", 0)
    n_misc = feature_counts.get("misc_feature", 0)
    n_exon = feature_counts.get("exon", 0)
    if n_mrna >= 1 and n_gene >= 1 and n_cds >= 1:
        votes["genomic"] += 1.5
        signals.append(f"feature table has gene+mRNA+CDS (n={n_gene},{n_mrna},{n_cds})")
    if n_misc >= 10 and n_misc > 3 * max(n_cds, 1):
        votes["plasmid"] += 1.0
        signals.append(f"feature table dominated by misc_feature (n={n_misc})")

    # Decide
    if votes["genomic"] > votes["plasmid"]:
        kind = "genomic"
        confidence = min(1.0, votes["genomic"] / (votes["genomic"] + votes["plasmid"] + 0.5))
    elif votes["plasmid"] > votes["genomic"]:
        kind = "plasmid"
        confidence = min(1.0, votes["plasmid"] / (votes["genomic"] + votes["plasmid"] + 0.5))
    else:
        kind = "unknown"
        confidence = 0.0

    return FileKind(
        kind=kind,
        topology=topology,
        organism=organism,
        division=division,
        accession=accession,
        n_genes=n_gene,
        n_mrna=n_mrna,
        n_cds=n_cds,
        n_exon=n_exon,
        n_misc_feature=n_misc,
        has_join_coords=has_join,
        has_refseq_keyword=has_refseq_keyword,
        confidence=round(confidence, 3),
        signals=signals,
    )
