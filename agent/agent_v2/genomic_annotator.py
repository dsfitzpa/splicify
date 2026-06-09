"""Annotate a genomic-kind GenBank file.

Distinct from the plasmid annotator (which is KB-driven for cloning). For
genomic gene records (e.g. KEAP1.gb), we trust the GenBank file's native
feature table because it carries authoritative gene / mRNA / CDS / exon
annotations from NCBI / RefSeq / GenBank submitters.

What this pass does:
  1. Parse the file with Biopython.
  2. Deduplicate identical features by (type, joined-coords, strand, gene).
  3. Reclassify ambiguous `misc_feature`s as CDS when the qualifier set
     strongly suggests it (UniProt / Ensembl xref, /note mentioning "exon"
     or "CDS" or a protein_id).
  4. Group features by `/gene` qualifier; collect each transcript's exon
     list via the mRNA feature's CompoundLocation.
  5. Translate every CDS — including ones spanning multiple exons (the
     joined coords are the spliced CDS sequence) and reverse-strand CDS.

Returns a digest the agent + feature_resolver can consume without ever
seeing the raw nucleotide string. The full translation string IS exposed
in the digest because it's the answer the agent needs to deliver.
"""
from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass
class GenomicFeature:
    type: str               # gene | mRNA | CDS | exon | misc_feature | ...
    gene: Optional[str]
    transcript_id: Optional[str]
    protein_id: Optional[str]
    strand: int             # +1 | -1
    start: int              # 0-indexed half-open lowest coord (across all parts)
    end: int                # 0-indexed half-open highest coord
    intervals: list[tuple[int, int]]  # ordered list of (start, end) for each join() part
    translation: Optional[str] = None  # only set for CDS / upgraded misc_feature
    qualifiers: dict[str, Any] = field(default_factory=dict)
    note: Optional[str] = None
    label: Optional[str] = None
    upgraded_from: Optional[str] = None  # set when misc_feature -> CDS


@dataclass
class GenomicAnnotation:
    organism: Optional[str]
    accession: Optional[str]
    chromosome: Optional[str]
    length_bp: int
    features: list[GenomicFeature]
    transcripts: dict[str, dict[str, Any]]  # transcript_id -> {gene, exons, protein_id?}
    genes: dict[str, dict[str, Any]]        # gene -> {start, end, strand, transcripts: [...]}


_UNIPROT_KEY = re.compile(r"uniprot", re.IGNORECASE)
_ENSEMBL_KEY = re.compile(r"ensembl|ENST\d", re.IGNORECASE)
_PROTEIN_ID_KEY = re.compile(r"^protein_id$|protein_id=", re.IGNORECASE)


def _qualifier_first(qualifiers: dict[str, list[str]], key: str) -> Optional[str]:
    v = qualifiers.get(key) if qualifiers else None
    if not v:
        return None
    return v[0] if isinstance(v, list) else v


def _flatten_intervals(location: Any) -> list[tuple[int, int]]:
    """Biopython FeatureLocation/CompoundLocation -> ordered list of (start, end)."""
    parts = getattr(location, "parts", None)
    if parts is None:
        return [(int(location.start), int(location.end))]
    return [(int(p.start), int(p.end)) for p in parts]


def _looks_like_cds(feat: Any) -> bool:
    """Should a `misc_feature` be upgraded to CDS?

    Conservative — only upgrade when the qualifier set carries strong
    CDS-specific evidence. Protein-domain misc_features (RefSeq's "Region:
    Kelch 1" / "Sensor for electrophilic agents" notes with UniProtKB
    xrefs) are NOT CDS — they're sub-regions of an already-annotated CDS.

    Strong indicators (any one is sufficient):
      - /protein_id qualifier present (NCBI CDS convention)
      - /codon_start qualifier present (only set on CDS features)
      - /translation qualifier present (pre-computed AA sequence)
      - /note explicitly contains the standalone word "exon" or the phrase
        "coding sequence" or "open reading frame" (hand-curated files
        sometimes annotate exons as misc_feature this way).

    UniProt / Ensembl xrefs alone are NOT enough — those are routinely
    attached to protein-domain region annotations.
    """
    q = feat.qualifiers or {}
    if "protein_id" in q or "codon_start" in q or "translation" in q:
        return True
    note = " ".join(q.get("note") or [])
    if not note:
        return False
    n_lower = note.lower()
    # Standalone-word check for "exon" to avoid matching e.g. "hexon".
    tokens = re.findall(r"[a-z]+", n_lower)
    if "exon" in tokens:
        return True
    return ("coding sequence" in n_lower) or ("open reading frame" in n_lower)


def _build_translation(seq: str, feat: Any) -> Optional[str]:
    """Translate a CDS feature against the parent record sequence.

    Uses the codon_start qualifier (1-indexed offset within the spliced CDS).
    Returns the protein sequence as a Python string. Stop codons render as '*'.
    Returns None on translation failure.
    """
    try:
        spliced = feat.extract(seq)  # Biopython does join() + reverse-complement.
    except Exception:
        return None
    try:
        codon_start = int((feat.qualifiers.get("codon_start") or ["1"])[0]) - 1
    except (ValueError, TypeError):
        codon_start = 0
    try:
        # Biopython's translate handles ambiguous bases / stop codons.
        from Bio.Seq import Seq
        spliced_str = str(spliced)
        if codon_start:
            spliced_str = spliced_str[codon_start:]
        # Truncate to a multiple of 3 — RefSeq CDS occasionally include the
        # stop codon; Biopython warns rather than failing, but be safe.
        spliced_str = spliced_str[: (len(spliced_str) // 3) * 3]
        return str(Seq(spliced_str).translate())
    except Exception:
        return None


def annotate_genomic_gb(gb_text: str) -> GenomicAnnotation:
    """Parse + dedup + translate a genomic-kind GenBank file."""
    from Bio import SeqIO

    handle = io.StringIO(gb_text)
    record = next(SeqIO.parse(handle, "genbank"))
    seq = str(record.seq)

    organism = None
    chromosome = None
    for f in record.features:
        if f.type == "source":
            organism = _qualifier_first(f.qualifiers, "organism")
            chromosome = _qualifier_first(f.qualifiers, "chromosome")
            break

    accession = (record.id or "").split()[0] if record.id else None
    if not accession:
        accession = record.annotations.get("accessions", [None])[0]

    # First pass — collect canonical features + upgrade candidates.
    seen: set[tuple] = set()
    features: list[GenomicFeature] = []
    for f in record.features:
        if f.type == "source":
            continue
        intervals = _flatten_intervals(f.location)
        gene = _qualifier_first(f.qualifiers, "gene")
        transcript_id = _qualifier_first(f.qualifiers, "transcript_id")
        protein_id = _qualifier_first(f.qualifiers, "protein_id")
        note = " ".join(f.qualifiers.get("note") or [])
        label = _qualifier_first(f.qualifiers, "label")
        loc_strand = getattr(f.location, "strand", None) or 1
        strand = 1 if loc_strand >= 0 else -1
        start = min(s for s, _ in intervals)
        end = max(e for _, e in intervals)

        ftype = f.type
        upgraded_from = None
        if ftype == "misc_feature" and _looks_like_cds(f):
            ftype = "CDS"
            upgraded_from = "misc_feature"

        # Dedup key — coordinate-identical features at the same gene/transcript
        # collapse to one row. This handles RefSeq files where the same CDS
        # is repeated across transcript variants.
        key = (ftype, tuple(intervals), strand, gene, transcript_id)
        if key in seen:
            continue
        seen.add(key)

        translation = None
        if ftype == "CDS":
            translation = _build_translation(seq, f)

        features.append(GenomicFeature(
            type=ftype,
            gene=gene,
            transcript_id=transcript_id,
            protein_id=protein_id,
            strand=strand,
            start=start,
            end=end,
            intervals=intervals,
            translation=translation,
            qualifiers={k: list(v) if isinstance(v, list) else v
                         for k, v in (f.qualifiers or {}).items()},
            note=note or None,
            label=label,
            upgraded_from=upgraded_from,
        ))

    # Build transcripts + genes index.
    transcripts: dict[str, dict[str, Any]] = {}
    genes: dict[str, dict[str, Any]] = {}
    for ft in features:
        if ft.transcript_id and ft.type in ("mRNA", "CDS"):
            entry = transcripts.setdefault(ft.transcript_id, {
                "gene": ft.gene, "exons": None, "cds_intervals": None,
                "protein_id": None, "strand": ft.strand,
            })
            if ft.type == "mRNA":
                entry["exons"] = ft.intervals
            elif ft.type == "CDS":
                entry["cds_intervals"] = ft.intervals
                entry["protein_id"] = ft.protein_id
                entry["translation"] = ft.translation
        if ft.type == "gene" and ft.gene:
            genes.setdefault(ft.gene, {
                "start": ft.start, "end": ft.end, "strand": ft.strand,
                "transcripts": [],
            })
        if ft.transcript_id and ft.gene:
            genes.setdefault(ft.gene, {
                "start": ft.start, "end": ft.end, "strand": ft.strand,
                "transcripts": [],
            })
            tlist = genes[ft.gene].setdefault("transcripts", [])
            if ft.transcript_id not in tlist:
                tlist.append(ft.transcript_id)

    return GenomicAnnotation(
        organism=organism,
        accession=accession,
        chromosome=chromosome,
        length_bp=len(seq),
        features=features,
        transcripts=transcripts,
        genes=genes,
    )
