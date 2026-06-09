"""Emit feature_db_data/ artefacts from classified GenoLIB parts.

Outputs:
  feature_reference.fna  + .csv  + _kb.json   (main tier, nt)
  feature_motifs.fna     + .csv  + _kb.json   (motif tier, nt)
  feature_protein.faa    + .csv  + _kb.json   (CDS translations, aa)
  feature_knowledge_base.json                 (legacy records-shape;
                                                 plasmid_analyzer reads it)

FASTAs and CSVs are written first; BLAST DB / MMseqs DB construction is
the caller's responsibility (see build_reference_db.py).
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Optional


_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_dna(dna: str) -> str:
    """Standard-table translation. Returns "X" for ambiguous/unknown codons,
    "*" for stops. Trims to the longest codon-aligned prefix."""
    d = (dna or "").upper().replace("U", "T")
    n = (len(d) // 3) * 3
    return "".join(_CODON.get(d[i:i + 3], "X") for i in range(0, n, 3))


def _record_to_legacy_shape(part: dict, *, idx: int) -> dict:
    """Convert a GenoLIB part to the legacy `records[]` JSON shape that
    plasmid_analyzer.KnowledgeBase expects (sseqid / feature_id /
    intrinsic_properties.sequence_derived.representative_sequence)."""
    display_id = part["displayId"]
    feature_type = part.get("feature_type") or "misc_feature"
    sequence = part.get("sequence") or ""
    return {
        "feature_id": f"GENOLIB_{feature_type.upper()}_{display_id.upper()}",
        "feature_name": part.get("name") or display_id,
        "normalized_feature_name": display_id,
        "feature_type": feature_type,
        "sseqid": display_id,
        "source": {
            "annotation_source": "GenoLIB",
            "source_dataset": "labhost_All.xml",
            "descriptions": [part.get("description") or ""] if part.get("description") else [],
            "name": "GenoLIB",
            "license": part.get("source", {}).get("license", "CC-BY-4.0"),
            "citation": part.get("source", {}).get("citation", ""),
        },
        "intrinsic_properties": {
            "feature_class": _coarse_feature_class(feature_type),
            "subclass": feature_type.lower(),
            "sequence_derived": {
                "representative_sequence": sequence,
                "length_bp": len(sequence),
            },
            "hosts": part.get("hosts") or [],
        },
        "curation": {
            "curation_status": "imported",
            "confidence": "high",
            "needs_manual_review": False,
        },
        "alternative_types": [],
        "alternative_names": [],
        "alternative_classes": [],
    }


def _coarse_feature_class(feature_type: str) -> str:
    return {
        "CDS": "cds_payload",
        "promoter": "promoter",
        "terminator": "terminator",
        "rep_origin": "replication_origin",
        "enhancer": "enhancer",
        "polyA_signal": "polyA_signal",
        "polyA_site": "polyA_signal",
        "intron": "intron",
        "exon": "exon",
        "LTR": "viral_element",
        "RBS": "ribosome_binding_site",
        "misc_recomb": "recombination_site",
        "gene": "gene",
        "mobile_element": "mobile_element",
        "regulatory": "regulatory",
        "oriT": "origin_of_transfer",
        "repeat_region": "repeat",
        "tRNA": "tRNA",
        "rRNA": "rRNA",
        "ncRNA": "ncRNA",
        "mRNA": "mRNA",
        "snoRNA": "snoRNA",
        "5'UTR": "UTR",
        "3'UTR": "UTR",
        "operon": "operon",
    }.get(feature_type, "other")


def emit_tier(
    *,
    out_dir: Path,
    tier: str,                # "feature_reference" | "feature_motifs" | "feature_protein"
    parts: list[dict],
    fasta_kind: str,          # "nt" or "aa"
) -> dict[str, int]:
    """Emit FASTA + CSV + _kb.json for one tier. Returns counts."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_ext = ".fna" if fasta_kind == "nt" else ".faa"
    fasta_path = out_dir / f"{tier}{fasta_ext}"
    csv_path = out_dir / f"{tier}.csv"
    kb_path = out_dir / f"{tier}_kb.json"

    wrote_fasta = 0
    rows: list[list[str]] = []
    kb_records: list[dict] = []

    with fasta_path.open("w") as fh:
        for idx, part in enumerate(parts):
            display_id = part["displayId"]
            if fasta_kind == "aa":
                seq = translate_dna(part["sequence"])
            else:
                seq = (part.get("sequence") or "").upper()
            if not seq:
                continue
            fh.write(f">{display_id}\n{seq}\n")
            wrote_fasta += 1
            rows.append([
                display_id,
                part.get("name") or display_id,
                part.get("feature_type") or "misc_feature",
                part.get("description") or "",
            ])
            kb_records.append(_record_to_legacy_shape(part, idx=idx))

    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sseqid", "Feature", "Type", "Description"])
        w.writerows(rows)

    kb_path.write_text(json.dumps(
        {"_meta": {"tier": tier, "source": "GenoLIB", "license": "CC-BY-4.0",
                   "record_count": len(kb_records)},
         "records": kb_records},
        indent=2,
    ))

    return {"fasta": wrote_fasta, "csv": len(rows), "kb": len(kb_records)}


def emit_legacy_kb(out_dir: Path, all_parts: list[dict]) -> int:
    """Write feature_knowledge_base.json — the legacy records-shape JSON
    that plasmid_analyzer.KnowledgeBase consumes."""
    out_dir = Path(out_dir)
    kb_path = out_dir / "feature_knowledge_base.json"
    records = [_record_to_legacy_shape(p, idx=i) for i, p in enumerate(all_parts)]
    payload = {
        "_meta": {
            "source": "GenoLIB (labhost_All.xml + RefSeq gap-closure to come)",
            "license": "CC-BY-4.0",
            "record_count": len(records),
        },
        "records": records,
    }
    kb_path.write_text(json.dumps(payload, indent=2))
    return len(records)
