"""Build the GenoLIB-seeded feature_db_data/ tree.

Phase 1 of the LLM_ANNOTATION_WORKFLOW.md spec. Produces:
    feature_db_data/
      feature_reference.{fna,csv}  + _kb.json   (main tier, nt)
      feature_motifs.{fna,csv}     + _kb.json   (motif tier, nt)
      feature_protein.{faa,csv}    + _kb.json   (CDS translations, aa)
      feature_knowledge_base.json                (legacy records-shape)
      BLAST_dbs/
          feature_reference.{nhr,nin,nsq,...}    (makeblastdb -dbtype nucl)
          feature_motifs.{nhr,nin,nsq,...}
          feature_protein.{phr,pin,psq,...}      (-dbtype prot)
      mmseqs_dbs/
          feature_protein_db  + index            (mmseqs createdb + createindex)
    databases.yml                                 (auto-written)

Run:
    python -m splicify_api.feature_db.build_reference_db

Phases 2-5 (RefSeq gap-closure, FPbase, SwissProt curation, Rfam curation)
attach to the same tree in follow-up scripts.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

from .sources.genolib import parse_labhost
from .filters import classify_part
from .emit import emit_tier, emit_legacy_kb


_SCRIPT_DIR = Path(__file__).resolve().parent
_SPLICIFY_API = _SCRIPT_DIR.parent
_FEATURE_DB_DIR = _SPLICIFY_API / "feature_db_data"
_BLAST_DIR = _FEATURE_DB_DIR / "BLAST_dbs"
_MMSEQS_DIR = _FEATURE_DB_DIR / "mmseqs_dbs"

_GENOLIB_XML = Path("/root/genolib_supplement/sbol_extracted/SBOL_files/labhost_All.xml")


def main() -> int:
    if not _GENOLIB_XML.exists():
        print(f"FAIL: GenoLIB XML missing at {_GENOLIB_XML}", file=sys.stderr)
        return 1

    _FEATURE_DB_DIR.mkdir(parents=True, exist_ok=True)
    _BLAST_DIR.mkdir(exist_ok=True)
    _MMSEQS_DIR.mkdir(exist_ok=True)

    print(f"[1/5] Parsing GenoLIB SBOL: {_GENOLIB_XML.name}")
    raw_parts = parse_labhost(_GENOLIB_XML)
    print(f"      unique displayIds: {len(raw_parts)}")

    print("[2/5] Classifying parts (main / motif / drop)...")
    main_parts: list[dict] = []
    motif_parts: list[dict] = []
    drop_reasons: dict[str, int] = {}
    for p in raw_parts:
        cls = classify_part(
            sequence=p.get("sequence", ""),
            feature_type=p.get("feature_type", ""),
            name=p.get("name", ""),
            description=p.get("description", ""),
            license=p.get("source", {}).get("license", "CC-BY-4.0"),
        )
        if cls.tier == "main":
            main_parts.append(p)
        elif cls.tier == "motif":
            motif_parts.append(p)
        else:
            drop_reasons[cls.reason] = drop_reasons.get(cls.reason, 0) + 1
    cds_parts = [p for p in main_parts if (p.get("feature_type") or "").lower() == "cds"]
    print(f"      main: {len(main_parts)}  motif: {len(motif_parts)}  "
          f"dropped: {sum(drop_reasons.values())} ({drop_reasons})")
    print(f"      cds (will be translated for feature_protein): {len(cds_parts)}")

    print("[3/5] Emitting FASTAs / CSVs / KBs...")
    ref_counts = emit_tier(
        out_dir=_FEATURE_DB_DIR, tier="feature_reference",
        parts=main_parts, fasta_kind="nt",
    )
    motif_counts = emit_tier(
        out_dir=_FEATURE_DB_DIR, tier="feature_motifs",
        parts=motif_parts, fasta_kind="nt",
    )
    protein_counts = emit_tier(
        out_dir=_FEATURE_DB_DIR, tier="feature_protein",
        parts=cds_parts, fasta_kind="aa",
    )
    legacy_count = emit_legacy_kb(_FEATURE_DB_DIR, raw_parts)
    print(f"      feature_reference: {ref_counts['fasta']} entries")
    print(f"      feature_motifs:    {motif_counts['fasta']} entries")
    print(f"      feature_protein:   {protein_counts['fasta']} entries")
    print(f"      feature_knowledge_base.json: {legacy_count} records (all parts)")

    print("[4/5] Building BLAST DBs (makeblastdb)...")
    for tier, ext, dbtype in (
        ("feature_reference", ".fna", "nucl"),
        ("feature_motifs", ".fna", "nucl"),
        ("feature_protein", ".faa", "prot"),
    ):
        src = _FEATURE_DB_DIR / f"{tier}{ext}"
        if src.stat().st_size == 0:
            print(f"      SKIP {tier}: empty FASTA")
            continue
        out = _BLAST_DIR / tier
        rc = subprocess.run(
            ["makeblastdb", "-in", str(src), "-dbtype", dbtype,
             "-out", str(out), "-parse_seqids"],
            capture_output=True, text=True,
        )
        if rc.returncode != 0:
            print(f"      FAIL {tier}: {rc.stderr.strip()[:200]}")
        else:
            print(f"      OK   {tier} ({dbtype})")

    print("[5/5] Building DIAMOND protein DB...")
    rc = subprocess.run(
        ["diamond", "makedb",
         "--in", str(_FEATURE_DB_DIR / "feature_protein.faa"),
         "-d", str(_BLAST_DIR / "feature_protein"),
         "--quiet"],
        capture_output=True, text=True,
    )
    if rc.returncode != 0:
        print(f"      FAIL: {rc.stderr.strip()[:200]}")
    else:
        print("      OK   feature_protein.faa -> BLAST_dbs/feature_protein.dmnd")
    for prior in _MMSEQS_DIR.glob("feature_protein_db*"):
        try: prior.unlink()
        except OSError: pass
    tmp_mm.mkdir()
    cd = subprocess.run(
        ["mmseqs", "createdb",
         str(_FEATURE_DB_DIR / "feature_protein.faa"),
         str(_MMSEQS_DIR / "feature_protein_db")],
        capture_output=True, text=True,
    )
    if cd.returncode != 0:
        print(f"      FAIL createdb: {cd.stderr.strip()[:200]}")
    else:
        ci = subprocess.run(
            ["mmseqs", "createindex",
             str(_MMSEQS_DIR / "feature_protein_db"),
             str(tmp_mm), "-k", "6", "-v", "1"],
            capture_output=True, text=True,
        )
        if ci.returncode != 0:
            print(f"      WARN createindex: {ci.stderr.strip()[:200]}")
        else:
            print(f"      OK   feature_protein.faa -> mmseqs_dbs/feature_protein_db")
    shutil.rmtree(tmp_mm, ignore_errors=True)

    # Write databases.yml — spec-tuned parameters per
    # LLM_ANNOTATION_WORKFLOW.md (line 31-36). Absolute paths only.
    print("[+] Writing databases.yml...")
    yaml_path = _FEATURE_DB_DIR / "databases.yml"
    yaml_path.write_text(
        f"""# feature_db_data/databases.yml — built by feature_db/build_reference_db.py
# Tiers + parameters per LLM_ANNOTATION_WORKFLOW.md spec.
# All paths absolute; no "Default" locations.

feature_reference:
  version: GenoLIB nt main tier (CC-BY-4.0)
  method: blastn
  location: {_BLAST_DIR}
  priority: 0
  parameters:
  - -perc_identity 90
  - -word_size 11
  - -max_target_seqs 5000
  - -culling_limit 25
  - -num_threads 1
  details:
    default_type: null
    location: {_FEATURE_DB_DIR / "feature_reference.csv"}
    compressed: false

feature_motifs:
  version: GenoLIB short motifs (6-19 bp, CC-BY-4.0)
  method: blastn
  location: {_BLAST_DIR}
  priority: 3
  parameters:
  - -perc_identity 95
  - -word_size 7
  - -max_target_seqs 500
  - -culling_limit 25
  - -num_threads 1
  details:
    default_type: null
    location: {_FEATURE_DB_DIR / "feature_motifs.csv"}
    compressed: false

feature_protein:
  version: GenoLIB CDS translations (CC-BY-4.0)
  method: mmseqs
  location: {_BLAST_DIR}
  priority: 1
  parameters:
  - -e 1e-5
  - --min-seq-id 0.5
  - -s 4.0
  details:
    default_type: CDS
    location: {_FEATURE_DB_DIR / "feature_protein.csv"}
    compressed: false
""")
    print(f"      databases.yml -> {yaml_path}")

    print()
    print(f"DONE — feature_db_data/ ready at {_FEATURE_DB_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
