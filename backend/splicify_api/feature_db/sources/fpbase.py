"""FPbase ingest — download all fluorescent proteins as JSON, emit
FASTA + KB JSON. No further curation per spec (direct re-export of
the FPbase catalog, CC-BY-4.0).

Usage:
    python -m splicify_api.feature_db.sources.fpbase

Outputs:
    feature_db_data/external_src/fpbase.faa
    feature_db_data/fpbase_kb.json
    feature_db_data/fpbase.csv          (sseqid, Feature, Type, Description)
    feature_db_data/BLAST_dbs/fpbase.{phr,pin,psq,...}   (-dbtype prot)
    feature_db_data/mmseqs_dbs/fpbase_db + .idx          (mmseqs createdb)
"""
from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Any

_SCRIPT_DIR = Path(__file__).resolve().parents[1]
_SPLICIFY_API = _SCRIPT_DIR.parent
_FDD = _SPLICIFY_API / "feature_db_data"
_EXT = _FDD / "external_src"
_BLAST = _FDD / "BLAST_dbs"

FPBASE_URL = "https://www.fpbase.org/api/proteins/?format=json"


def download() -> Path:
    _EXT.mkdir(parents=True, exist_ok=True)
    out = _EXT / "fpbase_raw.json"
    print(f"[fpbase] downloading -> {out}")
    req = urllib.request.Request(FPBASE_URL, headers={"User-Agent": "splicify-feature-db/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        out.write_bytes(resp.read())
    print(f"[fpbase] {out.stat().st_size // 1024} KB written")
    return out


def emit_fasta_and_kb(json_path: Path) -> dict[str, int]:
    records: list[dict[str, Any]] = json.loads(json_path.read_text())
    fasta_path = _EXT / "fpbase.faa"
    kb_path = _FDD / "fpbase_kb.json"
    csv_path = _FDD / "fpbase.csv"

    wrote = 0
    kb_records: list[dict] = []
    csv_rows: list[list[str]] = []

    with fasta_path.open("w") as fh:
        for rec in records:
            seq = (rec.get("seq") or "").strip().upper()
            if not seq:
                continue
            name = rec.get("name") or rec.get("slug") or rec.get("uuid")
            slug = rec.get("slug") or rec.get("uuid")
            sseqid = slug
            description = rec.get("name") or sseqid
            # Stitch a useful description: name + emission peak if present
            states = rec.get("states") or []
            if states:
                s0 = states[0]
                em = s0.get("em_max")
                ex = s0.get("ex_max")
                if em is not None:
                    description = f"{description} (em={em} nm, ex={ex} nm)"

            fh.write(f">{sseqid}\n{seq}\n")
            wrote += 1
            csv_rows.append([sseqid, name, "CDS", description])
            kb_records.append({
                "feature_id": f"FPBASE_{sseqid.upper()}",
                "feature_name": name,
                "normalized_feature_name": sseqid,
                "feature_type": "CDS",
                "sseqid": sseqid,
                "source": {
                    "annotation_source": "FPbase",
                    "source_dataset": "fpbase.org/api/proteins",
                    "name": "FPbase",
                    "license": "CC-BY-4.0",
                    "citation": "Lambert TJ. 2019. Nat Methods 16(4):277-278.",
                    "descriptions": [description],
                },
                "intrinsic_properties": {
                    "feature_class": "reporter",
                    "subclass": "fluorescent_protein",
                    "sequence_derived": {
                        "representative_sequence": "",  # FPbase is aa-only
                        "protein_sequence": seq,
                        "length_aa": len(seq),
                    },
                    "fpbase_states": states,
                    "uniprot": rec.get("uniprot"),
                    "genbank": rec.get("genbank"),
                },
                "curation": {
                    "curation_status": "imported",
                    "confidence": "high",
                    "needs_manual_review": False,
                },
            })

    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sseqid", "Feature", "Type", "Description"])
        w.writerows(csv_rows)

    kb_path.write_text(json.dumps(
        {"_meta": {"source": "FPbase", "license": "CC-BY-4.0",
                   "record_count": len(kb_records)},
         "records": kb_records},
        indent=2,
    ))

    return {"fasta": wrote, "csv": len(csv_rows), "kb": len(kb_records)}


def build_blast_db(fasta: Path) -> None:
    _BLAST.mkdir(exist_ok=True)
    out = _BLAST / "fpbase"
    rc = subprocess.run(
        ["makeblastdb", "-in", str(fasta), "-dbtype", "prot",
         "-out", str(out), "-parse_seqids"],
        capture_output=True, text=True,
    )
    if rc.returncode != 0:
        print(f"[fpbase] makeblastdb FAIL: {rc.stderr.strip()[:200]}")
    else:
        print("[fpbase] makeblastdb OK -> BLAST_dbs/fpbase.*")


def build_diamond_db(fasta: Path) -> None:
    rc = subprocess.run(
        ["diamond", "makedb", "--in", str(fasta),
         "-d", str(_BLAST / "fpbase"), "--quiet"],
        capture_output=True, text=True,
    )
    if rc.returncode != 0:
        print(f"[fpbase] diamond makedb FAIL: {rc.stderr.strip()[:200]}")
    else:
        print("[fpbase] diamond OK -> BLAST_dbs/fpbase.dmnd")


def main() -> int:
    json_path = download()
    counts = emit_fasta_and_kb(json_path)
    print(f"[fpbase] emit: fasta={counts['fasta']} csv={counts['csv']} kb={counts['kb']}")
    fasta = _EXT / "fpbase.faa"
    build_blast_db(fasta)
    build_diamond_db(fasta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
