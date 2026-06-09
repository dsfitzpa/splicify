"""SwissProt curation — filter the full PE<=3 dump down to the
plasmid-relevant ~66,000 subset per the LLM_ANNOTATION_WORKFLOW.md
spec.

Pipeline:
    external_src/swissprot_pe123.fasta.gz   (full PE<=3 dump, ~560k)
    -> Rule A: PE==1 AND OS in plasmid-relevant organism set
    -> Rule B: gene-name whitelist overlay (catches the ~4.4k entries
       that fail Rule A but are still relevant)
    -> external_src/swissprot_curated.faa   (~66k)
    -> BLAST_dbs/swissprot.{phr,pin,psq,...}
    -> mmseqs_dbs/swissprot_db + .idx

The download URL (UniProt REST stream):
    https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+(existence:1+OR+existence:2+OR+existence:3)&format=fasta&compressed=true

Run:
    python -m splicify_api.feature_db.sources.swissprot_curate
"""
from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parents[1]
_SPLICIFY_API = _SCRIPT_DIR.parent
_FDD = _SPLICIFY_API / "feature_db_data"
_EXT = _FDD / "external_src"
_BLAST = _FDD / "BLAST_dbs"

# ---------------------------------------------------------------------------
# Rule A — PE==1 AND OS in plasmid-relevant organism set.
# Match by case-insensitive substring against the OS= field.
# Hierarchy: organism-genus catches all strains under that genus.
# ---------------------------------------------------------------------------
_RULE_A_ORG_SUBSTRINGS = (
    # E. coli + close relatives (cloning hosts)
    "Escherichia coli",
    "Salmonella",
    "Shigella",
    # Yeasts
    "Saccharomyces cerevisiae",
    "Schizosaccharomyces pombe",
    "Pichia pastoris",
    "Komagataella",            # Pichia genus reorganised
    "Kluyveromyces lactis",
    "Hansenula polymorpha",
    "Yarrowia lipolytica",
    # Mammals (expression / model orgs)
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Cricetulus griseus",      # CHO
    "Bos taurus",
    "Sus scrofa",
    # Plants
    "Arabidopsis thaliana",
    "Nicotiana benthamiana",
    "Nicotiana tabacum",
    "Solanum lycopersicum",
    "Solanum tuberosum",
    "Zea mays",
    "Oryza sativa",
    "Glycine max",
    # Invertebrate model orgs sometimes used for protein expression
    "Drosophila melanogaster",
    "Caenorhabditis elegans",
    # Insect cells (Sf9 / Sf21 baculovirus host)
    "Spodoptera frugiperda",
    # Bacterial cloning hosts / industrial / Agrobacterium
    "Bacillus subtilis",
    "Bacillus licheniformis",
    "Streptomyces",
    "Agrobacterium tumefaciens",
    "Agrobacterium rhizogenes",
    "Lactococcus lactis",
    "Pseudomonas putida",
    "Pseudomonas fluorescens",
    "Corynebacterium glutamicum",
    # Viruses commonly used as payload backbones
    "Human immunodeficiency virus",
    "Lentivirus",
    "Adeno-associated",
    "Adenovirus",
    "Vaccinia virus",
    "Sindbis virus",
    "Baculovirus",
    "Autographa californica",
    "Simian virus 40",
    "Bacteriophage",
    "Phage",
    "Enterobacteria phage",
    "Escherichia virus T",     # T4 / T7 / etc.
    "Escherichia virus lambda",
    "Cyanophage",
)

_RULE_A_ORG_RE = re.compile(
    "|".join(re.escape(s) for s in _RULE_A_ORG_SUBSTRINGS),
    re.IGNORECASE,
)

# ---------------------------------------------------------------------------
# Rule B — gene-name whitelist (matched against the GN= field,
# case-insensitive whole-word).
# Categories: selection markers, recombinases, reporters, tags,
# transactivators, editing effectors, viral payload genes.
# ---------------------------------------------------------------------------
_RULE_B_GENE_NAMES = {
    # Selection markers (bacterial)
    "bla", "ampr", "ampc", "neo", "kan", "kanr", "aph", "cat", "cm", "cmr",
    "tet", "tetr", "tetA", "tetA(C)", "hph", "hyg", "hygr",
    "puromycin", "puroR", "pac",
    "zeo", "zeor", "ble", "bler", "sh ble", "blasticidin",
    "bsr", "BlaR",
    # Selection markers (mammalian / euk)
    "hisD", "hisG", "leu2", "trp1", "ura3",
    "tk", "thymidine kinase",
    "hgprt", "dhfr",
    # Recombinases / site-specific
    "Cre", "FLP", "FLPe", "FLPo",
    "Bxb1", "Bxb1-INT", "phiC31",
    "R", "lambda integrase", "Int", "Xis",
    # Reporters
    "gfp", "GFP", "EGFP", "eGFP", "yfp", "cfp", "mvenus", "mCerulean",
    "mCherry", "mPlum", "mOrange", "mApple", "rfp", "DsRed", "tdTomato",
    "mScarlet", "mEmerald", "Citrine",
    "luciferase", "luc", "lucF", "Renilla", "ren",
    "lacZ", "beta-galactosidase", "BLG",
    "SEAP", "alkaline phosphatase",
    # Epitope / affinity tags
    "FLAG", "HA", "Myc", "V5", "GST", "MBP", "His6", "HIS", "Strep",
    "SUMO", "SBP", "CBP", "AviTag", "AVI", "TAP",
    # Transactivators / regulators
    "VP16", "VP64", "p65", "Gal4", "LexA", "tTA", "rtTA", "TetR",
    "PIT", "PIP",
    # Editing effectors
    "Cas9", "Cas12", "Cas12a", "Cas13", "cas9", "dCas9",
    "APOBEC", "APOBEC1", "APOBEC3A", "PmCDA1", "TadA",
    "DNMT3A", "TET1", "p300", "KRAB",
    "TALE", "TALEN",
    # Viral payload (HIV / lenti / AAV / baculo)
    "gag", "pol", "env", "tat", "rev", "vif", "vpr", "vpu", "nef",
    "rep", "cap", "polyprotein",
    "VP1", "VP2", "VP3",
    # Core transcription / RNAP
    "T7 RNA polymerase",
    "SP6 RNA polymerase",
    # Inducible expression
    "AraC", "araB", "lacI", "rhaR", "rhaS",
}
_RULE_B_GENES_LOWER = {g.lower() for g in _RULE_B_GENE_NAMES}


_HEADER_RE = re.compile(
    r"^>sp\|([^|]+)\|([^\s]+)\s+(.*?)"
    r"(?:\s+OS=(?P<OS>.*?))?"
    r"(?:\s+OX=(?P<OX>\d+))?"
    r"(?:\s+GN=(?P<GN>\S+))?"
    r"(?:\s+PE=(?P<PE>\d))?"
    r"(?:\s+SV=(?P<SV>\d+))?\s*$"
)


def _parse_header(header: str) -> dict[str, str]:
    """Parse a SwissProt FASTA header into a dict.

    Approach: regex over the trailing OS=/OX=/GN=/PE=/SV= block. Stable
    on the standard UniProt format used since ~2011.
    """
    m = _HEADER_RE.match(header)
    if not m:
        return {"accession": "", "entry": "", "description": header,
                "OS": "", "OX": "", "GN": "", "PE": ""}
    return {
        "accession": m.group(1),
        "entry": m.group(2),
        "description": m.group(3),
        "OS": (m.group("OS") or "").strip(),
        "OX": (m.group("OX") or "").strip(),
        "GN": (m.group("GN") or "").strip(),
        "PE": (m.group("PE") or "").strip(),
    }


def _passes_rule_a(meta: dict[str, str]) -> bool:
    if meta["PE"] != "1":
        return False
    return bool(_RULE_A_ORG_RE.search(meta["OS"]))


def _passes_rule_b(meta: dict[str, str]) -> bool:
    gn = meta["GN"].lower()
    if not gn:
        return False
    return gn in _RULE_B_GENES_LOWER


def curate(
    src_gz: Path | None = None,
    out_curated: Path | None = None,
) -> dict[str, int]:
    src_gz = src_gz or (_EXT / "swissprot_pe123.fasta.gz")
    out_curated = out_curated or (_EXT / "swissprot_curated.faa")
    out_kb = _FDD / "swissprot_curated_kb.json"
    out_csv = _FDD / "swissprot.csv"

    n_total = 0
    n_rule_a = 0
    n_rule_b = 0
    kb_records: list[dict] = []
    csv_rows: list[list[str]] = []

    out_curated.parent.mkdir(parents=True, exist_ok=True)
    print(f"[swissprot] curating {src_gz}")
    with gzip.open(src_gz, "rt") as src, out_curated.open("w") as dst:
        cur_meta: dict[str, str] = {}
        cur_seq_buf: list[str] = []

        def _flush(meta: dict[str, str], seq_parts: list[str]) -> None:
            nonlocal n_rule_a, n_rule_b
            if not meta or not seq_parts:
                return
            a = _passes_rule_a(meta)
            b = (not a) and _passes_rule_b(meta)
            if not (a or b):
                return
            if a:
                n_rule_a += 1
            else:
                n_rule_b += 1
            sseqid = meta["entry"]
            seq = "".join(seq_parts).strip()
            dst.write(f">sp|{meta['accession']}|{sseqid} {meta['description']}\n{seq}\n")
            csv_rows.append([sseqid, meta["description"][:120],
                             "CDS", meta["description"]])
            kb_records.append({
                "feature_id": f"SWISSPROT_{sseqid.upper()}",
                "feature_name": meta["description"][:120],
                "normalized_feature_name": sseqid,
                "feature_type": "CDS",
                "sseqid": sseqid,
                "source": {
                    "annotation_source": "SwissProt",
                    "source_dataset": "UniProtKB reviewed (PE=1 + Rule B)",
                    "name": "UniProtKB/Swiss-Prot",
                    "license": "CC-BY-4.0",
                    "citation": "UniProt Consortium",
                    "descriptions": [meta["description"]],
                },
                "intrinsic_properties": {
                    "feature_class": "cds_payload",
                    "subclass": "protein",
                    "OS": meta["OS"],
                    "OX": meta["OX"],
                    "GN": meta["GN"],
                    "PE": meta["PE"],
                    "rule": "A" if a else "B",
                    "sequence_derived": {
                        "representative_sequence": "",
                        "protein_sequence": seq,
                        "length_aa": len(seq),
                    },
                },
                "curation": {
                    "curation_status": "imported",
                    "confidence": "high" if a else "medium",
                    "needs_manual_review": False,
                },
            })

        for line in src:
            if line.startswith(">"):
                _flush(cur_meta, cur_seq_buf)
                n_total += 1
                cur_meta = _parse_header(line.rstrip())
                cur_seq_buf = []
            else:
                cur_seq_buf.append(line.strip())
        _flush(cur_meta, cur_seq_buf)

    out_kb.write_text(json.dumps(
        {"_meta": {"source": "UniProtKB/Swiss-Prot", "license": "CC-BY-4.0",
                   "rule_a": n_rule_a, "rule_b": n_rule_b,
                   "total_curated": n_rule_a + n_rule_b,
                   "total_input": n_total},
         "records": kb_records},
        indent=2,
    ))
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sseqid", "Feature", "Type", "Description"])
        w.writerows(csv_rows)

    print(f"[swissprot] input PE<=3: {n_total}")
    print(f"[swissprot] kept: Rule A {n_rule_a} + Rule B {n_rule_b} = {n_rule_a + n_rule_b}")
    # Emit a short KB keyed by entry_name with display fields only —
    # ~15 MB vs ~110 MB for the full curated KB. Loaded by
    # plannotate_router._load_swissprot_kb at runtime.
    short = {}
    for rec in kb_records:
        eid = rec.get("sseqid")
        if not eid:
            continue
        props = rec.get("intrinsic_properties", {}) or {}
        descs = ((rec.get("source") or {}).get("descriptions") or [])
        full_desc = descs[0] if descs else ""
        protein_name = full_desc
        for marker in (" OS=", " OX=", " GN=", " PE="):
            if marker in protein_name:
                protein_name = protein_name.split(marker, 1)[0].strip()
        short[eid] = {
            "entry_name": eid,
            "gene_name": props.get("GN", "") or "",
            "protein_name": protein_name,
            "organism": props.get("OS", "") or "",
            "taxonomy_id": props.get("OX", "") or "",
            "protein_existence": props.get("PE", "") or "",
        }
    short_path = _FDD / "swissprot_short_kb.json"
    short_path.write_text(json.dumps(short))
    print(f"[swissprot] swissprot_short_kb.json: {len(short)} entries, {short_path.stat().st_size // 1024} KB")

    return {"total_in": n_total, "rule_a": n_rule_a, "rule_b": n_rule_b,
            "kept": n_rule_a + n_rule_b}


def build_indices(fasta: Path) -> None:
    _BLAST.mkdir(exist_ok=True)

    # BLAST
    out = _BLAST / "swissprot"
    rc = subprocess.run(
        ["makeblastdb", "-in", str(fasta), "-dbtype", "prot",
         "-out", str(out), "-parse_seqids"],
        capture_output=True, text=True,
    )
    if rc.returncode != 0:
        print(f"[swissprot] makeblastdb FAIL: {rc.stderr.strip()[:300]}")
    else:
        print(f"[swissprot] makeblastdb OK -> BLAST_dbs/swissprot.*")

    # DIAMOND
    rc = subprocess.run(
        ["diamond", "makedb", "--in", str(fasta),
         "-d", str(_BLAST / "swissprot"), "--quiet"],
        capture_output=True, text=True,
    )
    if rc.returncode != 0:
        print(f"[swissprot] diamond makedb FAIL: {rc.stderr.strip()[:300]}")
    else:
        print(f"[swissprot] diamond OK -> BLAST_dbs/swissprot.dmnd")


def main() -> int:
    counts = curate()
    fasta = _EXT / "swissprot_curated.faa"
    build_indices(fasta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
