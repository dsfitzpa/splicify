#!/usr/bin/env python3
"""
Fetch full GenBank files for the 41 RefSeq accessions already curated by
feature_db/sources/refseq_plasmids.py, caching them as individual .gb files
for tokenization. Public-domain (NCBI/NLM). Polite Entrez throttle (3 rps).
"""
from __future__ import annotations
import argparse
import logging
import os
import sys
import time
from pathlib import Path

logger = logging.getLogger("ncbi_fetch")

REFSEQ_ACCESSIONS: list[tuple[str, str]] = [
    ("L08752", "pBR322"), ("L09137", "pUC19"), ("L09136", "pUC18"),
    ("X06402", "pACYC184"), ("X52327", "pBluescriptII_KSp"),
    ("Y14837", "pBluescriptII_SKp"), ("X65324", "pGEM-3Z"), ("X65307", "pGEM-7Zfp"),
    ("M77789", "pET-3a"), ("U13858", "pET-21a"), ("U13859", "pET-22b"),
    ("U13860", "pET-28a"), ("U13861", "pET-30a"), ("U13862", "pET-32a"),
    ("U57607", "pGEX-4T-1"), ("M57964", "pMAL-c2"), ("U56255", "pBAD-HisA"),
    ("AF234296", "pACYC-Duet"), ("U55761", "pcDNA3"), ("AF009656", "pcDNA3.1p"),
    ("EU546821", "pEGFP-N1"), ("U55762", "pEGFP-C1"), ("AF013230", "pECFP-N1"),
    ("U57609", "pECFP-C1"), ("GU067380", "pAAV-MCS"), ("AY037764", "pLVX"),
    ("EU258679", "pBABE-puro"), ("AF105229", "pMSCV-puro"),
    ("KC195268", "lentiCRISPR_v2"), ("MH084623", "pX330_SpCas9"),
    ("U17145", "pDONR221"), ("U17142", "pDONR207"), ("EU547637", "pDEST12.2"),
    ("L29429", "pRS316"), ("L29428", "pRS313"), ("L29430", "pRS314"),
    ("L29431", "pRS315"), ("M77790", "pYEp13"),
    ("AJ235943", "pCAMBIA1300"), ("AJ235945", "pCAMBIA2300"),
    ("AF324462", "pBeloBAC11"), ("AY247204", "pFastBac1"),
]


def fetch_all(out_dir: Path, email: str, api_key: str | None = None) -> dict[str, bool]:
    out_dir.mkdir(parents=True, exist_ok=True)
    from Bio import Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    results: dict[str, bool] = {}
    for i, (acc, name) in enumerate(REFSEQ_ACCESSIONS):
        out_path = out_dir / f"{name}__{acc}.gb"
        if out_path.exists() and out_path.stat().st_size > 500:
            logger.debug("cache hit: %s", acc)
            results[acc] = True
            continue
        if i > 0:
            time.sleep(0.35 if api_key else 0.4)  # Entrez polite rate: ~3 rps
        try:
            h = Entrez.efetch(db="nuccore", id=acc, rettype="gbwithparts", retmode="text")
            body = h.read()
            h.close()
        except Exception as exc:
            logger.warning("fetch failed for %s: %s", acc, exc)
            results[acc] = False
            continue
        if not body.startswith("LOCUS"):
            logger.warning("non-GenBank body for %s", acc)
            results[acc] = False
            continue
        out_path.write_text(body)
        logger.info("fetched %s (%s) → %s (%d bytes)", acc, name, out_path.name, len(body))
        results[acc] = True
    return results


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, default=Path("/var/data/plasmid_dbs/ncbi_refseq/"))
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", "plasmid-research@example.org"))
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args(argv)
    logging.basicConfig(level=args.log_level, format="%(asctime)s - %(levelname)s - %(message)s")
    results = fetch_all(args.out_dir, args.email, args.api_key)
    n_ok = sum(1 for v in results.values() if v)
    logger.info("DONE: %d/%d succeeded; output dir %s", n_ok, len(results), args.out_dir)
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
