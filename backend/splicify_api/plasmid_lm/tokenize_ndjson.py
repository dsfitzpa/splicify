#!/usr/bin/env python3
"""
Tokenize NCBI records straight from the scraper's NDJSON output
(01_fetch_engineered_plasmids.py::raw_plasmid_records.ndjson).

Each record already has `sequence`, `accession`, `features_raw`, `definition`.
We bypass disk-materialization of GenBank files — POST the sequence directly to
/plannotate/annotate_sequence_llm and funnel the annotated result through
tokenize_plasmid(). Appends to the same JSONL corpus as tokenizer.py.

Supports resume-on-interrupt via --resume flag (reads existing JSONL, skips
already-tokenized accessions).
"""
from __future__ import annotations
import argparse
import dataclasses
import json
import logging
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from tokenizer import (
    DEFAULT_ANNOTATE_URL,
    TokenExample,
    tokenize_plasmid,
    write_vocabulary,
)
import requests
from collections import Counter

logger = logging.getLogger("tokenize_ndjson")


def infer_topology(record: dict) -> str:
    # NCBI flatfile topology not captured in record dict; default to circular
    # for engineered plasmids (the heuristic-filtered corpus is 99%+ circular).
    # Fallback to 'linear' only if clear textual signal.
    defn = (record.get("definition") or "").lower()
    if "linear" in defn and "circular" not in defn:
        return "linear"
    return "circular"


def annotate_sequence(seq: str, circular: bool, url: str, timeout: int = 180) -> dict | None:
    try:
        resp = requests.post(url,
                             json={"sequence": seq, "circular": circular, "detailed": True},
                             timeout=timeout)
        resp.raise_for_status()
        d = resp.json()
    except Exception as exc:
        logger.warning("annotate failed: %s", exc)
        return None
    flat = d.get("plannotate_annotations") or d.get("annotations") or []
    hier = d.get("hierarchical_annotations") or []
    hier_modules = [h for h in hier if h.get("layer") == "module"]
    hier_sub = [h for h in hier if h.get("layer") == "submodule"]
    if not hier_sub:
        hier_sub = d.get("cds_submodules") or d.get("cds_submodules_list") or []
    return {
        "features": flat,
        "hierarchical_annotations": hier_modules,
        "cds_submodules": hier_sub,
        "interactions": d.get("interactions") or [],
        "cloning_features": d.get("cloning_features") or {},
    }


def load_resume_set(jsonl_path: Path) -> set[str]:
    if not jsonl_path.exists():
        return set()
    done: set[str] = set()
    with open(jsonl_path) as fh:
        for line in fh:
            try:
                ex = json.loads(line)
                done.add(ex.get("plasmid_id", ""))
            except Exception:
                continue
    return done


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ndjson", type=Path, required=True)
    ap.add_argument("--output-jsonl", type=Path, required=True)
    ap.add_argument("--output-vocab", type=Path, required=True)
    ap.add_argument("--annotate-url", default=DEFAULT_ANNOTATE_URL)
    ap.add_argument("--k-rotations", type=int, default=6)
    ap.add_argument("--source-tag", default="ncbi")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--max-records", type=int, default=None)
    ap.add_argument("--progress-every", type=int, default=50)
    ap.add_argument("--min-length", type=int, default=500)
    ap.add_argument("--max-length", type=int, default=50000)
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args(argv)
    logging.basicConfig(level=args.log_level, format="%(asctime)s - %(levelname)s - %(message)s")

    done = load_resume_set(args.output_jsonl) if args.resume else set()
    logger.info("resume mode: %s; already-tokenized: %d", args.resume, len(done))

    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    fh_out = open(args.output_jsonl, "a" if args.resume else "w")

    n_done = n_skip = n_fail = n_total = 0
    vocab_counts: Counter = Counter()
    t0 = time.time()

    try:
        with open(args.ndjson) as fh_in:
            for line in fh_in:
                if args.max_records is not None and n_total >= args.max_records:
                    break
                n_total += 1
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                acc = rec.get("accession") or rec.get("accession_version") or ""
                if not acc:
                    continue
                if acc in done:
                    n_skip += 1
                    continue
                seq = (rec.get("sequence") or "").upper()
                if not seq or len(seq) < args.min_length or len(seq) > args.max_length:
                    n_skip += 1
                    continue

                circular = infer_topology(rec) == "circular"
                ann = annotate_sequence(seq, circular, args.annotate_url)
                if ann is None:
                    n_fail += 1
                    continue
                ann["_sequence"] = seq
                ann["_topology"] = "circular" if circular else "linear"
                ann["_plasmid_id"] = acc
                try:
                    examples = tokenize_plasmid(ann, args.source_tag, k_rotations=args.k_rotations)
                except Exception as exc:
                    logger.warning("tokenize failed for %s: %s", acc, exc)
                    n_fail += 1
                    continue

                for ex in examples:
                    for t in ex.tokens:
                        role = t.split(":", 1)[0] + ">" if ":" in t else t
                        vocab_counts[role] += 1
                    fh_out.write(json.dumps(dataclasses.asdict(ex), ensure_ascii=False) + "\n")
                fh_out.flush()
                n_done += 1

                if n_done and n_done % args.progress_every == 0:
                    rate = n_done / max(0.1, time.time() - t0)
                    logger.info("progress: done=%d fail=%d skip=%d total=%d rate=%.2f/s",
                                n_done, n_fail, n_skip, n_total, rate)
    finally:
        fh_out.close()

    logger.info("FINAL: done=%d fail=%d skip=%d total=%d", n_done, n_fail, n_skip, n_total)

    # Merge vocab with existing (so resume doesn't lose counts)
    if args.resume and args.output_vocab.exists():
        try:
            existing = json.load(open(args.output_vocab))
            for k, v in existing.get("role_counts", {}).items():
                vocab_counts[k] += v
        except Exception:
            pass
    write_vocabulary(args.output_vocab, vocab_counts)
    logger.info("wrote vocabulary → %s", args.output_vocab)
    return 0


if __name__ == "__main__":
    sys.exit(main())
