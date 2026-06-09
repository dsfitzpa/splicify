#!/usr/bin/env python3
"""
Ingest LLM-generated plasmid descriptions and part annotations into a
training-pair JSONL (mirrors the `design_training_pairs` postgres schema
row-for-row, so a single `\\copy` loads it into postgres later when the
schema is applied).

Inputs
------
/var/data/plasmid_lm_corpus/token_corpus.jsonl       # rotation-0 per plasmid gives target_strings
/var/data/plasmid_lm_corpus/plasmid_descriptions.jsonl
/var/data/plasmid_lm_corpus/part_annotations.jsonl
/var/data/plasmid_lm_corpus/vb_shorthand_pairs.jsonl # existing 26 VB pairs

Output
------
/var/data/plasmid_lm_corpus/design_training_pairs.jsonl
  Each row:
    natural_language  (str)
    target_strings    (list[str]   — the plasmid's tokens at rotation 0)
    target_tokens     (null        — populated later when plasmid_tokens vocab is built)
    source_type       (str         — 'llm_user_request_short' / 'llm_user_request_functional' /
                                     'llm_lab_slack_question' / 'llm_methods_spec' /
                                     'llm_part_short' / 'llm_part_long' /
                                     'vb_shorthand')
    source_plasmid_id (str)
    confidence        (float       — 1.0 for first-attempt passes, 0.8 for retries, 0.5 for part pairs)
    is_validated      (bool)
    split             (str         — 'train' | 'val' | 'test', deterministic on plasmid_id hash)
    metadata          (dict)

Splits:
  hash(plasmid_id) % 100:   0-4 -> val  (5 %)
                            5-9 -> test (5 %)
                           10+  -> train
"""
from __future__ import annotations
import argparse
import hashlib
import json
import logging
from collections import Counter, defaultdict
from pathlib import Path

logger = logging.getLogger("ingest")


def split_for(plasmid_id: str) -> str:
    h = int(hashlib.md5(plasmid_id.encode()).hexdigest(), 16) % 100
    if h < 5:
        return "val"
    if h < 10:
        return "test"
    return "train"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/token_corpus.jsonl"))
    ap.add_argument("--plasmid-desc", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/plasmid_descriptions.jsonl"))
    ap.add_argument("--part-annot", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/part_annotations.jsonl"))
    ap.add_argument("--vb-shorthand", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/vb_shorthand_pairs.jsonl"))
    ap.add_argument("--output", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/design_training_pairs.jsonl"))
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s - %(levelname)s - %(message)s")

    # 1. Build plasmid_id → rotation-0 token_strings map (once)
    logger.info("building plasmid_id → tokens map from %s", args.corpus)
    pid_tokens: dict[str, list[str]] = {}
    with open(args.corpus) as fh:
        for line in fh:
            ex = json.loads(line)
            if ex.get("rotation_idx", 0) != 0:
                continue
            if ex["plasmid_id"] in pid_tokens:
                continue
            pid_tokens[ex["plasmid_id"]] = ex["tokens"]
    logger.info("  %d unique rotation-0 plasmids", len(pid_tokens))

    out_fh = open(args.output, "w")
    n_out = 0
    by_source_type = Counter()
    by_split = Counter()

    # 2. Plasmid-level descriptions → 4 training pairs per successful plasmid
    logger.info("ingesting plasmid descriptions from %s", args.plasmid_desc)
    n_plasmid_in = n_plasmid_written = n_skipped_no_tokens = n_skipped_failed = 0
    with open(args.plasmid_desc) as fh:
        for line in fh:
            row = json.loads(line)
            n_plasmid_in += 1
            pid = row["plasmid_id"]
            if row.get("hallucination_check", {}).get("failed"):
                n_skipped_failed += 1
                continue
            tokens = pid_tokens.get(pid)
            if not tokens:
                n_skipped_no_tokens += 1
                continue
            split = split_for(pid)
            attempts = row.get("hallucination_check", {}).get("attempts", 1)
            confidence = {1: 1.0, 2: 0.85, 3: 0.7}.get(attempts, 0.6)
            for d in row.get("descriptions", []):
                text = (d.get("text") or "").strip()
                if not text:
                    continue
                source_type = f"llm_{d['style']}"
                pair = {
                    "natural_language": text,
                    "target_strings": tokens,
                    "target_tokens": None,
                    "source_type": source_type,
                    "source_plasmid_id": pid,
                    "confidence": confidence,
                    "is_validated": False,
                    "split": split,
                    "metadata": {
                        "host": row.get("host"),
                        "source_corpus": row.get("source_corpus"),
                        "tags": row.get("tags", []),
                        "llm_attempts": attempts,
                    },
                }
                out_fh.write(json.dumps(pair, ensure_ascii=False) + "\n")
                n_out += 1
                by_source_type[source_type] += 1
                by_split[split] += 1
                n_plasmid_written += 1
    logger.info("  plasmid desc: in=%d skipped_failed=%d skipped_no_tokens=%d pairs_written=%d",
                n_plasmid_in, n_skipped_failed, n_skipped_no_tokens, n_plasmid_written)

    # 3. Part annotations → 2 pairs each (short+long), target = the part's canonical_id as a
    #    one-element token list (this teaches description → single part fill-in behaviour)
    logger.info("ingesting part annotations from %s", args.part_annot)
    n_part_in = n_part_written = 0
    with open(args.part_annot) as fh:
        for line in fh:
            row = json.loads(line)
            n_part_in += 1
            if row.get("failed"):
                continue
            cid = row.get("canonical_id")
            if not cid:
                continue
            # Build a synthetic single-token "target" for part pairs
            target_strings = [f"<{cid.replace(':', ':')}>"]
            split = split_for(cid)
            for field, source_type in (("short", "llm_part_short"), ("long", "llm_part_long")):
                text = (row.get(field) or "").strip()
                if not text:
                    continue
                pair = {
                    "natural_language": text,
                    "target_strings": target_strings,
                    "target_tokens": None,
                    "source_type": source_type,
                    "source_plasmid_id": None,
                    "confidence": 0.6,   # parts are looser / single-token targets
                    "is_validated": False,
                    "split": split,
                    "metadata": {
                        "canonical_id": cid,
                        "role": row.get("role"),
                        "occurrence_count": row.get("occurrence_count"),
                        "typical_hosts": row.get("typical_hosts", []),
                        "tags": row.get("tags", []),
                    },
                }
                out_fh.write(json.dumps(pair, ensure_ascii=False) + "\n")
                n_out += 1
                n_part_written += 1
                by_source_type[source_type] += 1
                by_split[split] += 1
    logger.info("  part annot: in=%d pairs_written=%d", n_part_in, n_part_written)

    # 4. VectorBuilder shorthand → 1 pair each
    logger.info("ingesting VB shorthand from %s", args.vb_shorthand)
    n_vb = 0
    if args.vb_shorthand.exists():
        with open(args.vb_shorthand) as fh:
            for line in fh:
                row = json.loads(line)
                nl = row.get("natural_language", "")
                target = row.get("target_tokens", [])
                pid = row.get("source_plasmid_id", "") or "vb_shorthand"
                split = split_for(pid)
                pair = {
                    "natural_language": nl,
                    "target_strings": target,
                    "target_tokens": None,
                    "source_type": "vb_shorthand",
                    "source_plasmid_id": pid,
                    "confidence": 1.0,
                    "is_validated": True,
                    "split": split,
                    "metadata": {},
                }
                out_fh.write(json.dumps(pair, ensure_ascii=False) + "\n")
                n_out += 1
                n_vb += 1
                by_source_type["vb_shorthand"] += 1
                by_split[split] += 1
        logger.info("  vb_shorthand: pairs_written=%d", n_vb)

    out_fh.close()

    print()
    print(f"=== INGEST SUMMARY ===")
    print(f"total pairs written: {n_out}")
    print(f"output: {args.output}")
    print("by source_type:")
    for s, c in by_source_type.most_common():
        print(f"  {s:30s} {c}")
    print("by split:")
    for s, c in by_split.most_common():
        print(f"  {s:6s} {c}")

    # 5. Emit a companion README fragment
    summary_path = args.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps({
        "total_pairs": n_out,
        "by_source_type": dict(by_source_type),
        "by_split": dict(by_split),
        "n_skipped_failed_grounding": n_skipped_failed,
        "n_skipped_no_tokens": n_skipped_no_tokens,
        "output_file": str(args.output),
    }, indent=2))
    print(f"\nsummary: {summary_path}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
