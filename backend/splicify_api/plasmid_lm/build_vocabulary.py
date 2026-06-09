#!/usr/bin/env python3
"""
Build the LM tokenizer vocabulary: a single integer ID space covering both
plasmid role tokens and the English words used in LLM-generated descriptions.

Outputs:
  /var/data/plasmid_lm_corpus/lm_vocabulary.json
    {
      "version": 1,
      "specials": {"<PAD>":0, "<BOS>":1, "<EOS>":2, "<UNK>":3,
                    "<DESC_START>":4, "<DESC_END>":5, "<SEP>":6},
      "role_tokens":   {"<MOD_OPEN:...>": int, ...},
      "word_tokens":   {"lentivirus": int, ...},  # lowercased, punctuation-stripped
      "vocab_size":    int
    }

Usage:
  python build_vocabulary.py \
      --corpus /var/data/plasmid_lm_corpus/token_corpus.jsonl \
      --descriptions /var/data/plasmid_lm_corpus/design_training_pairs.jsonl \
      --output /var/data/plasmid_lm_corpus/lm_vocabulary.json \
      --min-role-count 3 \
      --word-vocab-size 8000
"""
from __future__ import annotations
import argparse
import json
import logging
import re
from collections import Counter
from pathlib import Path

logger = logging.getLogger("vocab")

SPECIALS = ["<PAD>", "<BOS>", "<EOS>", "<UNK>", "<DESC_START>", "<DESC_END>", "<SEP>"]

WORD_SPLIT = re.compile(r"[A-Za-z0-9\-]+")


def tokenize_description(text: str) -> list[str]:
    """Lowercased whitespace + hyphen preserving split; a tiny English tokenizer."""
    return [w.lower() for w in WORD_SPLIT.findall(text)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/token_corpus.jsonl"))
    ap.add_argument("--descriptions", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/design_training_pairs.jsonl"))
    ap.add_argument("--output", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/lm_vocabulary.json"))
    ap.add_argument("--min-role-count", type=int, default=3)
    ap.add_argument("--word-vocab-size", type=int, default=8000)
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s - %(levelname)s - %(message)s")

    role_counts: Counter = Counter()
    word_counts: Counter = Counter()

    logger.info("scanning token corpus: %s", args.corpus)
    with open(args.corpus) as fh:
        for line in fh:
            try:
                ex = json.loads(line)
            except Exception:
                continue
            for t in ex["tokens"]:
                role_counts[t] += 1

    logger.info("  %d unique role-tokens (pre-filter)", len(role_counts))
    role_vocab: dict[str, int] = {}
    for tok, c in role_counts.most_common():
        if c < args.min_role_count:
            continue
        role_vocab[tok] = len(SPECIALS) + len(role_vocab)
    logger.info("  %d role-tokens kept (min_count=%d)", len(role_vocab), args.min_role_count)

    logger.info("scanning descriptions: %s", args.descriptions)
    if args.descriptions.exists():
        with open(args.descriptions) as fh:
            for line in fh:
                try:
                    row = json.loads(line)
                except Exception:
                    continue
                nl = row.get("natural_language", "") or ""
                for w in tokenize_description(nl):
                    word_counts[w] += 1
    logger.info("  %d unique words", len(word_counts))

    word_vocab: dict[str, int] = {}
    word_start = len(SPECIALS) + len(role_vocab)
    for word, c in word_counts.most_common(args.word_vocab_size):
        word_vocab[word] = word_start + len(word_vocab)

    vocab = {
        "version": 1,
        "specials": {tok: i for i, tok in enumerate(SPECIALS)},
        "role_tokens": role_vocab,
        "word_tokens": word_vocab,
        "vocab_size": len(SPECIALS) + len(role_vocab) + len(word_vocab),
        "role_count_floor": args.min_role_count,
        "word_vocab_cap": args.word_vocab_size,
    }
    args.output.write_text(json.dumps(vocab, ensure_ascii=False))
    logger.info("wrote vocabulary (size=%d) to %s", vocab["vocab_size"], args.output)
    print(f"vocab_size={vocab['vocab_size']}  specials={len(SPECIALS)}  roles={len(role_vocab)}  words={len(word_vocab)}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
