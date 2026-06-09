#!/usr/bin/env python3
"""
Phase 1 (PLASMID_TOKEN_MODEL_PLAN.md §7): baseline plasmid retrieval.

Builds two artefacts:
  1. Plasmid-level embeddings — for each (plasmid, rotation=0) row in the
     merged token_corpus.jsonl, embed its module/feature token stream as a
     sentence via a sentence-transformers model. This is the pre-LM baseline:
     no training, no transformer from scratch, just ANN over a semantic
     embedding of the plasmid's role-token sequence.
  2. HNSW index (cosine) over the plasmid embeddings for "find similar plasmids
     to this spec" retrieval.

Evaluation hook (run with --eval): held-out-description recall@k using the
VectorBuilder shorthand pairs (natural description → target plasmid tokens).

Outputs go under /var/data/plasmid_lm_corpus/retrieval/:
  - plasmid_embeddings.npy   (N × D float32)
  - plasmid_ids.txt          (one id per row, same order as .npy)
  - plasmid_index.hnsw       (HNSW binary)
  - retrieval_meta.json      (model name, dim, N, build stamp)
  - recall_report.json       (written when --eval is passed)

Usage:
  python build_retrieval_index.py \
      --corpus /var/data/plasmid_lm_corpus/token_corpus.jsonl \
      --outdir /var/data/plasmid_lm_corpus/retrieval/ \
      --model sentence-transformers/all-MiniLM-L6-v2

  # add --eval to compute recall@k against vb_shorthand_pairs.jsonl
"""
from __future__ import annotations
import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np

logger = logging.getLogger("retrieval")

# Tokens to drop when flattening to a sentence (header/boilerplate that's the
# same on every plasmid and adds no semantic signal).
NOISE_PREFIXES = ("<BOS>", "<EOS>", "<PAD>", "<UNK>",
                  "<TOPOLOGY:", "<LEN_BIN:", "<ROTATION_IDX:",
                  "<MOD_CLOSE>", "<UPSTREAM_REG>", "</UPSTREAM_REG>",
                  "<CDS_MOD>", "</CDS_MOD>",
                  "<DOWNSTREAM_REG>", "</DOWNSTREAM_REG>",
                  "<VB_CAS>", "</VB_CAS>", "<DRIVES>")


def tokens_to_sentence(tokens: list[str]) -> str:
    parts: list[str] = []
    for t in tokens:
        if any(t.startswith(p) for p in NOISE_PREFIXES):
            continue
        # Strip angle brackets, map colons to spaces, keep role + payload
        s = t.strip("<>").replace(":", " ").replace("_", " ")
        parts.append(s)
    return " ".join(parts)


def embed_plasmids(model_name: str, texts: list[str], batch_size: int = 64) -> np.ndarray:
    from sentence_transformers import SentenceTransformer
    logger.info("loading %s", model_name)
    model = SentenceTransformer(model_name)
    logger.info("embedding %d documents", len(texts))
    t0 = time.time()
    emb = model.encode(texts, batch_size=batch_size, show_progress_bar=True,
                       normalize_embeddings=True)
    logger.info("embedded in %.1fs; shape=%s", time.time() - t0, emb.shape)
    return np.asarray(emb, dtype=np.float32)


def build_hnsw(embs: np.ndarray, out_path: Path,
               M: int = 32, ef_construction: int = 200) -> None:
    import hnswlib
    dim = embs.shape[1]
    index = hnswlib.Index(space="cosine", dim=dim)
    index.init_index(max_elements=embs.shape[0],
                     ef_construction=ef_construction, M=M)
    index.add_items(embs, ids=np.arange(embs.shape[0]))
    index.set_ef(100)
    index.save_index(str(out_path))
    logger.info("wrote HNSW index → %s", out_path)


def query_hnsw(index_path: Path, dim: int, query_embs: np.ndarray, k: int) -> np.ndarray:
    import hnswlib
    idx = hnswlib.Index(space="cosine", dim=dim)
    idx.load_index(str(index_path))
    idx.set_ef(max(100, k * 4))
    labels, _dists = idx.knn_query(query_embs, k=k)
    return labels


def evaluate_recall(corpus_path: Path,
                    shorthand_path: Path,
                    out_meta: dict,
                    index_path: Path,
                    dim: int,
                    plasmid_ids: list[str],
                    id_to_tokens: dict[str, list[str]],
                    model_name: str) -> dict:
    """Recall@k against VectorBuilder shorthand descriptions."""
    from sentence_transformers import SentenceTransformer
    if not shorthand_path.exists():
        return {"skipped": "no shorthand file"}

    pairs = []
    with open(shorthand_path) as fh:
        for line in fh:
            d = json.loads(line)
            pairs.append((d["natural_language"], d.get("source_plasmid_id", "")))

    # Only keep pairs whose target plasmid is present in the index
    id_set = set(plasmid_ids)
    kept = [(nl, pid) for nl, pid in pairs if pid in id_set]
    logger.info("eval pairs: %d (of %d with id-in-index)", len(kept), len(pairs))
    if not kept:
        return {"skipped": "no target id overlap"}

    model = SentenceTransformer(model_name)
    q_texts = [nl for nl, _pid in kept]
    q_embs = np.asarray(model.encode(q_texts, batch_size=32, normalize_embeddings=True),
                        dtype=np.float32)
    labels = query_hnsw(index_path, dim, q_embs, k=10)

    id_to_rank = {pid: r for r, pid in enumerate(plasmid_ids)}
    hits_at = {1: 0, 3: 0, 5: 0, 10: 0}
    for i, (nl, gold_pid) in enumerate(kept):
        gold_rank = id_to_rank[gold_pid]
        retrieved = labels[i]
        for k in hits_at:
            if gold_rank in retrieved[:k]:
                hits_at[k] += 1
    n = len(kept)
    return {
        "pairs_evaluated": n,
        "recall_at_1": hits_at[1] / n,
        "recall_at_3": hits_at[3] / n,
        "recall_at_5": hits_at[5] / n,
        "recall_at_10": hits_at[10] / n,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--model", default="sentence-transformers/all-MiniLM-L6-v2")
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--rotation", type=int, default=0,
                    help="Use this rotation's tokens per plasmid (default canonical=0)")
    ap.add_argument("--eval", action="store_true")
    ap.add_argument("--shorthand", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/vb_shorthand_pairs.jsonl"))
    ap.add_argument("--descriptions", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/plasmid_descriptions.jsonl"),
                    help="Haiku-generated natural-language descriptions per plasmid; "
                         "preferred document source for retrieval since user prompts are "
                         "in human language. Falls back to role-token sentences when a "
                         "plasmid has no description entry.")
    ap.add_argument("--source", choices=["descriptions", "tokens"], default="descriptions",
                    help="Document source for embeddings. 'descriptions' (default) embeds "
                         "the four-style Haiku descriptions per plasmid; 'tokens' falls back "
                         "to flattened role-token sentences (the legacy v1 path).")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s - %(levelname)s - %(message)s")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Always read the token corpus to enumerate plasmid_ids + keep tokens as a
    # fallback document source.
    plasmid_ids: list[str] = []
    id_to_tokens: dict[str, list[str]] = {}
    with open(args.corpus) as fh:
        for line in fh:
            ex = json.loads(line)
            if ex.get("rotation_idx", 0) != args.rotation:
                continue
            pid = ex["plasmid_id"]
            if pid in id_to_tokens:
                continue
            id_to_tokens[pid] = ex["tokens"]
            plasmid_ids.append(pid)
    logger.info("%d unique plasmids extracted at rotation=%d", len(plasmid_ids), args.rotation)

    docs: list[str] = []
    n_from_desc, n_from_tokens = 0, 0
    if args.source == "descriptions":
        # Load Haiku description bundles keyed by plasmid_id. Each plasmid has
        # multiple description styles (user_request_short, user_request_functional,
        # lab_slack_question, methods_spec); concatenate all styles per plasmid
        # into a single document so each style can match a different user prompt.
        id_to_desc: dict[str, str] = {}
        with open(args.descriptions) as fh:
            for line in fh:
                d = json.loads(line)
                pid = d.get("plasmid_id")
                if not pid:
                    continue
                styles = d.get("descriptions") or []
                texts = [str(s.get("text", "")).strip() for s in styles if s.get("text")]
                if texts:
                    id_to_desc[pid] = " ".join(texts)
        logger.info("%d plasmids have Haiku descriptions", len(id_to_desc))

        for pid in plasmid_ids:
            if pid in id_to_desc:
                docs.append(id_to_desc[pid])
                n_from_desc += 1
            else:
                # Fall back to role-token sentence so we still have a document
                # (better than dropping the plasmid from the index).
                docs.append(tokens_to_sentence(id_to_tokens[pid]))
                n_from_tokens += 1
        logger.info("documents: %d from descriptions, %d from tokens (fallback)",
                    n_from_desc, n_from_tokens)
    else:
        for pid in plasmid_ids:
            docs.append(tokens_to_sentence(id_to_tokens[pid]))
        n_from_tokens = len(docs)
        logger.info("documents: %d from tokens (legacy mode)", n_from_tokens)

    embs = embed_plasmids(args.model, docs, batch_size=args.batch_size)
    np.save(args.outdir / "plasmid_embeddings.npy", embs)
    (args.outdir / "plasmid_ids.txt").write_text("\n".join(plasmid_ids) + "\n")

    build_hnsw(embs, args.outdir / "plasmid_index.hnsw")

    meta = {
        "model": args.model,
        "dim": int(embs.shape[1]),
        "n_plasmids": int(embs.shape[0]),
        "rotation": args.rotation,
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "document_source": args.source,
        "n_from_descriptions": n_from_desc,
        "n_from_tokens_fallback": n_from_tokens,
    }
    (args.outdir / "retrieval_meta.json").write_text(json.dumps(meta, indent=2))
    logger.info("wrote meta: %s", meta)

    if args.eval:
        report = evaluate_recall(args.corpus, args.shorthand, meta,
                                 args.outdir / "plasmid_index.hnsw",
                                 meta["dim"], plasmid_ids, id_to_tokens, args.model)
        (args.outdir / "recall_report.json").write_text(json.dumps(report, indent=2))
        logger.info("recall report: %s", json.dumps(report))

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
