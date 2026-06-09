#!/usr/bin/env python3
"""
Phase-2 plasmid LM trainer.

Trains a small decoder-only transformer on two mixed objectives:
  A. plasmid-only LM (~42 k examples from token_corpus.jsonl)
  B. description-to-plasmid (~31 k pairs from design_training_pairs.jsonl)

Each step draws from a random mix (default 1:1). Loss is masked so B only gets
credit on plasmid-token positions.

CPU-friendly defaults (~2-3 M params, small batch, short schedule). Scale up
via flags on a GPU:
  --d-model 256 --n-layer 6 --n-head 8 --batch-size 32 --max-steps 20000

Checkpoints + training log at --out-dir (default /var/data/plasmid_lm_corpus/lm_ckpt/).
"""
from __future__ import annotations
import argparse
import hashlib
import json
import logging
import random
import sys
import time
from pathlib import Path

import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader

sys.path.insert(0, str(Path(__file__).parent))
from plasmid_lm_model import ModelConfig, PlasmidLM
from plasmid_lm_dataset import (Vocab, PlasmidOnlyDataset, DescriptionPairDataset,
                                 collate)

logger = logging.getLogger("train")


def split_for(plasmid_id: str) -> str:
    h = int(hashlib.md5(plasmid_id.encode()).hexdigest(), 16) % 100
    if h < 5: return "val"
    if h < 10: return "test"
    return "train"


def make_loader(ds, batch_size, vocab, shuffle=True):
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle,
                      collate_fn=lambda b: collate(b, pad_id=vocab.pad_id),
                      num_workers=0, drop_last=False)


def loss_over_batch(model, batch, device) -> torch.Tensor:
    input_ids = batch["input_ids"].to(device)
    target_ids = batch["target_ids"].to(device)
    loss_mask = batch["loss_mask"].to(device)
    attn_mask = batch["attn_mask"].to(device)
    logits = model(input_ids, attn_mask=attn_mask)
    B, T, V = logits.shape
    logp = F.log_softmax(logits, dim=-1)
    # Gather target log-probs
    tgt_lp = logp.gather(-1, target_ids.unsqueeze(-1)).squeeze(-1)
    n = loss_mask.sum().clamp_min(1.0)
    return -(tgt_lp * loss_mask).sum() / n


def evaluate(model, loader, device, max_batches: int = 50) -> dict:
    model.eval()
    with torch.no_grad():
        losses = []
        for i, batch in enumerate(loader):
            if i >= max_batches: break
            losses.append(loss_over_batch(model, batch, device).item())
    mean_loss = sum(losses) / max(1, len(losses))
    return {"loss": mean_loss, "perplexity": float(torch.exp(torch.tensor(mean_loss)).item())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/token_corpus.jsonl"))
    ap.add_argument("--pairs", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/design_training_pairs.jsonl"))
    ap.add_argument("--vocab", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/lm_vocabulary.json"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/lm_ckpt"))
    # Model size
    ap.add_argument("--d-model", type=int, default=128)
    ap.add_argument("--n-layer", type=int, default=4)
    ap.add_argument("--n-head", type=int, default=4)
    ap.add_argument("--max-seq-len", type=int, default=384)
    ap.add_argument("--dropout", type=float, default=0.1)
    # Optim
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--max-steps", type=int, default=500)
    ap.add_argument("--warmup-steps", type=int, default=50)
    ap.add_argument("--log-every", type=int, default=20)
    ap.add_argument("--eval-every", type=int, default=100)
    ap.add_argument("--mix-ratio", type=float, default=0.5,
                    help="Probability of drawing a desc-pair batch (else plasmid-only)")
    ap.add_argument("--grad-clip", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s - %(levelname)s - %(message)s")

    torch.manual_seed(args.seed)
    random.seed(args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "args.json").write_text(json.dumps({k: str(v) for k, v in vars(args).items()}, indent=2))

    device = torch.device(args.device)
    vocab = Vocab.load(args.vocab)
    logger.info("vocab_size=%d (specials=%d roles=%d words=%d)",
                vocab.size, len(vocab.specials), len(vocab.roles), len(vocab.words))

    cfg = ModelConfig(vocab_size=vocab.size, d_model=args.d_model,
                      n_layer=args.n_layer, n_head=args.n_head,
                      max_seq_len=args.max_seq_len, dropout=args.dropout,
                      pad_id=vocab.pad_id)
    model = PlasmidLM(cfg).to(device)
    logger.info("model params: %.2fM", model.num_params() / 1e6)

    # Datasets
    logger.info("building datasets…")
    plas_train = PlasmidOnlyDataset(args.corpus, vocab, max_len=args.max_seq_len,
                                    split="train", split_fn=split_for)
    plas_val = PlasmidOnlyDataset(args.corpus, vocab, max_len=args.max_seq_len,
                                  split="val", split_fn=split_for)
    pairs_train = DescriptionPairDataset(args.pairs, vocab, max_len=args.max_seq_len,
                                         split="train")
    pairs_val = DescriptionPairDataset(args.pairs, vocab, max_len=args.max_seq_len,
                                       split="val")
    logger.info("  plas_train=%d plas_val=%d pairs_train=%d pairs_val=%d",
                len(plas_train), len(plas_val), len(pairs_train), len(pairs_val))

    plas_loader = make_loader(plas_train, args.batch_size, vocab)
    pairs_loader = make_loader(pairs_train, args.batch_size, vocab)
    plas_val_loader = make_loader(plas_val, args.batch_size, vocab, shuffle=False)
    pairs_val_loader = make_loader(pairs_val, args.batch_size, vocab, shuffle=False)

    optim = torch.optim.AdamW(model.parameters(), lr=args.lr, betas=(0.9, 0.95),
                              weight_decay=0.01)

    def lr_sched(step: int) -> float:
        if step < args.warmup_steps:
            return step / max(1, args.warmup_steps)
        progress = (step - args.warmup_steps) / max(1, args.max_steps - args.warmup_steps)
        return 0.5 * (1.0 + torch.cos(torch.tensor(progress * 3.14159)).item())

    def lr_set(step: int) -> None:
        for pg in optim.param_groups:
            pg["lr"] = args.lr * lr_sched(step)

    log_path = args.out_dir / "train.log"
    log_fh = open(log_path, "w", buffering=1)

    def log_line(d: dict) -> None:
        log_fh.write(json.dumps(d) + "\n")

    plas_iter = iter(plas_loader)
    pairs_iter = iter(pairs_loader)
    step = 0
    t0 = time.time()
    model.train()

    for step in range(1, args.max_steps + 1):
        lr_set(step)
        # Objective mixture
        use_pairs = (random.random() < args.mix_ratio)
        try:
            batch = next(pairs_iter) if use_pairs else next(plas_iter)
        except StopIteration:
            if use_pairs:
                pairs_iter = iter(pairs_loader)
                batch = next(pairs_iter)
            else:
                plas_iter = iter(plas_loader)
                batch = next(plas_iter)

        loss = loss_over_batch(model, batch, device)
        optim.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), args.grad_clip)
        optim.step()

        if step % args.log_every == 0 or step == 1:
            rate = step / max(1e-5, time.time() - t0)
            msg = {"step": step, "loss": loss.item(),
                   "lr": optim.param_groups[0]["lr"],
                   "obj": "pairs" if use_pairs else "plas",
                   "rate_steps_per_s": rate}
            logger.info("step=%d loss=%.4f lr=%.2e obj=%s %.2f step/s",
                        step, loss.item(), optim.param_groups[0]["lr"],
                        msg["obj"], rate)
            log_line(msg)

        if step % args.eval_every == 0:
            plas_m = evaluate(model, plas_val_loader, device, max_batches=25)
            pairs_m = evaluate(model, pairs_val_loader, device, max_batches=25)
            eval_line = {"step": step,
                         "eval_plas_val_loss": plas_m["loss"],
                         "eval_plas_val_ppl": plas_m["perplexity"],
                         "eval_pairs_val_loss": pairs_m["loss"],
                         "eval_pairs_val_ppl": pairs_m["perplexity"]}
            logger.info("EVAL step=%d plas_ppl=%.2f pairs_ppl=%.2f",
                        step, plas_m["perplexity"], pairs_m["perplexity"])
            log_line(eval_line)
            model.train()

    # Final checkpoint + eval
    ckpt_path = args.out_dir / "ckpt.pt"
    torch.save({"model": model.state_dict(),
                "cfg": cfg.__dict__,
                "step": step,
                "vocab_path": str(args.vocab)},
               ckpt_path)
    logger.info("saved checkpoint → %s (%.1f MB)", ckpt_path, ckpt_path.stat().st_size / 1e6)
    plas_m = evaluate(model, plas_val_loader, device, max_batches=50)
    pairs_m = evaluate(model, pairs_val_loader, device, max_batches=50)
    final = {"final": True, "step": step,
             "plas_val_loss": plas_m["loss"], "plas_val_ppl": plas_m["perplexity"],
             "pairs_val_loss": pairs_m["loss"], "pairs_val_ppl": pairs_m["perplexity"]}
    logger.info("FINAL plas_ppl=%.2f pairs_ppl=%.2f", plas_m["perplexity"], pairs_m["perplexity"])
    log_line(final)
    log_fh.close()
    print(json.dumps(final, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
