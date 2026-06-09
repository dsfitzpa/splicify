"""
PyTorch Dataset + collator for the Phase-2 plasmid LM.

Two dataset views, same vocabulary:

1. PlasmidOnlyDataset   — iterates token_corpus.jsonl. Sample = BOS + tokens + EOS.
                           Target = shift-by-one; loss over every position.

2. DescriptionPairDataset — iterates design_training_pairs.jsonl. Sample =
                           DESC_START + word_tokens + DESC_END + SEP +
                           plasmid_tokens + EOS.
                           Loss masked over the description prefix so the model
                           only gets CE credit for predicting the plasmid tokens
                           from a desc conditioning.

Both datasets return dicts with keys input_ids, target_ids, loss_mask.
A single MixedCollator stacks them into padded batches.
"""
from __future__ import annotations
import json
import re
from dataclasses import dataclass
from pathlib import Path

import torch
from torch.utils.data import Dataset

WORD_SPLIT = re.compile(r"[A-Za-z0-9\-]+")


def tokenize_description(text: str) -> list[str]:
    return [w.lower() for w in WORD_SPLIT.findall(text)]


@dataclass
class Vocab:
    specials: dict[str, int]
    roles: dict[str, int]
    words: dict[str, int]
    size: int

    @classmethod
    def load(cls, path: Path) -> "Vocab":
        d = json.loads(Path(path).read_text())
        return cls(specials=d["specials"],
                   roles=d["role_tokens"],
                   words=d["word_tokens"],
                   size=d["vocab_size"])

    @property
    def pad_id(self) -> int: return self.specials["<PAD>"]
    @property
    def bos_id(self) -> int: return self.specials["<BOS>"]
    @property
    def eos_id(self) -> int: return self.specials["<EOS>"]
    @property
    def unk_id(self) -> int: return self.specials["<UNK>"]
    @property
    def desc_start(self) -> int: return self.specials["<DESC_START>"]
    @property
    def desc_end(self) -> int: return self.specials["<DESC_END>"]
    @property
    def sep_id(self) -> int: return self.specials["<SEP>"]

    def enc_role(self, tok: str) -> int:
        return self.roles.get(tok, self.unk_id)

    def enc_word(self, w: str) -> int:
        return self.words.get(w, self.unk_id)

    def encode_plasmid_tokens(self, tokens: list[str]) -> list[int]:
        return [self.enc_role(t) for t in tokens]

    def encode_description(self, text: str) -> list[int]:
        return [self.enc_word(w) for w in tokenize_description(text)]


class PlasmidOnlyDataset(Dataset):
    """Iterates per-example rows of token_corpus.jsonl; one (plasmid, rotation)
       per item. Optional split filter via plasmid_id hash."""

    def __init__(self, path: Path, vocab: Vocab, max_len: int = 512,
                 split: str | None = None,
                 split_fn=None):
        self.vocab = vocab
        self.max_len = max_len
        self.rows: list[dict] = []
        with open(path) as fh:
            for line in fh:
                ex = json.loads(line)
                if split and split_fn and split_fn(ex["plasmid_id"]) != split:
                    continue
                self.rows.append(ex)

    def __len__(self) -> int: return len(self.rows)

    def __getitem__(self, i: int) -> dict:
        ex = self.rows[i]
        ids = [self.vocab.bos_id] + self.vocab.encode_plasmid_tokens(ex["tokens"])[:self.max_len - 2] + [self.vocab.eos_id]
        input_ids = torch.tensor(ids[:-1], dtype=torch.long)
        target_ids = torch.tensor(ids[1:], dtype=torch.long)
        loss_mask = torch.ones_like(target_ids, dtype=torch.float)
        return {"input_ids": input_ids,
                "target_ids": target_ids,
                "loss_mask": loss_mask,
                "objective": "plasmid_only"}


class DescriptionPairDataset(Dataset):
    """Iterates design_training_pairs.jsonl. Optionally filter by split and/or
       source_type prefixes."""

    def __init__(self, path: Path, vocab: Vocab, max_len: int = 512,
                 split: str | None = None,
                 allow_source_types: set[str] | None = None):
        self.vocab = vocab
        self.max_len = max_len
        self.rows: list[dict] = []
        with open(path) as fh:
            for line in fh:
                row = json.loads(line)
                if split and row.get("split") != split:
                    continue
                if allow_source_types and row.get("source_type") not in allow_source_types:
                    continue
                if not row.get("target_strings"):
                    continue
                self.rows.append(row)

    def __len__(self) -> int: return len(self.rows)

    def __getitem__(self, i: int) -> dict:
        row = self.rows[i]
        v = self.vocab
        desc_ids = [v.desc_start] + v.encode_description(row["natural_language"])[:128] + [v.desc_end, v.sep_id]
        plas_ids = v.encode_plasmid_tokens(row["target_strings"])
        plas_ids = plas_ids[: max(1, self.max_len - len(desc_ids) - 1)] + [v.eos_id]
        ids = desc_ids + plas_ids
        input_ids = torch.tensor(ids[:-1], dtype=torch.long)
        target_ids = torch.tensor(ids[1:], dtype=torch.long)
        # Mask out loss on description prefix positions; keep loss on plasmid tokens
        loss_mask = torch.zeros_like(target_ids, dtype=torch.float)
        plas_offset = len(desc_ids) - 1  # first target position that lands on a plasmid token
        if plas_offset < loss_mask.shape[0]:
            loss_mask[plas_offset:] = 1.0
        return {"input_ids": input_ids,
                "target_ids": target_ids,
                "loss_mask": loss_mask,
                "objective": "desc_to_plasmid"}


def collate(batch: list[dict], pad_id: int = 0) -> dict:
    max_len = max(ex["input_ids"].shape[0] for ex in batch)
    B = len(batch)
    input_ids = torch.full((B, max_len), pad_id, dtype=torch.long)
    target_ids = torch.full((B, max_len), pad_id, dtype=torch.long)
    loss_mask = torch.zeros((B, max_len), dtype=torch.float)
    attn_mask = torch.zeros((B, max_len), dtype=torch.long)
    for i, ex in enumerate(batch):
        L = ex["input_ids"].shape[0]
        input_ids[i, :L] = ex["input_ids"]
        target_ids[i, :L] = ex["target_ids"]
        loss_mask[i, :L] = ex["loss_mask"]
        attn_mask[i, :L] = 1
    return {"input_ids": input_ids,
            "target_ids": target_ids,
            "loss_mask": loss_mask,
            "attn_mask": attn_mask}
