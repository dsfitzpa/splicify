"""
Small decoder-only transformer for the plasmid LM.

Two training objectives share the same decoder:
  Objective A (plasmid LM):   <BOS> + plasmid_tokens  ->  next-token CE
  Objective B (desc-to-plas): <DESC_START> + word_tokens + <DESC_END> + <SEP>
                              + plasmid_tokens         ->  next-token CE, loss
                              masked over the description prefix

One unified sequence model. At inference, feed the description prefix and let
it autoregressively generate the plasmid-token sequence.

Defaults are CPU-friendly (d_model=128, n_layer=4, n_head=4, ~2-3M params).
Scale up on GPU via the train script's CLI flags.
"""
from __future__ import annotations
import math
from dataclasses import dataclass

import torch
import torch.nn as nn
import torch.nn.functional as F


@dataclass
class ModelConfig:
    vocab_size: int
    d_model: int = 128
    n_layer: int = 4
    n_head: int = 4
    d_ff_mult: int = 4
    max_seq_len: int = 512
    dropout: float = 0.1
    pad_id: int = 0


class MultiHeadSelfAttention(nn.Module):
    def __init__(self, cfg: ModelConfig):
        super().__init__()
        assert cfg.d_model % cfg.n_head == 0
        self.cfg = cfg
        self.qkv = nn.Linear(cfg.d_model, 3 * cfg.d_model)
        self.out = nn.Linear(cfg.d_model, cfg.d_model)
        self.dropout = nn.Dropout(cfg.dropout)

    def forward(self, x: torch.Tensor, mask: torch.Tensor | None) -> torch.Tensor:
        B, T, D = x.shape
        H = self.cfg.n_head
        dh = D // H
        qkv = self.qkv(x).view(B, T, 3, H, dh).permute(2, 0, 3, 1, 4)  # 3,B,H,T,dh
        q, k, v = qkv[0], qkv[1], qkv[2]
        att = (q @ k.transpose(-2, -1)) / math.sqrt(dh)
        # Causal mask
        causal = torch.triu(torch.ones(T, T, device=x.device), diagonal=1).bool()
        att = att.masked_fill(causal, float("-inf"))
        if mask is not None:
            # mask: (B, T) with 1=keep, 0=pad → broadcast to (B, 1, 1, T)
            att = att.masked_fill(~mask.bool()[:, None, None, :], float("-inf"))
        att = F.softmax(att, dim=-1)
        att = self.dropout(att)
        out = (att @ v).transpose(1, 2).contiguous().view(B, T, D)
        return self.out(out)


class TransformerBlock(nn.Module):
    def __init__(self, cfg: ModelConfig):
        super().__init__()
        self.ln1 = nn.LayerNorm(cfg.d_model)
        self.attn = MultiHeadSelfAttention(cfg)
        self.ln2 = nn.LayerNorm(cfg.d_model)
        self.ff = nn.Sequential(
            nn.Linear(cfg.d_model, cfg.d_ff_mult * cfg.d_model),
            nn.GELU(),
            nn.Linear(cfg.d_ff_mult * cfg.d_model, cfg.d_model),
            nn.Dropout(cfg.dropout),
        )

    def forward(self, x: torch.Tensor, mask: torch.Tensor | None) -> torch.Tensor:
        x = x + self.attn(self.ln1(x), mask)
        x = x + self.ff(self.ln2(x))
        return x


class PlasmidLM(nn.Module):
    def __init__(self, cfg: ModelConfig):
        super().__init__()
        self.cfg = cfg
        self.token_emb = nn.Embedding(cfg.vocab_size, cfg.d_model, padding_idx=cfg.pad_id)
        self.pos_emb = nn.Embedding(cfg.max_seq_len, cfg.d_model)
        self.drop = nn.Dropout(cfg.dropout)
        self.blocks = nn.ModuleList([TransformerBlock(cfg) for _ in range(cfg.n_layer)])
        self.ln_final = nn.LayerNorm(cfg.d_model)
        self.head = nn.Linear(cfg.d_model, cfg.vocab_size, bias=False)
        # Tie weights
        self.head.weight = self.token_emb.weight
        self.apply(self._init)

    @staticmethod
    def _init(m: nn.Module) -> None:
        if isinstance(m, nn.Linear):
            nn.init.normal_(m.weight, mean=0.0, std=0.02)
            if m.bias is not None:
                nn.init.zeros_(m.bias)
        elif isinstance(m, nn.Embedding):
            nn.init.normal_(m.weight, mean=0.0, std=0.02)

    def forward(self, input_ids: torch.Tensor,
                attn_mask: torch.Tensor | None = None) -> torch.Tensor:
        B, T = input_ids.shape
        pos = torch.arange(T, device=input_ids.device).unsqueeze(0).expand(B, T)
        h = self.drop(self.token_emb(input_ids) + self.pos_emb(pos))
        for blk in self.blocks:
            h = blk(h, attn_mask)
        h = self.ln_final(h)
        logits = self.head(h)
        return logits

    @torch.no_grad()
    def generate(self, prefix_ids: torch.Tensor, max_new: int, eos_id: int,
                 temperature: float = 0.9, top_k: int | None = 40) -> torch.Tensor:
        self.eval()
        ids = prefix_ids
        for _ in range(max_new):
            if ids.shape[1] >= self.cfg.max_seq_len:
                break
            logits = self.forward(ids)[:, -1, :]
            logits = logits / max(1e-5, temperature)
            if top_k:
                v, _ = torch.topk(logits, top_k)
                logits[logits < v[:, [-1]]] = float("-inf")
            probs = F.softmax(logits, dim=-1)
            nxt = torch.multinomial(probs, num_samples=1)
            ids = torch.cat([ids, nxt], dim=1)
            if (nxt == eos_id).all():
                break
        return ids

    def num_params(self) -> int:
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
