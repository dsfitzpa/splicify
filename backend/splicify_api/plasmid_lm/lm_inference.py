"""
PlasmidLMService — loads the trained plasmid LM checkpoint and exposes a simple
`generate(description) → {tokens, structure, ...}` call.

Used by `lm_router.py` (HTTP endpoint) and by the `lm_design` intent branch in
`chat.py`. Singleton loader: the checkpoint + vocab are read once per worker.
"""
from __future__ import annotations
import json
import logging
import threading
from pathlib import Path
from typing import Any

import torch

from .plasmid_lm_model import ModelConfig, PlasmidLM  # type: ignore
from .plasmid_lm_dataset import Vocab  # type: ignore

logger = logging.getLogger("plasmid_lm.inference")

_DEFAULT_CKPT = Path("/var/data/plasmid_lm_corpus/lm_ckpt/ckpt.pt")
_SERVICE_LOCK = threading.Lock()
_SERVICE: "PlasmidLMService | None" = None


class PlasmidLMService:
    def __init__(self, ckpt_path: Path = _DEFAULT_CKPT, device: str = "cpu"):
        logger.info("loading LM checkpoint from %s", ckpt_path)
        ckpt = torch.load(str(ckpt_path), map_location=device, weights_only=False)
        self.vocab = Vocab.load(Path(ckpt["vocab_path"]))
        cfg = ModelConfig(**ckpt["cfg"])
        self.model = PlasmidLM(cfg).to(device)
        self.model.load_state_dict(ckpt["model"])
        self.model.eval()
        self.device = device
        self.ckpt_path = str(ckpt_path)
        self.training_step = int(ckpt.get("step", 0))
        # Reverse vocab for token_id → token_string
        full = {**self.vocab.specials, **self.vocab.roles, **self.vocab.words}
        self._id2tok = {v: k for k, v in full.items()}
        logger.info("LM loaded: %.2fM params, vocab_size=%d, trained %d steps",
                    self.model.num_params() / 1e6, self.vocab.size,
                    self.training_step)

    @torch.no_grad()
    def generate(self,
                 description: str,
                 max_new_tokens: int = 256,
                 temperature: float = 0.8,
                 top_k: int = 40,
                 seed: int | None = None) -> dict[str, Any]:
        if seed is not None:
            torch.manual_seed(int(seed))
        v = self.vocab
        desc_ids = ([v.desc_start]
                    + v.encode_description(description)[:128]
                    + [v.desc_end, v.sep_id])
        prefix = torch.tensor([desc_ids], dtype=torch.long, device=self.device)
        out = self.model.generate(prefix, max_new=max_new_tokens,
                                  eos_id=v.eos_id,
                                  temperature=temperature, top_k=top_k)
        full_ids = out[0].tolist()
        full_toks = [self._id2tok.get(i, "<?>") for i in full_ids]

        # Split description prefix from generated plasmid tokens
        try:
            sep_pos = full_toks.index("<SEP>") + 1
        except ValueError:
            sep_pos = len(desc_ids)
        plasmid_toks = full_toks[sep_pos:]
        plasmid_ids = full_ids[sep_pos:]

        parsed = self._parse_structure(plasmid_toks)
        stop_reason = "eos" if plasmid_toks and plasmid_toks[-1] == "<EOS>" else "max_tokens"

        return {
            "description": description,
            "token_count": len(plasmid_toks),
            "token_strings": plasmid_toks,
            "token_ids": plasmid_ids,
            "header": parsed["header"],
            "modules": parsed["modules"],
            "interactions": parsed["interactions"],
            "cloning_features": parsed["cloning_features"],
            "orphan_features": parsed["orphan_features"],
            "stop_reason": stop_reason,
            "validation": self._validate(plasmid_toks),
            "model_info": {
                "training_step": self.training_step,
                "params_M": round(self.model.num_params() / 1e6, 2),
                "ckpt": self.ckpt_path,
                "note": "CPU smoke-training checkpoint — structural grammar learned; description-conditioning weak. Full GPU training pending.",
            },
        }

    @staticmethod
    def _parse_structure(toks: list[str]) -> dict:
        header: dict[str, str] = {}
        modules: list[dict] = []
        interactions: list[dict] = []
        cloning_features: list[dict] = []
        orphan_features: list[dict] = []

        current_module: dict | None = None
        depth = 0
        for t in toks:
            if not t.startswith("<"):
                continue
            body = t.strip("<>")
            parts = body.split(":", 2)
            role = parts[0]

            # Header tokens
            if role in ("BOS", "EOS", "PAD", "UNK"):
                continue
            if role in ("TOPOLOGY", "LEN_BIN", "HOST", "SOURCE", "ROTATION_IDX") and len(parts) > 1:
                header[role.lower()] = parts[1]
                continue
            if role == "MOD_OPEN":
                mt = parts[1] if len(parts) > 1 else "unknown"
                current_module = {"module_type": mt, "features": []}
                modules.append(current_module)
                depth += 1
                continue
            if role == "MOD_CLOSE":
                depth = max(0, depth - 1)
                current_module = None
                continue
            if role == "INT":
                interactions.append({"rule_id": parts[1] if len(parts) > 1 else "",
                                     "participant": parts[2] if len(parts) > 2 else None})
                continue
            if role == "CLN":
                sub = parts[1] if len(parts) > 1 else ""
                val = parts[2] if len(parts) > 2 else ""
                cloning_features.append({"subtype": sub, "value": val})
                continue
            # Feature tokens: live inside a module or as orphans
            feat = {"role": role,
                    "payload": parts[1] if len(parts) > 1 else None}
            if current_module is not None:
                current_module["features"].append(feat)
            else:
                orphan_features.append(feat)
        return {"header": header, "modules": modules,
                "interactions": interactions,
                "cloning_features": cloning_features,
                "orphan_features": orphan_features}

    @staticmethod
    def _validate(toks: list[str]) -> dict:
        stack: list[str] = []
        imbalance = False
        for t in toks:
            if t.startswith("<MOD_OPEN:"):
                stack.append(t)
            elif t == "<MOD_CLOSE>":
                if not stack:
                    imbalance = True
                else:
                    stack.pop()
        return {
            "bracket_balance": not imbalance and not stack,
            "unclosed_modules": len(stack),
            "has_eos": toks and toks[-1] == "<EOS>",
        }


def get_lm_service(ckpt_path: Path = _DEFAULT_CKPT,
                   device: str = "cpu") -> PlasmidLMService:
    global _SERVICE
    if _SERVICE is not None:
        return _SERVICE
    with _SERVICE_LOCK:
        if _SERVICE is None:
            _SERVICE = PlasmidLMService(ckpt_path=ckpt_path, device=device)
        return _SERVICE


def format_lm_reply(result: dict) -> str:
    """Markdown summary of a generation result for chat display."""
    hdr = result.get("header", {})
    mods = result.get("modules", [])
    ints = result.get("interactions", [])
    clns = result.get("cloning_features", [])
    val = result.get("validation", {})
    info = result.get("model_info", {})

    lines: list[str] = []
    lines.append(f"**Generated plasmid design** (from LM checkpoint at step {info.get('training_step','?')})")
    lines.append("")
    lines.append("> ⚠ *Experimental — CPU smoke-training checkpoint (~2 M params). Structural grammar is learned; description-conditioning is weak. Full GPU training pending.*")
    lines.append("")
    if hdr:
        bits = []
        if "host" in hdr: bits.append(f"host `{hdr['host']}`")
        if "topology" in hdr: bits.append(f"topology `{hdr['topology']}`")
        if "len_bin" in hdr: bits.append(f"length `{hdr['len_bin']}`")
        if "source" in hdr: bits.append(f"source `{hdr['source']}`")
        lines.append("**Header:** " + " · ".join(bits))
        lines.append("")
    lines.append(f"**Modules ({len(mods)}):**")
    for m in mods[:20]:
        feats = m.get("features", [])
        feat_str = ", ".join(f"`{f['role']}:{f.get('payload') or ''}`" for f in feats[:6])
        more = f" … (+{len(feats)-6} more)" if len(feats) > 6 else ""
        lines.append(f"- `{m['module_type']}` — {feat_str}{more}")
    if len(mods) > 20:
        lines.append(f"… and {len(mods)-20} more modules")
    lines.append("")
    if ints:
        lines.append(f"**Interactions emitted ({len(ints)}):**")
        seen = set()
        for ix in ints:
            rid = ix.get("rule_id", "")
            if rid in seen: continue
            seen.add(rid)
            lines.append(f"- `{rid}`")
        lines.append("")
    if clns:
        uniq = sorted(set(c["value"] for c in clns if c.get("subtype") == "unique_cutter"))
        noncut = sorted(set(c["value"] for c in clns if c.get("subtype") == "nonCutter"))
        if uniq:
            lines.append(f"**Unique cutters** ({len(uniq)}): " + ", ".join(f"`{e}`" for e in uniq[:12]) + ("…" if len(uniq) > 12 else ""))
        if noncut:
            lines.append(f"**Non-cutters** ({len(noncut)}): " + ", ".join(f"`{e}`" for e in noncut[:12]) + ("…" if len(noncut) > 12 else ""))
        lines.append("")
    v_bits = []
    v_bits.append("✓ balanced" if val.get("bracket_balance") else f"✗ {val.get('unclosed_modules',0)} unclosed")
    v_bits.append("✓ EOS" if val.get("has_eos") else "✗ no EOS (hit max_tokens)")
    lines.append("**Structural validation:** " + " · ".join(v_bits))
    lines.append(f"**Token count:** {result.get('token_count')}  ·  **Stop reason:** {result.get('stop_reason')}")
    return "\n".join(lines)
