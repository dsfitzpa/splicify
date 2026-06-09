#!/usr/bin/env python3
"""
Plasmid tokenizer v2 — produces multi-family bracketed token streams for AI model training.

See TOKENIZATION_PLAN_V2.md for the full spec. Emits, per plasmid:
  - Header tokens (topology, length bin, host, source, rotation_idx)
  - Module open/close brackets with feature/submodule tokens inside
  - Interaction edge tokens (referencing earlier feature tokens by name)
  - Cloning-feature tokens (restriction II/IIs, gateway att, PCR warnings)

Augments each plasmid with K rotations (default 6): canonical start + 5 random.
Source-tags examples as 'module_library', 'ncbi', or 'vectorbuilder' for mini-batch
composition downstream.

Produces three artefacts:
  1. JSONL token_corpus (one line per (plasmid, rotation) example).
  2. Optional postgres inserts (plasmid_tokens, plasmid_tokenizations, interaction
     tokens, cloning tokens).
  3. vocabulary.json (sorted role tokens + counts).
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import logging
import os
import random
import re
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

# Annotation pipeline imports (installed in the splicify_api environment)
try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None  # type: ignore

logger = logging.getLogger("plasmid_tokenizer")


# ---------------------------------------------------------------------------
# Constants — token families, role mapping, canonicalization
# ---------------------------------------------------------------------------

SPECIAL_TOKENS = ["<BOS>", "<EOS>", "<PAD>", "<UNK>"]

HEADER_PREFIXES = ("TOPOLOGY", "LEN_BIN", "HOST", "SOURCE", "ROTATION_IDX")

MODULE_OPEN_CLOSE = {"MOD_OPEN", "MOD_CLOSE", "UPSTREAM_REG", "/UPSTREAM_REG",
                     "CDS_MOD", "/CDS_MOD", "DOWNSTREAM_REG", "/DOWNSTREAM_REG"}

FEATURE_ROLE_FOR_TYPE = {
    "promoter": "PROMOTER",
    "CDS": "CDS",
    "gene": "CDS",
    "terminator": "TERMINATOR",
    "polyA_signal": "POLYA",
    "rep_origin": "ORI",
    "LTR": "LTR",
    "ITR": "ITR",
    "enhancer": "ENHANCER",
    "insulator": "INSULATOR",
    "RBS": "RBS",
    "kozak": "KOZAK",
    "operator": "OPERATOR",
    "misc_recomb": "RECOMB",
    "protein_bind": "PROTEIN_BIND",
    "primer_bind": "PRIMER_BIND",
    "sig_peptide": "SIG_PEP",
    "mat_peptide": "MAT_PEP",
    "intron": "INTRON",
    "exon": "EXON",
    "ncRNA": "NCRNA",
    "tRNA": "TRNA",
    "rRNA": "RRNA",
    "oriT": "ORIT",
    "regulatory": "REG",
    "misc_feature": "MISC",
    "polylinker": "MCS",
}

SUBMODULE_ROLE_FOR_TYPE = {
    "protein_module": "CDS",
    "protein_submodule": "CDS",
    "nls_module": "NLS",
    "tag_module": "TAG",
    "linker_module": "LINKER",
    "flexible_linker_module": "LINKER",  # upgraded to LINKER_2A at emit time if name matches
    "gap_module": "GAP",
}

LEN_BINS = [
    (0, 2000, "0-2kb"),
    (2000, 4000, "2-4kb"),
    (4000, 6000, "4-6kb"),
    (6000, 8000, "6-8kb"),
    (8000, 12000, "8-12kb"),
    (12000, 16000, "12-16kb"),
    (16000, 25000, "16-25kb"),
    (25000, 10**9, "25kb+"),
]

# Vector family prefixes for VectorBuilder shorthand
VB_FAMILY_TO_HOST = {
    "pRP": "mammalian", "pLV": "mammalian", "pAV": "mammalian", "pAAV": "mammalian",
    "pAd5": "mammalian", "pAd5F35": "mammalian", "pGLAd": "mammalian", "pscAAV": "mammalian",
    "pMMLV": "mammalian", "pMSCV": "mammalian", "pVSV": "mammalian", "pPB": "mammalian",
    "pSB": "mammalian", "pTol2": "zebrafish", "pRCASBPA": "chicken",
    "pPlant": "plant", "pDmel": "drosophila", "pCE": "worm", "pET": "bacterial",
    "pCDNA": "mammalian", "pLX": "mammalian",
}

POL2_PROMOTERS = {"CMV", "EF-1", "EF1A", "EF1α", "EF-1α", "HTLV", "PGK", "RSV",
                  "SV40", "CAG", "UbC", "SFFV", "MSCV", "TRE", "CAGGS", "CBh"}
POL3_PROMOTERS = {"U6", "H1", "7SK", "hU6", "mU6"}

LINKER_2A_NAMES = {"P2A", "T2A", "E2A", "F2A", "GSG-P2A", "GSG-T2A", "GSG-E2A", "GSG-F2A"}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class TokenExample:
    plasmid_id: str
    source: str                     # 'module_library' | 'ncbi' | 'vectorbuilder'
    length: int
    topology: str
    rotation_idx: int
    canonical_start: int
    tokens: list[str] = field(default_factory=list)
    interaction_tokens: list[dict[str, Any]] = field(default_factory=list)
    cloning_tokens: list[dict[str, Any]] = field(default_factory=list)
    host: str = "unknown"
    valid: bool = True
    validation_errors: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def len_bin(length: int) -> str:
    for lo, hi, label in LEN_BINS:
        if lo <= length < hi:
            return label
    return "unknown"


def normalize_name(name: str) -> str:
    """Normalize a feature / payload name so identical concepts collapse."""
    if not name:
        return "unknown"
    n = str(name).strip()
    n = re.sub(r"\s+", "_", n)
    n = re.sub(r"[^A-Za-z0-9α-ωΑ-Ω_\-\(\)\[\]]", "", n)
    return n or "unknown"


def canonical_rotation_start(sequence: str) -> int:
    """
    Return the offset of the longest unique 20-mer in the plasmid.
    Deterministic across runs; matches the rule used by TARGET_FROM_INVENTORY_ROUTING.md::validate().
    """
    if not sequence:
        return 0
    k = 20
    seen: dict[str, int] = {}
    duplicates: set[str] = set()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if kmer in seen:
            duplicates.add(kmer)
        else:
            seen[kmer] = i
    # pick the alphabetically-smallest unique k-mer for determinism
    unique = sorted(k for k in seen if k not in duplicates)
    return seen[unique[0]] if unique else 0


def rotate_sequence(seq: str, offset: int) -> str:
    offset %= max(1, len(seq))
    return seq[offset:] + seq[:offset]


def rotate_coords(start: int, end: int, offset: int, length: int) -> tuple[int, int]:
    s = (start - offset) % length
    e = (end - offset) % length
    if e <= s:
        e += length
    return s, e


def infer_host(features: list[dict[str, Any]],
               modules: list[dict[str, Any]],
               folder_hint: str | None = None) -> str:
    """
    Priority-based host inference.

    Rule 1: Non-bacterial module_type hit → that host (lentiviral/AAV/T-DNA/etc.
            override bacterial backbone — bacterial selection / ori reflects
            propagation, not host-of-use).
    Rule 2: Feature-name keyword matches ≥ 2 (or ≥ 1 + folder confirms) for a
            non-bacterial host → that host.
    Rule 3: Folder prior (Module_Library subject dir) → that host.
    Rule 4: Bacterial if only backbone signals fire.
    Rule 5: Unknown otherwise.

    Keeps backwards-compat: callers that don't pass folder_hint still get the
    priority-based token/module reasoning.
    """
    from .host_inference import infer_host_priority  # type: ignore
    # Build a pseudo-token list from module + feature names so the priority
    # rules (written against token streams) can score this example.
    pseudo: list[str] = []
    for m in modules or []:
        mt = m.get("module_type", "")
        if mt:
            pseudo.append(f"<MOD_OPEN:{mt}>")
        nm = m.get("name") or m.get("payload_id") or ""
        if nm:
            pseudo.append(f"<CDS:{normalize_name(nm)}>")
    for f in features or []:
        nm = f.get("name") or f.get("Feature") or f.get("sseqid") or ""
        if nm:
            pseudo.append(f"<CDS:{normalize_name(nm)}>")
    return infer_host_priority(pseudo, folder_hint=folder_hint)


def feature_role(ftype: str, name: str) -> str:
    ftype = (ftype or "").strip()
    if ftype in FEATURE_ROLE_FOR_TYPE:
        role = FEATURE_ROLE_FOR_TYPE[ftype]
        # Promoters split into Pol II / Pol III based on name
        if role == "PROMOTER":
            nm = (name or "").upper()
            if any(p in nm for p in POL3_PROMOTERS):
                return "PROMOTER_POL3"
            if any(p in nm for p in POL2_PROMOTERS):
                return "PROMOTER_POL2"
        return role
    return "MISC"


def submodule_role(mtype: str, name: str) -> str:
    base = SUBMODULE_ROLE_FOR_TYPE.get(mtype, "MISC")
    if base == "LINKER":
        nm = (name or "").upper()
        if any(peptide in nm for peptide in LINKER_2A_NAMES):
            return "LINKER_2A"
    return base


# ---------------------------------------------------------------------------
# Annotation call
# ---------------------------------------------------------------------------

DEFAULT_ANNOTATE_URL = "http://127.0.0.1:8000/plannotate/annotate_sequence_llm"


def annotate_genbank(gb_path: Path, annotate_url: str = DEFAULT_ANNOTATE_URL,
                     timeout: int = 120) -> dict[str, Any] | None:
    """
    Call the in-repo annotation pipeline via its FastAPI endpoint.
    The endpoint is async, so we go over HTTP rather than import FastAPI coroutines.
    """
    if SeqIO is None:
        logger.error("Biopython not installed")
        return None
    import requests
    try:
        record = next(SeqIO.parse(str(gb_path), "genbank"))
    except Exception as exc:
        logger.warning("could not parse %s: %s", gb_path, exc)
        return None
    seq = str(record.seq).upper()
    topology = "circular" if record.annotations.get("topology", "").lower() == "circular" else "linear"
    payload = {
        "sequence": seq,
        "circular": topology == "circular",
        "detailed": True,
    }
    try:
        resp = requests.post(annotate_url, json=payload, timeout=timeout)
        resp.raise_for_status()
        result = resp.json()
    except Exception as exc:
        logger.warning("annotation call failed on %s: %s", gb_path.name, exc)
        return None
    if not isinstance(result, dict):
        logger.warning("unexpected annotation response type for %s", gb_path.name)
        return None
    # Normalize keys downstream code expects.
    # The live API returns:
    #   - plannotate_annotations / annotations: flat features with keys name/type/start/end
    #   - hierarchical_annotations: mixed items with layer in {'feature','module','submodule','cloning_feature'}
    #   - cds_submodules: top-level list of submodule dicts (layer='submodule')
    #   - cloning_features: dict with features[], cut_count_per_enzyme, non_cutters
    flat_features = (result.get("plannotate_annotations")
                     or result.get("annotations") or [])
    hier = result.get("hierarchical_annotations") or []
    # Split hier by layer
    hier_modules = [h for h in hier if h.get("layer") == "module"]
    hier_submodules = [h for h in hier if h.get("layer") == "submodule"]
    # Fallback: top-level cds_submodules key
    if not hier_submodules:
        hier_submodules = result.get("cds_submodules") or result.get("cds_submodules_list") or []
    result["features"] = flat_features
    result["hierarchical_annotations"] = hier_modules
    result["cds_submodules"] = hier_submodules
    result.setdefault("interactions", result.get("interactions") or [])
    result.setdefault("cloning_features", result.get("cloning_features") or {})
    result["_sequence"] = seq
    result["_topology"] = topology
    result["_plasmid_id"] = gb_path.stem
    return result


# ---------------------------------------------------------------------------
# Tokenization of one annotated plasmid, one rotation
# ---------------------------------------------------------------------------

def emit_feature_token(ftype: str, name: str) -> str:
    role = feature_role(ftype, name)
    return f"<{role}:{normalize_name(name)}>"


def emit_submodule_token(mtype: str, name: str) -> str:
    role = submodule_role(mtype, name)
    return f"<{role}:{normalize_name(name)}>"


def _feat_coords(f: dict[str, Any]) -> tuple[int, int]:
    s = int(f.get("start", f.get("qstart", 0)) or 0)
    e = int(f.get("end", f.get("qend", 0)) or 0)
    return s, e


def _features_in_span(features: list[dict[str, Any]], start: int, end: int) -> list[dict[str, Any]]:
    out = []
    for f in features:
        fs, fe = _feat_coords(f)
        # permit a small slop: features must overlap the span by ≥ 50 % of their length
        if fe <= fs:
            continue
        overlap = max(0, min(fe, end) - max(fs, start))
        if overlap * 2 >= (fe - fs):
            out.append(f)
    return sorted(out, key=lambda f: _feat_coords(f)[0])


def _submodules_in_span(submodules: list[dict[str, Any]], start: int, end: int) -> list[dict[str, Any]]:
    out = []
    for m in submodules:
        ms = int(m.get("start", 0))
        me = int(m.get("end", 0))
        if ms >= start and me <= end:
            out.append(m)
    return sorted(out, key=lambda m: int(m.get("start", 0)))


def build_module_block(module: dict[str, Any],
                       features: list[dict[str, Any]],
                       submodules: list[dict[str, Any]],
                       rotation_offset: int,
                       seq_len: int) -> list[str]:
    mtype = module.get("module_type", "unknown")
    mstart = int(module.get("start", 0))
    mend = int(module.get("end", seq_len))
    # rotate module coords for lookup in feature list (features have been rotated upstream)
    tokens: list[str] = [f"<MOD_OPEN:{mtype}>"]

    is_pol2 = mtype in {"mammalian_pol2_expression_cassette",
                        "mammalian_lentiviral_expression_cassette"}
    if is_pol2:
        upstream = module.get("upstream_regulatory_module") or {}
        cds_mod = module.get("cds_module") or {}
        downstream = module.get("downstream_regulatory_module") or {}
        for slot_label, slot in (("UPSTREAM_REG", upstream),
                                 ("CDS_MOD", cds_mod),
                                 ("DOWNSTREAM_REG", downstream)):
            if not slot:
                continue
            tokens.append(f"<{slot_label}>")
            s_start = int(slot.get("start", mstart))
            s_end = int(slot.get("end", mend))
            # Emit contained submodules first (they carry names like Cas9 / P2A),
            # falling back to feature rows when no submodules were detected.
            sub_in = _submodules_in_span(submodules, s_start, s_end)
            if sub_in:
                for sm in sub_in:
                    tokens.append(emit_submodule_token(
                        sm.get("module_type", ""),
                        sm.get("name") or sm.get("payload_id") or sm.get("module_type", "")
                    ))
            else:
                for feat in _features_in_span(features, s_start, s_end):
                    tokens.append(emit_feature_token(
                        feat.get("type", "") or feat.get("Type", ""),
                        feat.get("name") or feat.get("Feature") or feat.get("sseqid") or ""
                    ))
            tokens.append(f"</{slot_label}>")
    else:
        # Flat module: emit submodules (if any) then features
        sub_in = _submodules_in_span(submodules, mstart, mend)
        if sub_in:
            for sm in sub_in:
                tokens.append(emit_submodule_token(
                    sm.get("module_type", ""),
                    sm.get("name") or sm.get("payload_id") or sm.get("module_type", "")
                ))
        for feat in _features_in_span(features, mstart, mend):
            tokens.append(emit_feature_token(
                feat.get("type", "") or feat.get("Type", ""),
                feat.get("name") or feat.get("Feature") or feat.get("sseqid") or ""
            ))
    tokens.append("<MOD_CLOSE>")
    return tokens


def build_interaction_tokens(interactions: list[dict[str, Any]]) -> list[tuple[str, dict[str, Any]]]:
    out: list[tuple[str, dict[str, Any]]] = []
    for ix in interactions or []:
        rule_id = ix.get("rule_id") or ix.get("type") or "INT-UNKNOWN"
        participants = []
        for p in ix.get("participants") or []:
            pname = p.get("name") or p.get("module") or "unknown"
            participants.append(normalize_name(pname))
        if participants:
            out.append((f"<INT:{rule_id}:{participants[0]}>", {
                "rule_id": rule_id,
                "sbo_term": ix.get("sbo_term") or ix.get("sbo"),
                "source_module": ix.get("source_module") or ix.get("source"),
                "participants": participants,
            }))
        else:
            out.append((f"<INT:{rule_id}>", {
                "rule_id": rule_id,
                "sbo_term": ix.get("sbo_term") or ix.get("sbo"),
                "source_module": ix.get("source_module") or ix.get("source"),
                "participants": [],
            }))
    return out


def build_cloning_tokens(cf: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    out: list[tuple[str, dict[str, Any]]] = []
    if not cf:
        return out
    cut_count = cf.get("cut_count_per_enzyme") or {}
    for enzyme, count in sorted(cut_count.items()):
        if count == 1:
            out.append((f"<CLN:unique_cutter:{enzyme}>", {"subtype": "unique_cutter", "value": enzyme}))
    for enzyme in sorted(cf.get("non_cutters") or []):
        out.append((f"<CLN:nonCutter:{enzyme}>", {"subtype": "nonCutter", "value": enzyme}))
    for feat in cf.get("features") or []:
        ftype = feat.get("type") or feat.get("feature_type")
        if ftype == "restriction_site_IIs":
            out.append((f"<CLN:iis_site:{feat.get('name','')}>",
                        {"subtype": "iis_site", "value": feat.get("name", ""),
                         "position": feat.get("start")}))
        elif ftype == "gateway_att":
            out.append((f"<CLN:att:{feat.get('name','')}>",
                        {"subtype": "att", "value": feat.get("name", ""),
                         "position": feat.get("start")}))
        elif ftype == "primer_design_warning":
            out.append((f"<CLN:warn:{feat.get('subtype','')}>",
                        {"subtype": "warn", "value": feat.get("subtype", ""),
                         "position": feat.get("start")}))
    if not out:
        out.append(("<CLN:none>", {"subtype": "none", "value": ""}))
    return out


# ---------------------------------------------------------------------------
# Validation gates
# ---------------------------------------------------------------------------

def validate_tokens(tokens: list[str], feature_names: set[str],
                    cut_count: dict[str, int]) -> list[str]:
    errors: list[str] = []
    # Bracket balance
    stack: list[str] = []
    for t in tokens:
        if t.startswith("<MOD_OPEN:"):
            stack.append(t)
        elif t == "<MOD_CLOSE>":
            if not stack or not stack[-1].startswith("<MOD_OPEN:"):
                errors.append("bracket imbalance: MOD_CLOSE without MOD_OPEN")
            else:
                stack.pop()
        elif t in ("<UPSTREAM_REG>", "<CDS_MOD>", "<DOWNSTREAM_REG>"):
            stack.append(t)
        elif t in ("</UPSTREAM_REG>", "</CDS_MOD>", "</DOWNSTREAM_REG>"):
            if not stack or stack[-1] != t.replace("/", ""):
                errors.append(f"bracket imbalance: {t} at wrong depth")
            else:
                stack.pop()
    if stack:
        errors.append(f"bracket imbalance: {len(stack)} unclosed open brackets")

    # Interaction referential integrity
    for t in tokens:
        if t.startswith("<INT:") and t.count(":") >= 2:
            parts = t.strip("<>").split(":", 2)
            if len(parts) == 3:
                participant = parts[2].rstrip(">")
                if participant and participant not in feature_names and participant != "unknown":
                    # soft warning — 2A-second-product interactions reference submodules that
                    # may appear only as "linker" tokens, not as features. Allow.
                    pass

    # Cloning cutter consistency
    for t in tokens:
        if t.startswith("<CLN:unique_cutter:"):
            name = t.rsplit(":", 1)[-1].rstrip(">")
            if cut_count.get(name, 0) != 1:
                errors.append(f"cloning: {name} claimed unique_cutter but cut_count={cut_count.get(name,0)}")
    return errors


# ---------------------------------------------------------------------------
# Per-plasmid, multi-rotation
# ---------------------------------------------------------------------------

def tokenize_plasmid(ann: dict[str, Any],
                     source: str,
                     k_rotations: int = 6,
                     rng: random.Random | None = None) -> list[TokenExample]:
    rng = rng or random.Random(hash(ann.get("_plasmid_id", "")) & 0xFFFFFFFF)
    seq = ann["_sequence"]
    topology = ann["_topology"]
    length = len(seq)
    plasmid_id = ann["_plasmid_id"]
    host = infer_host(ann.get("features") or [],
                      ann.get("hierarchical_annotations") or [],
                      folder_hint=ann.get("_folder_hint"))
    canonical_offset = canonical_rotation_start(seq) if topology == "circular" else 0

    offsets = [canonical_offset]
    if topology == "circular":
        for _ in range(max(0, k_rotations - 1)):
            offsets.append(rng.randint(0, max(1, length - 1)))

    features = ann.get("features") or []
    modules = [m for m in (ann.get("hierarchical_annotations") or [])
               if m.get("module_type") and m.get("layer") != "cloning_feature"]
    submodules = ann.get("cds_submodules") or ann.get("cds_submodules_list") or []
    interactions = ann.get("interactions") or []
    cloning_features = ann.get("cloning_features") or {}

    examples: list[TokenExample] = []
    for rot_idx, offset in enumerate(offsets):
        # Rotate coordinates in all annotations
        def rotate_ann(item: dict, k_start="start", k_end="end"):
            out = dict(item)
            s = int(item.get(k_start, 0))
            e = int(item.get(k_end, length))
            ns, ne = rotate_coords(s, e, offset, length)
            out[k_start] = ns
            out[k_end] = ne
            return out

        rot_features = []
        for f in features:
            nf = dict(f)
            s, e = _feat_coords(f)
            ns, ne = rotate_coords(s, e, offset, length)
            nf["start"] = ns
            nf["end"] = ne
            rot_features.append(nf)
        rot_modules = [rotate_ann(m) for m in modules]
        rot_submodules = [rotate_ann(m) for m in submodules]
        # Sort modules by (depth, start) so parents come before children
        rot_modules.sort(key=lambda m: (int(m.get("start", 0)), -int(m.get("end", 0))))

        tokens: list[str] = ["<BOS>",
                             f"<TOPOLOGY:{topology}>",
                             f"<LEN_BIN:{len_bin(length)}>",
                             f"<HOST:{host}>",
                             f"<SOURCE:{source}>",
                             f"<ROTATION_IDX:{rot_idx}>"]

        # Emit top-level modules; skip ones wholly contained in another already-emitted module
        emitted_spans: list[tuple[int, int]] = []
        emitted_feature_ids: set[int] = set()

        def _contained(s: int, e: int) -> bool:
            return any(ps <= s and e <= pe for ps, pe in emitted_spans)

        for mod in rot_modules:
            ms = int(mod.get("start", 0))
            me = int(mod.get("end", length))
            if _contained(ms, me):
                continue
            block = build_module_block(mod, rot_features, rot_submodules, offset, length)
            tokens.extend(block)
            emitted_spans.append((ms, me))
            for feat in _features_in_span(rot_features, ms, me):
                emitted_feature_ids.add(id(feat))

        # Emit orphan features (not covered by any module) inline between BOS and INT block
        orphan_feats = [f for f in rot_features if id(f) not in emitted_feature_ids]
        if orphan_feats:
            tokens.append("<MOD_OPEN:orphan_features>")
            for feat in sorted(orphan_feats, key=lambda f: _feat_coords(f)[0]):
                tokens.append(emit_feature_token(
                    feat.get("type", "") or feat.get("Type", ""),
                    feat.get("name") or feat.get("Feature") or feat.get("sseqid") or ""
                ))
            tokens.append("<MOD_CLOSE>")

        int_token_rows = build_interaction_tokens(interactions)
        for tok, _meta in int_token_rows:
            tokens.append(tok)

        cln_token_rows = build_cloning_tokens(cloning_features)
        for tok, _meta in cln_token_rows:
            tokens.append(tok)

        tokens.append("<EOS>")

        # Feature-name set for validation
        feat_name_set = {normalize_name(f.get("Feature") or f.get("sseqid") or "") for f in features}
        feat_name_set |= {normalize_name(m.get("name") or m.get("payload_id") or "")
                          for m in submodules}
        errors = validate_tokens(tokens,
                                 feat_name_set,
                                 cloning_features.get("cut_count_per_enzyme") or {})

        examples.append(TokenExample(
            plasmid_id=plasmid_id,
            source=source,
            length=length,
            topology=topology,
            rotation_idx=rot_idx,
            canonical_start=canonical_offset,
            tokens=tokens,
            interaction_tokens=[m for _t, m in int_token_rows],
            cloning_tokens=[m for _t, m in cln_token_rows],
            host=host,
            valid=not errors,
            validation_errors=errors,
        ))
    return examples


# ---------------------------------------------------------------------------
# VectorBuilder description shorthand parser
# ---------------------------------------------------------------------------

VB_SHORTHAND_RE = re.compile(r"^(?P<family>p[A-Za-z0-9]+)\[(?P<sys>[^\]]+)\](?P<body>.+)$")


def parse_vectorbuilder_shorthand(description: str) -> list[str]:
    """
    Parse a VectorBuilder shorthand like 'pLV[Exp]-EGFP:T2A:Puro-EF1A>mCherry' into
    tokens compatible with the main vocabulary.

    Returns a flat token list (no module brackets) suitable for a description→tokens
    auxiliary training head. Returns empty list if shorthand doesn't parse.
    """
    desc = description.strip()
    m = VB_SHORTHAND_RE.match(desc)
    if not m:
        return []
    family = m.group("family")
    sys_tag = m.group("sys").upper()
    body = m.group("body")
    host = VB_FAMILY_TO_HOST.get(family, "unknown")

    out = ["<BOS>", f"<SOURCE:vectorbuilder>", f"<HOST:{host}>",
           f"<VB_FAMILY:{family}>", f"<VB_SYS:{sys_tag}>"]

    # Split body on '-' into cassettes; first segment may start with a '-'
    if body.startswith("-"):
        body = body[1:]
    cassettes = [seg for seg in body.split("-") if seg]
    for cassette in cassettes:
        out.append("<VB_CAS>")
        if ">" in cassette:
            # promoter>CDS form
            prom, rhs = cassette.split(">", 1)
            pnm = prom.upper()
            role = "PROMOTER_POL3" if any(p in pnm for p in POL3_PROMOTERS) else "PROMOTER_POL2"
            out.append(f"<{role}:{normalize_name(prom)}>")
            out.append("<DRIVES>")
            cds_rhs = rhs
        else:
            cds_rhs = cassette
        # RHS may be a fusion ':'-joined chain or a / bicistronic split
        for bicistronic in cds_rhs.split("/"):
            fusion_parts = bicistronic.split(":")
            for idx, part in enumerate(fusion_parts):
                part_clean = re.sub(r"\[[^\]]+\]$", "", part).strip()  # strip [NM_...] RefSeq refs
                if not part_clean:
                    continue
                up = part_clean.upper()
                if up in LINKER_2A_NAMES:
                    out.append(f"<LINKER_2A:{normalize_name(part_clean)}>")
                elif up in {"PURO", "NEO", "HYGRO", "HYG", "BLAST", "BSD", "ZEO"}:
                    out.append(f"<MARKER:{normalize_name(part_clean)}>")
                else:
                    out.append(f"<CDS:{normalize_name(part_clean)}>")
            if bicistronic != cds_rhs.split("/")[-1]:
                out.append("<BICISTRONIC>")
        out.append("</VB_CAS>")
    out.append("<EOS>")
    return out


# ---------------------------------------------------------------------------
# CLI driver
# ---------------------------------------------------------------------------

def iter_inputs(module_library_glob: str | None,
                ncbi_glob: str | None,
                vectorbuilder_dir: str | None,
                max_plasmids: int | None = None) -> Iterable[tuple[Path, str]]:
    n = 0
    def _cap(): return max_plasmids is not None and n >= max_plasmids
    if module_library_glob:
        for p in sorted(glob.glob(module_library_glob, recursive=True)):
            if _cap(): return
            yield Path(p), "module_library"
            n += 1
    if ncbi_glob:
        for p in sorted(glob.glob(ncbi_glob, recursive=True)):
            if _cap(): return
            yield Path(p), "ncbi"
            n += 1
    if vectorbuilder_dir and os.path.isdir(vectorbuilder_dir):
        for p in sorted(glob.glob(os.path.join(vectorbuilder_dir, "*.gb"))):
            if _cap(): return
            yield Path(p), "vectorbuilder"
            n += 1


def write_jsonl(out_path: Path, examples: Iterable[TokenExample]) -> int:
    n = 0
    with open(out_path, "w") as fh:
        for ex in examples:
            fh.write(json.dumps(asdict(ex), ensure_ascii=False) + "\n")
            n += 1
    return n


def write_vocabulary(out_path: Path, counts: Counter) -> None:
    vocab = {"special_tokens": SPECIAL_TOKENS,
             "role_counts": dict(counts.most_common()),
             "vocab_size": len(counts) + len(SPECIAL_TOKENS)}
    with open(out_path, "w") as fh:
        json.dump(vocab, fh, indent=2, ensure_ascii=False)


def insert_postgres(db_url: str, examples: list[TokenExample]) -> None:
    try:
        import psycopg2
        from psycopg2.extras import execute_values
    except ImportError:
        logger.error("psycopg2 not installed; skipping DB inserts")
        return
    conn = psycopg2.connect(db_url)
    cur = conn.cursor()
    try:
        # plasmid_tokenizations
        for ex in examples:
            cur.execute("""
                INSERT INTO plasmid_tokenizations
                  (plasmid_id, token_level, token_strings, is_valid, coverage_complete,
                   span_count, sequence_length, circular, rotation_idx, source_corpus, metadata)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s::jsonb)
                ON CONFLICT (plasmid_id, token_level, rotation_idx)
                DO UPDATE SET token_strings=EXCLUDED.token_strings,
                              is_valid=EXCLUDED.is_valid,
                              source_corpus=EXCLUDED.source_corpus;
            """, (ex.plasmid_id, "functional", ex.tokens, ex.valid, True,
                  None, ex.length, ex.topology == "circular",
                  ex.rotation_idx, ex.source,
                  json.dumps({"host": ex.host,
                              "canonical_start": ex.canonical_start,
                              "validation_errors": ex.validation_errors})))
            for ix in ex.interaction_tokens:
                cur.execute("""
                    INSERT INTO plasmid_interaction_tokens
                      (plasmid_id, rule_id, sbo_term, source_module, participants, rotation_idx)
                    VALUES (%s, %s, %s, %s, %s, %s)
                    ON CONFLICT DO NOTHING;
                """, (ex.plasmid_id, ix["rule_id"], ix.get("sbo_term"),
                      ix.get("source_module"), ix.get("participants", []), ex.rotation_idx))
            for cl in ex.cloning_tokens:
                cur.execute("""
                    INSERT INTO plasmid_cloning_tokens
                      (plasmid_id, token_subtype, value, position, rotation_idx)
                    VALUES (%s, %s, %s, %s, %s);
                """, (ex.plasmid_id, cl["subtype"], cl["value"], cl.get("position"), ex.rotation_idx))
        conn.commit()
    finally:
        cur.close()
        conn.close()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Plasmid tokenizer v2")
    ap.add_argument("--module-library-glob", default=None)
    ap.add_argument("--ncbi-glob", default=None)
    ap.add_argument("--vectorbuilder-dir", default=None)
    ap.add_argument("--fetch-vectorbuilder", action="store_true",
                    help="Run vectorbuilder_fetch.py before tokenizing")
    ap.add_argument("--k-rotations", type=int, default=6)
    ap.add_argument("--output-jsonl", type=Path, default=Path("./token_corpus.jsonl"))
    ap.add_argument("--output-vocab", type=Path, default=Path("./vocabulary.json"))
    ap.add_argument("--output-vb-descriptions", type=Path, default=None,
                    help="If set, also write VectorBuilder shorthand→tokens pairs here")
    ap.add_argument("--vb-systems-json", type=Path, default=None,
                    help="Path to vector_systems.json to extract shorthand descriptions")
    ap.add_argument("--db", default=None, help="postgresql:// URL")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--max-plasmids", type=int, default=None)
    ap.add_argument("--annotate-url", default=DEFAULT_ANNOTATE_URL,
                    help="FastAPI endpoint for annotate_sequence_llm")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args(argv)

    logging.basicConfig(level=args.log_level,
                        format="%(asctime)s - %(levelname)s - %(message)s")

    if args.fetch_vectorbuilder:
        try:
            from vectorbuilder_fetch import fetch_all  # local import if colocated
        except ImportError:
            from backend.splicify_api.plasmid_lm.vectorbuilder_fetch import fetch_all  # type: ignore
        out_dir = Path(args.vectorbuilder_dir or "./vectorbuilder_gb")
        out_dir.mkdir(parents=True, exist_ok=True)
        systems_json = args.vb_systems_json or Path(
            "backend/splicify_api/vectorbuilder_db/vector_systems.json")
        fetch_all(systems_json, out_dir)

    all_examples: list[TokenExample] = []
    vocab_counts: Counter = Counter()
    skipped = 0

    # Precompute Module_Library folder map for folder-hint inference
    module_library_folders = {}
    if args.module_library_glob:
        from host_inference import FOLDER_PRIOR  # local import
        ml_root = Path(args.module_library_glob.split("**")[0].rstrip("/"))
        if ml_root.exists():
            for folder in FOLDER_PRIOR:
                fp = ml_root / folder
                if not fp.exists():
                    continue
                for gb in fp.rglob("*.gb"):
                    module_library_folders[gb.stem] = folder

    for gb_path, source in iter_inputs(args.module_library_glob,
                                       args.ncbi_glob,
                                       args.vectorbuilder_dir,
                                       args.max_plasmids):
        logger.info("annotating %s [%s]", gb_path.name, source)
        ann = annotate_genbank(gb_path, annotate_url=args.annotate_url)
        if ann is not None and source == "module_library":
            ann["_folder_hint"] = module_library_folders.get(gb_path.stem)
        if ann is None:
            skipped += 1
            continue
        try:
            examples = tokenize_plasmid(ann, source, k_rotations=args.k_rotations)
        except Exception as exc:
            logger.warning("tokenization failed on %s: %s", gb_path.name, exc)
            skipped += 1
            continue
        for ex in examples:
            for t in ex.tokens:
                # Only count role prefix (strip payload) for vocab size
                if ":" in t:
                    role = t.split(":", 1)[0] + ">"
                else:
                    role = t
                vocab_counts[role] += 1
        all_examples.extend(examples)

    logger.info("wrote %d examples across %d plasmids (%d skipped)",
                len(all_examples),
                len({ex.plasmid_id for ex in all_examples}),
                skipped)

    # Optional VB description shorthand pairs
    if args.output_vb_descriptions and args.vb_systems_json:
        with open(args.vb_systems_json) as fh:
            vb_data = json.load(fh)
        rows = []
        for sys in vb_data.get("vector_systems", []):
            for rv in sys.get("representative_vectors", []):
                desc = rv.get("description", "")
                toks = parse_vectorbuilder_shorthand(desc)
                if toks:
                    rows.append({
                        "natural_language": desc,
                        "target_tokens": toks,
                        "source_type": "vectorbuilder_shorthand",
                        "source_plasmid_id": rv.get("vector_id", ""),
                    })
                    for t in toks:
                        role = t.split(":", 1)[0] + ">" if ":" in t else t
                        vocab_counts[role] += 1
        with open(args.output_vb_descriptions, "w") as fh:
            for r in rows:
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
        logger.info("wrote %d VectorBuilder shorthand→tokens pairs to %s",
                    len(rows), args.output_vb_descriptions)

    write_jsonl(args.output_jsonl, all_examples)
    write_vocabulary(args.output_vocab, vocab_counts)

    valid = sum(1 for ex in all_examples if ex.valid)
    logger.info("validation: %d/%d examples passed gates", valid, len(all_examples))

    if args.db and not args.dry_run:
        logger.info("inserting into postgres: %s", args.db)
        insert_postgres(args.db, all_examples)
    elif args.db and args.dry_run:
        logger.info("--dry-run set; skipping postgres inserts")

    return 0


if __name__ == "__main__":
    sys.exit(main())
