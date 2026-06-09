"""
describe_plasmid handler — semantic-search-and-edit replacement for the old
LLM-driven plasmid_design path.

Flow (no LLM in v1):
  1. Build a PlasmidSpec from the prompt + KB-resolved features + any
     uploaded modules.
  2. Embed the spec as a sentence (same model used to build the corpus
     index) and query the HNSW retrieval index for the top-k plasmid hits.
  3. Resolve the top-1 hit's source sequence (from Module_Library_gb or the
     NCBI raw NDJSON) and run it through the full annotation pipeline →
     foundation.gb (LLM-annotated, replacing the original annotation).
  4. Diff foundation.modules vs spec.modules_required → edit_ops.
     Each missing module is classified as a primer-tail edit (<40 bp insert)
     or a synthesis fragment (>=40 bp). KB sequence is used when available.
  5. Apply edit_ops to the foundation sequence (deterministic; for v1 we
     append synthesis fragments and emit primers for tail edits — the
     orchestrator-driven splice-by-position is a Phase 4 follow-up).
  6. Re-annotate the edited sequence at modules_only depth → designed.gb.
  7. Return: foundation.gb, designed.gb, viz, workflow_trace.json.
"""
from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .workflow_input import WorkflowInput

import io
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from . import _data
from typing import Any, Dict, List, Optional, Tuple

from .annotation_cache import annotate_cached
from .plasmid_spec import PlasmidSpec, SpecModule, diff_spec_against_target

logger = logging.getLogger("describe_plasmid")

# Cached singletons --------------------------------------------------------
_RETRIEVAL_DIR = Path("/var/data/plasmid_lm_corpus/retrieval")
_TOKEN_CORPUS  = Path("/var/data/plasmid_lm_corpus/token_corpus.jsonl")
_MODULE_LIB    = _data.data_path("Module_Library_gb")
_NCBI_RAW      = Path("/var/data/plasmid_lm_corpus/raw_plasmid_records_refresh.ndjson")

_PRIMER_TAIL_THRESHOLD = 40  # bp — at or above this we suggest synthesis

_RETRIEVAL: Optional[Dict[str, Any]] = None  # {ids, embeddings|None, sentences, mode}
_PLASMID_SOURCES: Optional[Dict[str, Dict[str, Any]]] = None  # plasmid_id -> {tokens, source}


def _tokens_to_sentence(tokens: List[str]) -> str:
    """Flatten role tokens to a space-separated sentence — matches the format
    used by build_retrieval_index.py so embedding-based retrieval is comparable
    and Jaccard fallback uses the same vocabulary."""
    NOISE = ("<BOS>", "<EOS>", "<PAD>", "<UNK>",
            "<TOPOLOGY:", "<LEN_BIN:", "<ROTATION_IDX:",
            "<MOD_CLOSE>", "<UPSTREAM_REG>", "</UPSTREAM_REG>",
            "<CDS_MOD>", "</CDS_MOD>",
            "<DOWNSTREAM_REG>", "</DOWNSTREAM_REG>",
            "<VB_CAS>", "</VB_CAS>", "<DRIVES>")
    out: List[str] = []
    for t in tokens:
        if any(t.startswith(p) for p in NOISE):
            continue
        out.append(t.strip("<>").replace(":", " ").replace("_", " "))
    return " ".join(out)


# ---------------------------------------------------------------------------
# Edit operations
# ---------------------------------------------------------------------------
@dataclass
class EditOp:
    op:       str  # "add" | "remove" | "replace"
    target:   str  # name/role being added/removed
    role:     Optional[str] = None
    sequence: Optional[str] = None
    length:   Optional[int] = None
    strategy: str = "synthesis"  # "primer_tail" | "synthesis" | "none"
    rationale: str = ""

    def to_dict(self) -> Dict[str, Any]:
        d = dict(self.__dict__)
        if d["sequence"] and len(d["sequence"]) > 80:
            d["sequence_preview"] = d["sequence"][:60] + "…"
        return d


# ---------------------------------------------------------------------------
# Retrieval (HNSW over MiniLM embeddings — same model used at build time)
# ---------------------------------------------------------------------------
def _load_retrieval() -> Dict[str, Any]:
    """Load retrieval data. Prefers the precomputed MiniLM embeddings + a
    sentence-transformer query encoder when both numpy and sentence-transformers
    are available; falls back to Jaccard over flattened role-token sentences
    otherwise. No dependence on hnswlib in either path — brute-force cosine on
    7k × 384 is microseconds.
    """
    global _RETRIEVAL
    if _RETRIEVAL is not None:
        return _RETRIEVAL

    meta = {}
    meta_path = _RETRIEVAL_DIR / "retrieval_meta.json"
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text())
        except Exception:
            meta = {}

    ids_path = _RETRIEVAL_DIR / "plasmid_ids.txt"
    if not ids_path.exists():
        raise FileNotFoundError(f"missing {ids_path}")
    ids = ids_path.read_text().splitlines()

    # Build the per-plasmid sentence used for both Jaccard fallback and
    # (optionally) for re-embedding queries on the fly.
    sources = _load_plasmid_sources()
    sentences = [
        _tokens_to_sentence(sources.get(pid, {}).get("tokens") or [])
        for pid in ids
    ]

    state: Dict[str, Any] = {"ids": ids, "sentences": sentences, "meta": meta,
                             "embeddings": None, "encoder": None, "mode": "jaccard"}

    # Try the embedding path
    try:
        import numpy as np
        emb_path = _RETRIEVAL_DIR / "plasmid_embeddings.npy"
        if emb_path.exists():
            state["embeddings"] = np.load(str(emb_path))
        try:
            from sentence_transformers import SentenceTransformer  # type: ignore
            state["encoder"] = SentenceTransformer(meta.get("model") or
                                                    "sentence-transformers/all-MiniLM-L6-v2")
            state["mode"] = "embedding"
        except Exception as exc:
            logger.info("sentence-transformers unavailable (%s) — falling back to Jaccard", exc)
    except Exception as exc:
        logger.info("numpy unavailable (%s) — falling back to Jaccard", exc)

    _RETRIEVAL = state
    logger.info("retrieval loaded: %d plasmids, mode=%s", len(ids), state["mode"])
    return state


def _load_plasmid_sources() -> Dict[str, Dict[str, Any]]:
    """Map plasmid_id → {tokens, source} so we know how to resolve the
    foundation sequence (Module_Library_gb file vs NCBI sequence)."""
    global _PLASMID_SOURCES
    if _PLASMID_SOURCES is not None:
        return _PLASMID_SOURCES
    out: Dict[str, Dict[str, Any]] = {}
    with open(_TOKEN_CORPUS) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("rotation_idx", 0) != 0:
                continue
            pid = r["plasmid_id"]
            if pid in out:
                continue
            out[pid] = {"source": r.get("source"), "tokens": r.get("tokens", [])}
    _PLASMID_SOURCES = out
    return out


def _spec_to_query_sentence(spec: PlasmidSpec) -> str:
    """Build a retrieval query string from the spec, matching the format used
    when the index was built (tokens flattened to space-separated words)."""
    parts: List[str] = [spec.description]
    if spec.host_hint:
        parts.append(f"host {spec.host_hint}")
    parts.append(f"topology {spec.topology}")
    for m in spec.modules_required + spec.modules_present:
        if m.role:
            parts.append(f"{m.role} {m.name}")
        else:
            parts.append(m.name)
    return " ".join(p for p in parts if p)


def _retrieve_top_k(spec: PlasmidSpec, k: int = 5) -> List[Dict[str, Any]]:
    r = _load_retrieval()
    q_text = _spec_to_query_sentence(spec)

    if r["mode"] == "embedding" and r["embeddings"] is not None and r["encoder"] is not None:
        import numpy as np
        q = r["encoder"].encode([q_text], normalize_embeddings=True)
        sims = (r["embeddings"] @ np.asarray(q, dtype=np.float32).T).ravel()
        order = np.argsort(-sims)[:k]
        return [
            {"plasmid_id": r["ids"][i], "score": float(sims[i]), "rank": rk + 1}
            for rk, i in enumerate(order)
        ]

    # Jaccard fallback: bag-of-tokens over the flattened sentence.
    q_tokens = set(q_text.lower().split())
    if not q_tokens:
        return []
    scored: List[Tuple[float, int]] = []
    for i, sent in enumerate(r["sentences"]):
        s = set(sent.lower().split())
        if not s:
            continue
        inter = len(q_tokens & s)
        if inter == 0:
            continue
        union = len(q_tokens | s)
        scored.append((inter / union, i))
    scored.sort(reverse=True)
    return [
        {"plasmid_id": r["ids"][i], "score": float(score), "rank": rk + 1}
        for rk, (score, i) in enumerate(scored[:k])
    ]


# ---------------------------------------------------------------------------
# Resolve foundation sequence
# ---------------------------------------------------------------------------
def _resolve_foundation_sequence(plasmid_id: str) -> Optional[Tuple[str, str]]:
    """Return (sequence, source_label) for a plasmid_id, or None."""
    sources = _load_plasmid_sources()
    info = sources.get(plasmid_id) or {}
    source = info.get("source")

    if source == "module_library":
        # Match basename (case sensitive) under Module_Library_gb/*/
        for gb in _MODULE_LIB.rglob(f"{plasmid_id}.gb"):
            try:
                from Bio import SeqIO
                rec = SeqIO.read(gb, "genbank")
                return str(rec.seq).upper(), f"module_library:{gb.relative_to(_MODULE_LIB)}"
            except Exception as exc:
                logger.warning("could not read %s: %s", gb, exc)
                return None
    if source in ("ncbi", "ncbi_engineered") and _NCBI_RAW.exists():
        # NCBI rows are stored without plasmid_id; we can't directly match.
        # In this v1 we only resolve module_library hits; NCBI hits return
        # None and the caller falls back to the next ranked hit.
        return None
    return None


# ---------------------------------------------------------------------------
# Diff → edit ops
# ---------------------------------------------------------------------------
def _classify_edit_strategy(seq: Optional[str]) -> str:
    if not seq:
        return "none"
    if len(seq) < _PRIMER_TAIL_THRESHOLD:
        return "primer_tail"
    return "synthesis"


def _build_edit_ops(spec: PlasmidSpec, target_modules: List[Dict[str, Any]]) -> List[EditOp]:
    diff = diff_spec_against_target(spec, target_modules)
    spec_by_name = {(m.name or "").lower(): m for m in spec.modules_required}

    ops: List[EditOp] = []
    for entry in diff["missing"]:
        nm = (entry.get("name") or "").lower()
        spec_mod = spec_by_name.get(nm)
        seq = spec_mod.sequence if spec_mod else None
        strategy = _classify_edit_strategy(seq)
        rationale = (
            f"spec requires {entry.get('role') or 'module'} '{entry.get('name')}' "
            f"but it is absent from the foundation. "
        )
        if strategy == "primer_tail":
            rationale += f"Insert is {len(seq)} bp — embed in primer tail."
        elif strategy == "synthesis":
            rationale += f"Insert is {len(seq)} bp — order as gBlock / synthesis fragment."
        else:
            rationale += "Sequence is not in the KB — manual sourcing required."
        ops.append(EditOp(
            op="add",
            target=entry.get("name") or "",
            role=entry.get("role"),
            sequence=seq,
            length=len(seq) if seq else None,
            strategy=strategy,
            rationale=rationale,
        ))
    for entry in diff["unexpected"]:
        ops.append(EditOp(
            op="remove",
            target=str(entry.get("name") or ""),
            role=entry.get("role"),
            strategy="primer_tail",
            rationale="present on foundation but not in user spec — drop or flag for review.",
        ))
    return ops


def _apply_edit_ops(
    foundation_seq: str,
    foundation_topology: str,
    ops: List[EditOp],
) -> Tuple[str, str, List[Dict[str, Any]]]:
    """Apply edit_ops deterministically.

    v1 behaviour: append synthesis-class additions to the foundation sequence
    (preserving topology). Primer-tail and remove ops are recorded but not
    applied to the sequence; their handling lives in the cloning workflow
    primer designer downstream. The orchestrator (Phase 4) will splice by
    position once it knows where each module belongs.
    """
    edited = foundation_seq
    notes: List[Dict[str, Any]] = []
    for op in ops:
        if op.op == "add" and op.strategy == "synthesis" and op.sequence:
            edited += op.sequence.upper()
            notes.append({
                "op": "add",
                "target": op.target,
                "applied": True,
                "method": "appended synthesis fragment",
                "delta_bp": len(op.sequence),
            })
        else:
            notes.append({
                "op": op.op,
                "target": op.target,
                "applied": False,
                "method": op.strategy,
                "deferred_to": "primer designer / orchestrator",
            })
    return edited, foundation_topology, notes


# ---------------------------------------------------------------------------
# .gb + viz + trace builders
# ---------------------------------------------------------------------------
def _annotations_from_modules(mods: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for m in mods or []:
        try:
            start = int(m.get("start", 0))
            end = int(m.get("end", start + 1))
        except (TypeError, ValueError):
            continue
        out.append({
            "name": m.get("name") or m.get("module_type") or "feature",
            "start": start,
            "end": end,
            "direction": -1 if m.get("strand") == -1 else 1,
            "color": m.get("color") or "#6B7280",
            "feat_type": m.get("type") or m.get("module_type") or "misc_feature",
        })
    return out


def _build_response_files(
    foundation_seq: str, foundation_modules: List[Dict[str, Any]],
    designed_seq: str, designed_modules: List[Dict[str, Any]],
    foundation_topology: str, designed_topology: str,
    plasmid_id: str,
    trace: Dict[str, Any],
) -> List[Dict[str, str]]:
    from .files import _make_genbank, _file
    files: List[Dict[str, str]] = []
    files.append(_file(
        f"{plasmid_id}_foundation.gb",
        "application/octet-stream",
        _make_genbank(
            seq=foundation_seq,
            name=plasmid_id[:16] or "foundation",
            description=f"Top-hit foundation plasmid (LLM annotation)",
            annotations=_annotations_from_modules(foundation_modules),
            topology=foundation_topology,
        ),
    ))
    files.append(_file(
        "designed_plasmid.gb",
        "application/octet-stream",
        _make_genbank(
            seq=designed_seq,
            name="designed",
            description="Edited foundation matching the user's described plasmid",
            annotations=_annotations_from_modules(designed_modules),
            topology=designed_topology,
        ),
    ))
    files.append(_file(
        "workflow_trace.json",
        "application/json",
        json.dumps(trace, indent=2, default=str).encode("utf-8"),
    ))
    return files


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
async def describe_plasmid(
    workflow_input: "WorkflowInput",
) -> Dict[str, Any]:
    """Step 3f of the workflow_input migration: this handler now consumes a
    WorkflowInput. The PlasmidSpec is pre-built by
    `build_for_plasmid_design` and stashed on
    `workflow_input.provenance["plasmid_spec_object"]`."""
    from .plasmid_spec import build_plasmid_spec
    session_id = workflow_input.session_id
    spec = workflow_input.provenance.get("plasmid_spec_object")
    if spec is None:
        # Defensive fallback: construct the spec from message + intent_result
        # if the adapter wasn't used. Path normally unreachable but keeps the
        # handler self-sufficient for direct callers.
        message = workflow_input.provenance.get("message", "")
        intent_result = workflow_input.provenance.get("intent_result") or {}
        spec = build_plasmid_spec(message=message, intent_result=intent_result)

    # 1. Retrieve top-k plasmids from the index
    try:
        hits = _retrieve_top_k(spec, k=5)
    except Exception as exc:
        logger.error("retrieval failed: %s", exc, exc_info=True)
        return {
            "reply": (
                "## Plasmid description\n\n"
                f"Retrieval index unavailable: `{exc}`. "
                "Make sure `/var/data/plasmid_lm_corpus/retrieval/` is populated."
            ),
            "viz": None,
            "files": [],
        }

    # 2. Resolve foundation sequence — walk the ranked hits until one resolves
    foundation_seq = None
    foundation_topology = "circular"
    chosen_hit: Optional[Dict[str, Any]] = None
    for h in hits:
        seq_pair = _resolve_foundation_sequence(h["plasmid_id"])
        if seq_pair:
            foundation_seq, source_label = seq_pair
            chosen_hit = {**h, "source": source_label}
            break
    if not foundation_seq or not chosen_hit:
        return {
            "reply": (
                "## Plasmid description\n\n"
                "Top retrieval hits could not be resolved to a sequence "
                "(NCBI-only hits not yet supported in v1). "
                f"Top-{len(hits)} ids: " + ", ".join(h['plasmid_id'] for h in hits)
            ),
            "viz": None,
            "files": [],
        }

    # 3. Annotate foundation at full depth (replacing any original annotation).
    # v2 (2026-04-28): the handler stops here. We surface the foundation as
    # the answer — annotated by our pipeline with features, modules, and
    # interactions — and describe it back to the user. No edit_ops, no
    # designed plasmid, no orchestrator. The user picks whether the
    # foundation is a good starting point.
    foundation_ann = await annotate_cached(
        foundation_seq, circular=True, depth="full",
    )
    foundation_modules = (
        foundation_ann.get("module_annotations") or foundation_ann.get("modules") or []
    )
    foundation_features = foundation_ann.get("annotations") or []
    foundation_interactions = foundation_ann.get("interactions") or []

    trace = {
        "session_id": session_id,
        "spec": spec.to_dict(),
        "retrieval": hits,
        "foundation": {
            "plasmid_id": chosen_hit["plasmid_id"],
            "source": chosen_hit["source"],
            "score": chosen_hit["score"],
            "length": len(foundation_seq),
            "topology": foundation_topology,
            "module_count": len(foundation_modules),
            "feature_count": len(foundation_features),
            "interaction_count": len(foundation_interactions),
        },
        "modules": foundation_modules,
        "interactions": foundation_interactions,
    }

    files = _build_foundation_files(
        foundation_seq, foundation_modules, foundation_topology,
        chosen_hit["plasmid_id"], trace,
    )

    viz = {
        "type": "design",
        "title": chosen_hit["plasmid_id"],
        "sequence": foundation_seq,
        "topology": foundation_topology,
        "total_length": len(foundation_seq),
        "annotations": _annotations_from_modules(foundation_modules),
        "method": "describe_plasmid",
    }

    reply = _format_foundation_reply(
        spec, chosen_hit, foundation_modules, foundation_interactions,
        foundation_features,
    )

    return {"reply": reply, "viz": viz, "files": files, "trace": trace}


def _build_foundation_files(
    foundation_seq: str,
    foundation_modules: List[Dict[str, Any]],
    foundation_topology: str,
    plasmid_id: str,
    trace: Dict[str, Any],
) -> List[Dict[str, str]]:
    """Foundation-only files: annotated .gb of the top-hit plasmid + trace.

    No designed_plasmid.gb because v2 does not edit the foundation.
    """
    from .files import _make_genbank, _file
    files: List[Dict[str, str]] = []
    files.append(_file(
        f"{plasmid_id}_foundation.gb",
        "application/octet-stream",
        _make_genbank(
            seq=foundation_seq,
            name=plasmid_id[:16] or "foundation",
            description="Foundation plasmid (clean-room annotation pipeline)",
            annotations=_annotations_from_modules(foundation_modules),
            topology=foundation_topology,
        ),
    ))
    files.append(_file(
        "workflow_trace.json",
        "application/json",
        json.dumps(trace, indent=2, default=str).encode("utf-8"),
    ))
    return files


def _format_foundation_reply(
    spec: PlasmidSpec,
    hit: Dict[str, Any],
    modules: List[Dict[str, Any]],
    interactions: List[Dict[str, Any]],
    features: List[Dict[str, Any]],
) -> str:
    """Markdown reply describing what's on the foundation plasmid.

    Groups annotations by role/module type so the user can quickly see what's
    on the recommended plasmid. No edit suggestions; no design output.
    """
    L: List[str] = []
    L.append("## Recommended foundation plasmid")
    L.append("")
    L.append(
        f"**`{hit['plasmid_id']}`** — retrieval score {hit['score']:.3f}, "
        f"source `{hit['source']}`"
    )
    L.append("")
    L.append(
        f"This is the closest match in the LLM-annotated plasmid corpus "
        f"({len(modules)} module(s), {len(features)} feature(s), "
        f"{len(interactions)} interaction(s)). Use it as a starting point — "
        "we'll iterate on edits in a follow-up turn if you want to modify it."
    )
    L.append("")

    # Modules grouped by type
    if modules:
        L.append("### Modules")
        by_type: Dict[str, List[str]] = {}
        for m in modules:
            mt = str(m.get("module_type") or m.get("type") or "module")
            label = str(m.get("name") or m.get("label") or m.get("module_type") or "—")
            start = m.get("start")
            end = m.get("end")
            loc = f" ({start}–{end})" if start is not None and end is not None else ""
            by_type.setdefault(mt, []).append(f"{label}{loc}")
        for mt in sorted(by_type):
            entries = by_type[mt]
            L.append(f"- **{mt}** ({len(entries)}): " + ", ".join(entries[:6])
                     + ("…" if len(entries) > 6 else ""))
        L.append("")

    # Features grouped by feat_type / role
    if features:
        L.append("### Features")
        by_role: Dict[str, List[str]] = {}
        for f in features:
            role = str(f.get("feat_type") or f.get("type") or f.get("role") or "feature")
            name = str(f.get("name") or f.get("feature_name") or "—")
            start = f.get("start")
            end = f.get("end")
            loc = f" ({start}–{end})" if start is not None and end is not None else ""
            by_role.setdefault(role, []).append(f"{name}{loc}")
        for role in sorted(by_role):
            entries = by_role[role]
            L.append(f"- **{role}** ({len(entries)}): " + ", ".join(entries[:8])
                     + ("…" if len(entries) > 8 else ""))
        L.append("")

    # Interactions
    if interactions:
        L.append("### Interactions")
        seen: Dict[str, int] = {}
        for ix in interactions:
            rid = str(ix.get("rule_id") or ix.get("rule") or "interaction")
            seen[rid] = seen.get(rid, 0) + 1
        for rid in sorted(seen):
            L.append(f"- `{rid}` × {seen[rid]}")
        L.append("")

    # What the user asked for, for context
    if spec.modules_required:
        L.append("### You asked for")
        for m in spec.modules_required:
            L.append(f"- {m.name} `{m.role or '?'}`")
        L.append("")

    L.append("Files: `<plasmid_id>_foundation.gb`, `workflow_trace.json`.")
    return "\n".join(L)
