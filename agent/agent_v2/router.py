"""POST /agent_v2/chat (multipart or JSON), GET /agent_v2/health."""
import os
from fastapi import APIRouter
import asyncio
import json
from typing import Any, Optional

from fastapi import Request
from fastapi.responses import StreamingResponse

from agent_v2.orchestrator import OrchestratorDeps, run_orchestrator

router = APIRouter()


@router.get("/health")
async def health() -> dict:
    return {
        "ok": True,
        "service": "agent_v2_api",
        "model": os.getenv("AGENT_MODEL", "claude-sonnet-4-6"),
        "anthropic_sdk_installed": _anthropic_installed(),
        "api_key_set": bool(os.getenv("ANTHROPIC_API_KEY")),
        "redis_url_set": bool(os.getenv("REDIS_URL")),
    }


@router.post("/chat")
async def chat(request: Request) -> dict:
    """Non-streaming chat endpoint. Multipart or JSON body; JSON envelope response.

    Multipart fields: message, session_id?, file (target .gb), inventory_files[].
    JSON fields:      {message, session_id?, target_genbank?, inventory_genbank?}.
    """
    message, session_id, registry = await _parse_chat_body(request)
    try:
        envelope = await run_orchestrator(message, registry, session_id=session_id)
    except Exception as e:
        return {
            "ok": False,
            "reply": "",
            "files": None,
            "viz": None,
            "agent_trace": [],
            "n_tool_calls": 0,
            "error": f"{type(e).__name__}: {e}",
            "session_id": session_id,
        }
    return envelope


def _anthropic_installed() -> bool:
    try:
        import anthropic  # noqa: F401
        return True
    except ImportError:
        return False

from fastapi import File, UploadFile


@router.post("/annotate-on-upload")
async def annotate_on_upload(file: UploadFile = File(...)) -> dict:
    """Eagerly annotate an uploaded .gb file so the frontend can paint the
    viewer immediately, before the chat call.

    Routes by kind:
      - plasmid (synthetic construct, circular)  -> v1 annotate_llm_cached
        (KB-driven feature_annotator + Pol2 + interactions + cloning_features)
      - genomic (natural organism, linear, NCBI/RefSeq) -> v2 genomic_annotator
        (native GenBank gene/mRNA/CDS/exon/misc_feature, no KB tampering)

    Without this branching, KEAP1.gb-style records get hit with the plasmid
    KB scanner and surface false-positive features (tetracysteine tag, hGH
    polyA signal, etc.) — wrong by construction for a chromosomal slice.

    No agent involvement, no LLM cost. Returns viz + the FileKind digest so
    the frontend can adjust its UI per kind.
    """
    raw = (await file.read()).decode("utf-8", errors="replace")
    if not raw.strip():
        return {"ok": False, "error": "empty file"}

    import agent_v2  # noqa: F401 — triggers path shim
    from splicify_api.agent.agent_tools import extract_seq_from_genbank
    from agent_v2.file_kind import classify_genbank

    seq = extract_seq_from_genbank(raw)
    if not seq or len(seq) < 50:
        return {"ok": False, "error": "no plasmid sequence found in upload"}

    title = (file.filename or "uploaded").rsplit(".", 1)[0] or "uploaded"
    classification = classify_genbank(raw)
    is_circular = (classification.topology == "circular") or (classification.kind == "plasmid")

    if classification.kind == "genomic":
        try:
            from agent_v2.genomic_annotator import annotate_genomic_gb
            ann = annotate_genomic_gb(raw)
        except Exception as e:
            return {"ok": False, "error": f"genomic annotation failed: {type(e).__name__}: {e}"}

        # Shape genomic features into SeqViz-compatible annotation rows so the
        # viewer can render them without learning a new schema. type maps:
        # gene/mRNA -> the same names; CDS -> CDS (with /translation in
        # description); exon/misc_feature pass through. We expose only fields
        # the viewer reads (name, type, start, end, direction).
        annotations = []
        for ft in ann.features:
            # Detect features the CRISPR / cloning pipeline injected on a
            # prior run (qualifier /added_by="crispr_v2" written by
            # emit_guides_gb). Surfaced as added_by_design so the viewer
            # colors them distinctly from native genomic annotations.
            ft_q = ft.qualifiers or {}
            added_by_val = ft_q.get("added_by")
            if isinstance(added_by_val, list):
                added_by_val = added_by_val[0] if added_by_val else None
            added_by_design = bool(added_by_val and (
                "crispr" in str(added_by_val).lower()
                or "design" in str(added_by_val).lower()))
            for (seg_start, seg_end) in (ft.intervals or [(ft.start, ft.end)]):
                annotations.append({
                    "name": ft.label or ft.gene or ft.transcript_id or ft.protein_id or ft.type,
                    "type": ft.type,
                    "start": int(seg_start),
                    "end": int(seg_end),
                    "direction": ft.strand,
                    "strand": ft.strand,
                    "source": "genbank_native",
                    "description": _genomic_feature_description(ft),
                    "added_by_design": added_by_design or None,
                })

        # Emit per-exon translation annotations so the viewer's AA-strip
        # (layer="translation") renders codon-by-codon across each CDS exon
        # segment — same shape the plasmid pipeline emits from orf_detection
        # but with one entry per exon to handle complement(join(...)) coords.
        for ft in ann.features:
            if ft.type != "CDS" or not ft.translation:
                continue
            annotations.extend(_genomic_translation_annotations(ft))

        viz = {
            "type": "genomic",
            "title": title,
            "sequence": seq,
            "circular": False,   # genomic is always linear
            "annotations": annotations,
            "modules": [],
            "cloning_features": [],
            "interactions": [],
            "hierarchical_annotations": [],
            "genomic_summary": {
                "organism": ann.organism,
                "accession": ann.accession,
                "chromosome": ann.chromosome,
                "length_bp": ann.length_bp,
                "n_genes": len(ann.genes),
                "n_transcripts": len(ann.transcripts),
                "n_features": len(ann.features),
                "genes": list(ann.genes.keys()),
                "transcripts": list(ann.transcripts.keys()),
            },
        }
        return {
            "ok": True, "viz": viz, "length_bp": len(seq), "title": title,
            "kind": "genomic",
            "kind_confidence": classification.confidence,
            "kind_signals": classification.signals,
        }

    # Plasmid kind (or unknown — assume plasmid for the safe default).
    try:
        from splicify_api.annotation_cache import annotate_llm_cached
        ann = await annotate_llm_cached(seq, circular=is_circular)
    except Exception as e:
        return {"ok": False, "error": f"annotation failed: {type(e).__name__}: {e}"}

    viz = {
        "type": "plasmid",
        "title": title,
        "sequence": seq,
        "circular": is_circular,
        "annotations": ann.get("annotations") or ann.get("features") or [],
        "modules": ann.get("modules") or [],
        "cloning_features": ann.get("cloning_features") or [],
        "interactions": ann.get("interactions") or [],
        "hierarchical_annotations": ann.get("hierarchical_annotations") or [],
    }
    return {
        "ok": True, "viz": viz, "length_bp": len(seq), "title": title,
        "kind": classification.kind,
        "kind_confidence": classification.confidence,
        "kind_signals": classification.signals,
    }


def _genomic_translation_annotations(ft) -> list[dict[str, object]]:
    """Walk a CDS's exon intervals in spliced reading order, emit one
    translation annotation per exon carrying that exon's slice of the
    spliced AA sequence in metadata.aa_sequence.

    For + strand: exons read left-to-right on the genome.
    For - strand: exons read right-to-left on the genome (rightmost first).
    cumulative_bp tracks position within the SPLICED CDS so codon
    boundaries map cleanly to AA indices, even when exons split a codon.
    """
    intervals = sorted(ft.intervals or [], key=lambda p: p[0])
    if ft.strand < 0:
        intervals = list(reversed(intervals))
    if not intervals:
        return []
    total_aa = len(ft.translation)
    n_exons = len(intervals)
    label_name = ft.gene or ft.transcript_id or ft.protein_id or "CDS"

    out: list[dict[str, object]] = []
    cumulative_bp = 0
    for exon_idx, (seg_start, seg_end) in enumerate(intervals, start=1):
        exon_len = seg_end - seg_start
        if exon_len <= 0:
            continue
        aa_start = (cumulative_bp // 3) + 1
        # Ceiling division for end so partial-codon exons still cover the AA they contribute to.
        aa_end = min((cumulative_bp + exon_len + 2) // 3, total_aa)
        aa_end = max(aa_end, aa_start)
        aa_slice = ft.translation[aa_start - 1: aa_end]
        suffix = (f" exon {exon_idx}/{n_exons}" if n_exons > 1 else "")
        out.append({
            "name": f"Translation {label_name}{suffix} ({aa_start}-{aa_end} aa)",
            "type": "translation",
            "start": int(seg_start),
            "end": int(seg_end),
            "direction": ft.strand,
            "strand": ft.strand,
            "color": "#673AB7",
            "layer": "translation",
            "module_type": "translation",
            "source": "orf_detection",
            "metadata": {
                "aa_length": len(aa_slice),
                "aa_sequence": aa_slice,
                "aa_start_global": aa_start,
                "aa_end_global": aa_end,
                "cds_total_aa": total_aa,
                "exon_idx": exon_idx,
                "n_exons": n_exons,
                "gene": ft.gene,
                "transcript_id": ft.transcript_id,
                "protein_id": ft.protein_id,
                "feature_regions": [],
                "orf_detected": True,
            },
        })
        cumulative_bp += exon_len
    return out


def _genomic_feature_description(ft) -> str:
    """Build a one-line description for a genomic feature so the viewer's
    annotation card shows useful metadata."""
    parts = []
    if ft.gene:
        parts.append(f"gene={ft.gene}")
    if ft.transcript_id:
        parts.append(f"transcript={ft.transcript_id}")
    if ft.protein_id:
        parts.append(f"protein={ft.protein_id}")
    if ft.type == "CDS" and ft.translation:
        parts.append(f"len={len(ft.translation)} aa")
    elif ft.intervals and len(ft.intervals) > 1:
        parts.append(f"n_exons={len(ft.intervals)}")
    if ft.upgraded_from:
        parts.append(f"upgraded_from={ft.upgraded_from}")
    if ft.note:
        n = ft.note.strip().replace("\n", " ")
        if len(n) > 120:
            n = n[:117] + "..."
        parts.append(n)
    return " | ".join(parts)

async def _add_genbank_to_registry(registry: Any, raw: str, fallback_name: str) -> bool:
    """Extract a sequence from GenBank text, classify plasmid vs genomic,
    register with the right topology, and stash the classification + raw
    GenBank text for downstream genomic-aware tools.

    For plasmid uploads, eagerly annotate via annotate_llm_cached and
    seed the per-attachment annotation cache so find_features /
    graft_parts can read feature coordinates without re-annotating.
    Annotation failures are non-fatal (lazy retry happens on first
    tool use).
    """
    import agent_v2  # noqa: F401 — triggers path shim
    from splicify_api.agent.agent_tools import extract_seq_from_genbank
    from agent_v2 import attachment_kinds
    from agent_v2.file_kind import classify_genbank

    seq = extract_seq_from_genbank(raw or "")
    if not seq or len(seq) < 50:
        return False
    name = fallback_name.rsplit(".", 1)[0] if "." in fallback_name else fallback_name

    classification = classify_genbank(raw or "")
    is_circular = (classification.topology == "circular") or (classification.kind == "plasmid")
    aid = registry.register_product(name or "attachment", seq, circular=is_circular)
    attachment_kinds.stash_kind(aid, classification, gb_text=raw)

    # Eager annotation for plasmid uploads. Genomic records go through
    # a different pipeline (NCBI annotator) and skip this branch.
    if classification.kind != "genomic":
        try:
            from splicify_api.annotation_cache import annotate_llm_cached
            from agent_v2.tools import cache_annotation
            envelope = await annotate_llm_cached(seq, circular=is_circular)
            if isinstance(envelope, dict):
                envelope["sequence"] = seq
                cache_annotation(aid, envelope)
        except Exception as e:
            import logging as _l
            _l.getLogger(__name__).warning(
                "upload annotation failed for %s (%s); tools will retry lazily",
                aid, f"{type(e).__name__}: {e}",
            )
    return True


async def _build_registry_from_json(body: dict) -> Any:
    """Build an AttachmentRegistry from JSON fields target_genbank /
    inventory_genbank. Async because per-attachment annotation is
    eager (see _add_genbank_to_registry)."""
    import agent_v2  # noqa: F401 — triggers path shim
    from splicify_api.agent.agent_tools import AttachmentRegistry

    reg = AttachmentRegistry()
    target = body.get("target_genbank")
    if isinstance(target, str) and target.strip():
        await _add_genbank_to_registry(reg, target, "target")
    inventory = body.get("inventory_genbank") or []
    if isinstance(inventory, str):
        inventory = [inventory]
    for i, gb in enumerate(inventory):
        if isinstance(gb, str) and gb.strip():
            await _add_genbank_to_registry(reg, gb, f"inventory_{i}")
    return reg


async def _parse_chat_body(request: Request) -> tuple[str, Optional[str], Any]:
    """Unified body parser. Returns (message, session_id, registry).

    Supports multipart/form-data with `message`, `session_id`, `file` (target),
    `inventory_files[]`, and application/json with the same fields plus
    `target_genbank` / `inventory_genbank`.
    """
    import agent_v2  # noqa: F401 — triggers path shim
    from splicify_api.agent.agent_tools import AttachmentRegistry

    content_type = request.headers.get("content-type", "")
    if content_type.startswith("multipart/form-data"):
        form = await request.form()
        message = str(form.get("message") or "")
        session_id = form.get("session_id")
        if session_id is not None:
            session_id = str(session_id)
        registry = AttachmentRegistry()

        target = form.get("file")
        if target is not None and hasattr(target, "read"):
            raw = (await target.read()).decode("utf-8", errors="replace")
            await _add_genbank_to_registry(
                registry, raw, target.filename or "target",
            )

        for f in form.getlist("inventory_files"):
            if hasattr(f, "read"):
                raw = (await f.read()).decode("utf-8", errors="replace")
                await _add_genbank_to_registry(
                    registry, raw,
                    f.filename or f"inventory_{len(registry.public_summary())}",
                )
        return message, session_id, registry

    # JSON fallback (also handles unknown / no content-type gracefully)
    try:
        body = await request.json()
    except Exception:
        body = {}
    body = body or {}
    return str(body.get("message") or ""), body.get("session_id"), await _build_registry_from_json(body)


@router.post("/interpret")
async def interpret(request: Request) -> dict:
    """Ask a natural-language question over an inventory of plasmids.

    Body:
        question: str
        plasmids: [{plasmid_id, name, sequence}]   # required
        plasmid_id: str                            # optional, scopes the answer

    Returns:
        {ok, answer, citations: [...], trace: [...]}
    """
    import agent_v2  # noqa: F401 — path shim
    from splicify_api.annotation_cache import annotate_llm_cached
    from agent_v2.interpreter.plasmid_registry import PlasmidRegistry
    from agent_v2.interpreter.agent import run_interpreter

    body = await request.json()
    question = (body.get("question") or "").strip()
    plasmids = body.get("plasmids") or []
    if not question:
        return {"ok": False, "error": "question is required"}
    if not isinstance(plasmids, list) or not plasmids:
        return {"ok": False, "error": "plasmids list is required"}

    registry = PlasmidRegistry()
    for entry in plasmids:
        if not isinstance(entry, dict):
            continue
        pid = entry.get("plasmid_id") or entry.get("id")
        seq = entry.get("sequence")
        if not pid or not seq:
            continue
        try:
            envelope = await annotate_llm_cached(seq, circular=bool(entry.get("circular", True)))
            envelope.setdefault("sequence", seq)
            registry.register(pid, envelope, name=entry.get("name") or pid)
        except Exception as e:
            return {"ok": False, "error": f"annotation failed for {pid}: {type(e).__name__}: {e}"}

    if registry.n() == 0:
        return {"ok": False, "error": "no valid plasmid entries in body"}

    try:
        result = await run_interpreter(question, registry)
    except Exception as e:
        return {"ok": False, "error": f"interpreter failed: {type(e).__name__}: {e}"}

    return {
        "ok": True,
        "answer": result.answer,
        "citations": result.citations,
        "trace": result.trace,
        "n_tool_calls": result.n_tool_calls,
        "n_plasmids": registry.n(),
    }


@router.post("/chat-stream")
async def chat_stream(request: Request) -> StreamingResponse:
    """SSE streaming chat: emits `event: shorthand` as soon as the triage
    classifier finishes, then `event: envelope` with the full chat envelope.

    Body (JSON):
      {message, session_id?, target_genbank?, inventory_genbank?}
    """
    message, session_id, registry = await _parse_chat_body(request)

    queue: asyncio.Queue = asyncio.Queue()

    async def on_triage(tr: Any) -> None:
        await queue.put(("shorthand", {
            "intent":           tr.intent,
            "shorthand":        tr.shorthand,
            "is_new_topic":     tr.is_new_topic,
            "rejection_reason": tr.rejection_reason,
        }))

    async def on_tool_event(evt: dict) -> None:
        await queue.put(("tool_call", evt))

    async def runner() -> None:
        try:
            envelope = await run_orchestrator(
                message, registry,
                session_id=session_id,
                on_triage=on_triage,
                on_tool_event=on_tool_event,
            )
            await queue.put(("envelope", envelope))
        except Exception as e:
            await queue.put(("error", {"error": f"{type(e).__name__}: {e}"}))
        finally:
            await queue.put(None)

    task = asyncio.create_task(runner())

    async def event_stream():
        # First write: a comment line + 2 KB padding to push past
        # Cloudflare / Vercel buffer thresholds so the very first real
        # event reaches the browser without delay. CDN buffering
        # otherwise holds events until ~4 KB of body has accumulated,
        # which can lock the frontend on "Thinking…" for minutes.
        yield ":start\n" + (" " * 2048) + "\n\n"
        try:
            while True:
                # Heartbeat every 10 s if no agent event has fired yet.
                # Comment lines (lines starting with ':') are ignored by
                # the SSE spec but force flush past CDN buffers and keep
                # the connection alive past idle-timeout cuts.
                try:
                    item = await asyncio.wait_for(queue.get(), timeout=10.0)
                except asyncio.TimeoutError:
                    yield ":heartbeat\n\n"
                    continue
                if item is None:
                    break
                event_name, data = item
                yield f"event: {event_name}\ndata: {json.dumps(data)}\n\n"
        finally:
            await task  # propagate any unawaited exceptions

    return StreamingResponse(event_stream(), media_type="text/event-stream",
                              headers={
                                  "Cache-Control": "no-cache, no-transform",
                                  "Connection": "keep-alive",
                                  "X-Accel-Buffering": "no",
                              })
