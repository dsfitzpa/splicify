"""Shared ReAct loop + types for the three Explore-style subagents.

Each subagent (PartScout, TargetBuilder, MethodRouter) supplies its own
narrow tool roster + system prompt; the loop body, JSON parsing, and
no-raw-DNA discipline live here.

Mirrors the Claude Code blog-post pattern: each Explore subagent returns a
digested ExploreFinding — never raw DNA, never full annotation lists. The
Plan agent later receives only these digests + the user prompt.
"""
from __future__ import annotations
import asyncio

import json
import os
from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass
class ExploreFinding:
    role: str  # "part_scout" | "target_builder" | "method_router"
    summary_md: str  # human-readable digest; <=400 tokens by convention
    key_facts: dict[str, Any] = field(default_factory=dict)
    references: list[str] = field(default_factory=list)  # attachment_ids
    trace: list[dict[str, Any]] = field(default_factory=list)


def _block_to_dict(b: Any) -> dict[str, Any]:
    t = getattr(b, "type", None)
    if t == "text":
        return {"type": "text", "text": getattr(b, "text", "")}
    if t == "tool_use":
        return {"type": "tool_use", "id": b.id, "name": b.name, "input": b.input}
    return {}


def _strip_fences(text: str) -> str:
    text = text.strip()
    if text.startswith("```"):
        nl = text.find("\n")
        if nl != -1:
            text = text[nl + 1:]
        if text.endswith("```"):
            text = text[:-3]
    return text.strip()


def _summarize_args(args: dict[str, Any]) -> str:
    return ", ".join(f"{k}={str(v)[:30]}" for k, v in (args or {}).items())[:120]


def _digest_annotation(att_result: dict[str, Any]) -> dict[str, Any]:
    """Compact annotation summary for embedding in a finding's
    key_facts. Drops bulk sequence + feature-list bodies; keeps the
    counts, module-type set, top feature names, and length so the
    Main agent can read it directly without calling annotate_attachment."""
    feats = att_result.get("features") or att_result.get("annotations") or []
    modules = att_result.get("modules") or []
    cf = att_result.get("cloning_features") or []
    if isinstance(cf, dict):
        cf = cf.get("features") or []
    feature_names = [
        (f.get("name") or "")
        for f in feats[:30]
        if isinstance(f, dict) and f.get("name")
    ]
    module_types = sorted({
        (m.get("module_type") or m.get("type"))
        for m in modules
        if isinstance(m, dict) and (m.get("module_type") or m.get("type"))
    })
    return {
        "attachment_id": att_result.get("attachment_id"),
        "name": att_result.get("name"),
        "length_bp": att_result.get("length_bp") or att_result.get("length"),
        "n_features": len(feats),
        "n_modules": len(modules),
        "n_cloning_features": len(cf) if isinstance(cf, list) else 0,
        "module_types": module_types,
        "feature_names_top": feature_names,
    }


def _finalise_finding(*, role: str, text: str, registry, trace,
                      annotation_index: dict[str, Any],
                      empty_fact_fallback: bool = False) -> ExploreFinding:
    cleaned = _strip_fences(text or "")
    refs = [a["attachment_id"] for a in registry.public_summary()]
    summaries = [_digest_annotation(v) for v in annotation_index.values()]
    try:
        data = json.loads(cleaned) if cleaned else {}
        kf = dict(data.get("key_facts", {}) or {})
        # Always merge in the deterministic annotation summaries so the
        # Main agent doesn't have to re-run annotate_attachment.
        if summaries:
            kf.setdefault("annotation_summaries", summaries)
        return ExploreFinding(
            role=role,
            summary_md=str(data.get("summary_md", ""))[:8000] or (text.strip()[:8000] if text else "(no summary)"),
            key_facts=kf,
            references=refs,
            trace=trace,
        )
    except (json.JSONDecodeError, ValueError):
        kf = {}
        if summaries:
            kf["annotation_summaries"] = summaries
        return ExploreFinding(
            role=role,
            summary_md=(text or "").strip()[:8000] or "(no summary)",
            key_facts=kf,
            references=refs,
            trace=trace,
        )


def _build_finding(role: str, text: str, registry, trace) -> ExploreFinding:
    """Backwards-compatible shim used by callers that didn't pass an
    annotation_index. Just forwards to _finalise_finding with an empty
    index."""
    return _finalise_finding(role=role, text=text, registry=registry,
                              trace=trace, annotation_index={})


async def run_explore_subagent(
    *,
    role: str,
    user_message: str,
    registry: Any,
    tools: list[dict[str, Any]],
    system_prompt: str,
    client: Any = None,
    dispatch_fn: Any = None,
    max_iters: int = 8,
    model: Optional[str] = None,
    on_tool_event: Any = None,
) -> ExploreFinding:
    """Generic Anthropic ReAct loop for an Explore subagent.

    Each subagent passes its narrow `tools` + `system_prompt`; the loop body
    is shared. `dispatch_fn` defaults to v1's async `dispatch_tool` which
    already runs `_strip_sequences` on every result.
    """
    if client is None:
        import anthropic
        client = anthropic.Anthropic()
    if dispatch_fn is None:
        import agent_v2  # noqa: F401 — triggers path shim
        from splicify_api.agent.agent_tools import dispatch_tool as _v1_dispatch
        dispatch_fn = _v1_dispatch
    model = model or os.getenv("AGENT_MODEL", "claude-sonnet-4-6")

    messages: list[dict[str, Any]] = [{
        "role": "user",
        "content": (
            f"User message: {user_message}\n\n"
            f"Registered attachments: {json.dumps(registry.public_summary())}"
        ),
    }]
    trace: list[dict[str, Any]] = []
    # annotation_index[attachment_id] = digested annotate_attachment result.
    # Surfaced in the finding so the Main agent skips re-annotation.
    annotation_index: dict[str, Any] = {}

    cached_system = [{"type": "text", "text": system_prompt,
                       "cache_control": {"type": "ephemeral"}}]
    cached_tools = [{**t} for t in tools]
    if cached_tools:
        cached_tools[-1] = {**cached_tools[-1],
                              "cache_control": {"type": "ephemeral"}}
    if messages and messages[0]["role"] == "user" and isinstance(messages[0].get("content"), str):
        messages[0] = {
            "role": "user",
            "content": [{"type": "text", "text": messages[0]["content"],
                          "cache_control": {"type": "ephemeral"}}],
        }

    for i in range(max_iters):
        response = client.messages.create(
            model=model,
            max_tokens=4096,
            system=cached_system,
            tools=cached_tools,
            messages=messages,
        )
        tool_uses = [b for b in response.content
                     if getattr(b, "type", None) == "tool_use"]
        if not tool_uses:
            text = "".join(getattr(b, "text", "") for b in response.content
                           if getattr(b, "type", None) == "text")
            return _finalise_finding(role=role, text=text, registry=registry,
                                      trace=trace, annotation_index=annotation_index)

        messages.append({
            "role": "assistant",
            "content": [_block_to_dict(b) for b in response.content],
        })

        # Globally-unique iteration index across this subagent's loop so the
        # frontend can correlate start/end events.
        next_iter = len(trace)
        # Fire `start` events so the frontend can render a `running` row for
        # any subagent tool call (especially find_genomic_record, which lives
        # in TargetLocator). Iteration numbers are offset by -1000 to avoid
        # colliding with Main agent indices.
        for offset, b in enumerate(tool_uses):
            if on_tool_event is not None:
                try:
                    res_evt = on_tool_event({
                        "phase": "start",
                        "iteration": -1000 - (next_iter + offset),
                        "tool": b.name,
                        "input": b.input or {},
                        "subagent_role": role,
                        "t_unix": __import__("time").time(),
                    })
                    if hasattr(res_evt, "__await__"):
                        await res_evt
                except Exception:
                    pass

        # Parallel tool dispatch — see the matching block in main_agent.py.
        async def _dispatch_one(block, idx):
            t0 = __import__("time").monotonic()
            result = await dispatch_fn(block.name, block.input or {}, registry)
            return block, idx, result, int((__import__("time").monotonic() - t0) * 1000)

        dispatched = await asyncio.gather(*(_dispatch_one(b, idx) for idx, b in enumerate(tool_uses)))

        tool_results: list[dict[str, Any]] = []
        for b, off, result, elapsed_ms in dispatched:
            iter_idx = -1000 - (next_iter + off)
            tool_results.append({
                "type": "tool_result",
                "tool_use_id": b.id,
                "content": json.dumps(result),
            })
            trace.append({
                "iteration": i,
                "tool": b.name,
                "args_summary": _summarize_args(b.input or {}),
                "result_keys": list(result.keys())[:8] if isinstance(result, dict) else [],
            })
            # Capture annotation results so the Main agent receives an
            # annotation digest without having to re-run annotate_attachment.
            if b.name == "annotate_attachment" and isinstance(result, dict):
                aid = result.get("attachment_id") or (b.input or {}).get("attachment_id")
                if aid:
                    annotation_index[aid] = result
            # Fire `end` event — include any file envelope (e.g.
            # find_genomic_record's retrieved .gb) so the frontend renders
            # it mid-stream without waiting for the full pipeline to finish.
            if on_tool_event is not None:
                try:
                    evt_payload = {
                        "phase": "end",
                        "iteration": iter_idx,
                        "tool": b.name,
                        "elapsed_ms": elapsed_ms,
                        "ok": (result.get("ok") if isinstance(result, dict) else None),
                        "subagent_role": role,
                        "t_unix": __import__("time").time(),
                    }
                    if isinstance(result, dict) and isinstance(result.get("file"), dict):
                        evt_payload["file"] = result["file"]
                    res_evt = on_tool_event(evt_payload)
                    if hasattr(res_evt, "__await__"):
                        await res_evt
                except Exception:
                    pass
        messages.append({"role": "user", "content": tool_results})

    return _finalise_finding(
        role=role, text="(max iterations reached without final summary)",
        registry=registry, trace=trace, annotation_index=annotation_index,
        empty_fact_fallback=True,
    )
