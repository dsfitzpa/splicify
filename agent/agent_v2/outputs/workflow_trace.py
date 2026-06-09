"""emit_workflow_trace — fourth and final output emitter.

Writes a plain-text workflow_trace.txt with these sections:
  Header (session, turn, method, product, timestamp)
  --- AGENT TRACE ---       per-iteration tool / args / result_keys
  --- PLAN.MD ---           final plan with checked items
  --- DECISIONS LEDGER ---  chosen vs alternative, with reason
  --- DESIGN VERIFICATION ---  (optional; from verifier output)
  --- EXPLORE FINDINGS ---  (optional; the 3 subagent summaries)
"""
from __future__ import annotations

import base64
import json
import pathlib
import time
from typing import Any, Optional
from agent_v2.outputs import prefixed_filename, derive_descriptor


def _hr(label: str) -> str:
    return f"\n--- {label} ---\n"


def _format_trace(trace: list[dict[str, Any]]) -> str:
    if not trace:
        return "(no tool calls)\n"
    lines = []
    for entry in trace:
        i = entry.get("iteration", "?")
        tool = entry.get("tool", "?")
        args = entry.get("args_summary", "")
        keys = entry.get("result_keys", [])
        keys_str = ", ".join(keys) if keys else "(no keys)"
        lines.append(f"iter {i}: {tool}({args}) -> {keys_str}")
    return "\n".join(lines) + "\n"


def _format_decisions(decisions: list[dict[str, Any]]) -> str:
    if not decisions:
        return "(no decisions recorded)\n"
    lines = []
    for d in decisions:
        choice = d.get("choice", "?")
        alt = d.get("alternative")
        reason = d.get("reason", "")
        if alt:
            lines.append(f"- Chose {choice}; runner-up {alt}. Reason: {reason}")
        else:
            lines.append(f"- {choice}. Reason: {reason}")
    return "\n".join(lines) + "\n"


def _format_verifier(verifier: dict[str, Any]) -> str:
    if not verifier:
        return ""
    out = []
    passed = verifier.get("passed")
    if passed is True:
        out.append("Passed.")
    elif passed is False:
        out.append("FAILED.")
    else:
        out.append("(verifier not run)")
    warnings = verifier.get("warnings") or []
    if warnings:
        out.append("Warnings:")
        for w in warnings:
            if isinstance(w, dict):
                feat = w.get("feature_name") or w.get("feature") or "?"
                msg = w.get("remediation") or w.get("message") or json.dumps(w)
                out.append(f"  - {feat}: {msg}")
            else:
                out.append(f"  - {w}")
    else:
        out.append("Warnings: (none)")
    return "\n".join(out) + "\n"


def _format_findings(findings: list[dict[str, Any]]) -> str:
    if not findings:
        return ""
    lines = []
    for f in findings:
        role = f.get("role", "?")
        summary = f.get("summary_md", "").strip()
        key_facts = f.get("key_facts") or {}
        lines.append(f"[{role}] {summary}")
        if key_facts:
            lines.append(f"key_facts: {json.dumps(key_facts, default=str)}")
        lines.append("")
    return "\n".join(lines)


def _build_text(args: dict[str, Any]) -> str:
    session_id = args.get("session_id", "(no session)")
    turn_id = args.get("turn_id", "(no turn)")
    method = args.get("assembly_method", "(unknown)")
    product = args.get("product_attachment_id") or "(no product)"
    timestamp = args.get("timestamp") or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    parts: list[str] = []
    parts.append("=== agent_v2 workflow trace ===")
    parts.append(f"session: {session_id}")
    parts.append(f"turn:    {turn_id}")
    parts.append(f"written: {timestamp}")
    parts.append(f"method:  {method}")
    parts.append(f"product: {product}")

    parts.append(_hr("AGENT TRACE"))
    parts.append(_format_trace(args.get("agent_trace") or []))

    parts.append(_hr("PLAN.MD"))
    plan_md = args.get("plan_md") or "(no plan)"
    parts.append(plan_md.rstrip() + "\n")

    parts.append(_hr("DECISIONS LEDGER"))
    parts.append(_format_decisions(args.get("decisions") or []))

    verifier_block = _format_verifier(args.get("verifier") or {})
    if verifier_block:
        parts.append(_hr("DESIGN VERIFICATION"))
        parts.append(verifier_block)

    findings_block = _format_findings(args.get("findings") or [])
    if findings_block:
        parts.append(_hr("EXPLORE FINDINGS"))
        parts.append(findings_block)

    return "".join(parts).rstrip() + "\n"


async def emit_workflow_trace(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    text = _build_text(args)

    _descriptor = derive_descriptor(args)
    file_envelope = {
        "fileName": prefixed_filename("workflow_trace.txt", _descriptor),
        "dataBase64": base64.b64encode(text.encode("utf-8")).decode("ascii"),
    }

    written_path: Optional[str] = None
    if output_dir is not None:
        out = pathlib.Path(output_dir) / prefixed_filename("workflow_trace.txt", _descriptor)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
        written_path = str(out)

    return {
        "ok": True,
        "file": file_envelope,
        "n_chars": len(text),
        "n_trace_entries": len(args.get("agent_trace") or []),
        "n_decisions": len(args.get("decisions") or []),
        "written_path": written_path,
    }

# ---------------------------------------------------------------------------
# Server-side auto-emission (no LLM involvement)
# ---------------------------------------------------------------------------
async def auto_emit_workflow_trace(
    *,
    session_id: str,
    turn_id: str,
    main_result,
    findings,
    plan_md: str,
    output_dir=None,
    registry=None,
    assembly_method: str = "",
    product_attachment_id: str = "",
) -> dict:
    """Build emit_workflow_trace's args dict from the post-main-agent state
    and call the emitter directly. Used by the orchestrator after the
    Main agent returns so the LLM doesn't have to compose this big dict
    in a final tool_use turn (a known cause of the run pushing past
    Vercel's 5-minute SSE-proxy timeout).
    """
    import time as _time
    # Strip the heavy `file` blob from trace entries — emit_workflow_trace
    # only needs the iteration/tool/args/result_keys shape.
    trace_compact = []
    for entry in (main_result.trace or []):
        trace_compact.append({
            "iteration": entry.get("iteration"),
            "tool": entry.get("tool"),
            "args_summary": entry.get("args_summary", ""),
            "result_keys": entry.get("result_keys", []) or [],
        })

    # Pull a cached descriptor ONLY when it's bound to THIS session +
    # turn. _EMITTER_ARG_CACHE is process-wide; pulling "the most
    # recent" descriptor leaks CGAS-style names from a CRISPR session
    # into an unrelated cloning session's workflow_trace filename.
    cached_descriptor = None
    try:
        from agent_v2.tools import _EMITTER_ARG_CACHE
        # The cache key shape is session_id:turn_id (see tools.py's
        # _emitter_cache_key). Match exactly.
        bucket = _EMITTER_ARG_CACHE.get(f"{session_id}:{turn_id}")
        if isinstance(bucket, dict) and bucket.get("descriptor"):
            cached_descriptor = bucket["descriptor"]
    except Exception:
        pass

    args = {
        "session_id": session_id,
        "turn_id": turn_id,
        "assembly_method": assembly_method,
        "product_attachment_id": product_attachment_id,
        "descriptor": cached_descriptor,
        "timestamp": _time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime()),
        "agent_trace": trace_compact,
        "plan_md": plan_md or "",
        "decisions": [],
        "findings": [
            {"role": f.role, "summary_md": f.summary_md,
              "key_facts": dict(f.key_facts) if f.key_facts else {}}
            for f in (findings or [])
        ],
    }
    return await emit_workflow_trace(args, registry, output_dir=output_dir)
