"""
Anthropic Claude tool-use loop for AIPlasmidDesign (v2).

The agent NEVER receives raw DNA. Each attached plasmid is registered in
an AttachmentRegistry; for MCQ benchmarks the unredacted choice texts can
also be registered server-side so compare_to_choice can hash-match without
exposing the underlying sequence to Claude.
"""
from __future__ import annotations

import json
import logging
import os
from typing import Any, Dict, List, Optional

from .tool_schemas import ALL_TOOLS, SYSTEM_PROMPT
from .agent_tools import (
    AttachmentRegistry, Attachment, dispatch_tool,
    extract_seq_from_genbank, extract_name_from_genbank,
)

logger = logging.getLogger("agent.loop")

DEFAULT_MODEL  = os.environ.get("AGENT_MODEL", "claude-opus-4-7")
MAX_ITERATIONS = int(os.environ.get("AGENT_MAX_ITERATIONS", "24"))
MAX_TOKENS     = int(os.environ.get("AGENT_MAX_TOKENS", "16384"))


def _build_registry(
    target_gb: Optional[str], inventory_gbs: List[str],
    choices: Optional[List[Dict[str, str]]] = None,
) -> AttachmentRegistry:
    reg = AttachmentRegistry()
    next_id = 1
    if target_gb:
        seq = extract_seq_from_genbank(target_gb)
        if seq:
            aid = f"att_{next_id}"; next_id += 1
            reg.add(Attachment(
                attachment_id=aid,
                name=extract_name_from_genbank(target_gb, aid),
                sequence=seq, circular=True, role="target",
            ))
    for gb in inventory_gbs:
        seq = extract_seq_from_genbank(gb)
        if not seq:
            continue
        aid = f"att_{next_id}"; next_id += 1
        reg.add(Attachment(
            attachment_id=aid,
            name=extract_name_from_genbank(gb, aid),
            sequence=seq, circular=True, role="inventory",
        ))
    if choices:
        for c in choices:
            letter = (c or {}).get("letter")
            text = (c or {}).get("text", "")
            if letter:
                reg.add_choice(letter, text)
    return reg


def _build_user_message(user_text: str, reg: AttachmentRegistry) -> str:
    summary = reg.public_summary()
    parts = [user_text or ""]
    if summary:
        lines = ["", "Attached plasmids (referenced by attachment_id; raw sequences NOT shown):"]
        for s in summary:
            lines.append(f"  - {s['attachment_id']}  name={s['name']}  "
                         f"length={s['length_bp']} bp  {s['topology']}  role={s['role']}")
        parts.append("\n".join(lines))
    if reg.choices:
        parts.append("\nMCQ choices are registered server-side; use compare_to_choice with a letter (A/B/C/D...) to hash-match an attachment against a choice.")
    return "\n".join(p for p in parts if p).strip()


async def run_agent(
    user_message: str,
    target_gb: Optional[str] = None,
    inventory_gbs: Optional[List[str]] = None,
    choices: Optional[List[Dict[str, str]]] = None,
    model: str = DEFAULT_MODEL,
) -> Dict[str, Any]:
    inventory_gbs = inventory_gbs or []

    try:
        import anthropic
    except ImportError:
        return {"reply": "AI Agent unavailable: `anthropic` not installed.",
                "trace": [], "n_tool_calls": 0, "error": "anthropic_not_installed"}

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        return {"reply": "AI Agent unavailable: ANTHROPIC_API_KEY not set.",
                "trace": [], "n_tool_calls": 0, "error": "missing_api_key"}

    client = anthropic.AsyncAnthropic(api_key=api_key)

    registry = _build_registry(target_gb, inventory_gbs, choices)
    initial_user = _build_user_message(user_message, registry)
    messages: List[Dict[str, Any]] = [{"role": "user", "content": initial_user}]
    trace: List[Dict[str, Any]] = []
    n_tool_calls = 0

    for iteration in range(MAX_ITERATIONS):
        try:
            resp = await client.messages.create(
                model=model, max_tokens=MAX_TOKENS,
                system=SYSTEM_PROMPT, tools=ALL_TOOLS, messages=messages,
            )
        except Exception as e:
            logger.exception("anthropic call failed")
            return {"reply": f"AI Agent error: {type(e).__name__}: {e}",
                    "trace": trace, "n_tool_calls": n_tool_calls,
                    "error": "anthropic_call_failed", "_registry": registry}

        messages.append({"role": "assistant", "content": resp.content})

        if resp.stop_reason == "end_turn":
            final_text = "".join(b.text for b in resp.content
                                 if getattr(b, "type", None) == "text")
            trace.append({"iteration": iteration, "stop_reason": "end_turn",
                          "text": final_text[:500]})
            return {"reply": final_text.strip(), "trace": trace,
                    "n_tool_calls": n_tool_calls, "_registry": registry}

        if resp.stop_reason != "tool_use":
            trace.append({"iteration": iteration, "stop_reason": resp.stop_reason})
            return {"reply": f"AI Agent stopped unexpectedly (stop_reason={resp.stop_reason}).",
                    "trace": trace, "n_tool_calls": n_tool_calls, "_registry": registry}

        tool_results: List[Dict[str, Any]] = []
        for block in resp.content:
            if getattr(block, "type", None) != "tool_use":
                continue
            tool_name = block.name
            tool_args = block.input or {}
            n_tool_calls += 1
            logger.info("agent tool call: %s args_keys=%s",
                        tool_name, list(tool_args.keys()))
            try:
                result = await dispatch_tool(tool_name, tool_args, registry)
            except Exception as e:
                logger.exception("tool %s raised", tool_name)
                result = {"error": f"{type(e).__name__}: {e}"}

            trace.append({
                "iteration": iteration, "tool": tool_name,
                "args_summary": {
                    k: (v[:60] + "..." if isinstance(v, str) and len(v) > 60 else v)
                    for k, v in tool_args.items()},
                "result_keys": list(result.keys()) if isinstance(result, dict) else None,
            })
            tool_results.append({
                "type": "tool_result", "tool_use_id": block.id,
                "content": json.dumps(result, default=str)[:50000],
            })

        if not tool_results:
            return {"reply": "AI Agent stopped: tool_use stop_reason but no tool_use blocks.",
                    "trace": trace, "n_tool_calls": n_tool_calls, "_registry": registry}
        messages.append({"role": "user", "content": tool_results})

    return {"reply": f"AI Agent reached max iterations ({MAX_ITERATIONS}).",
            "trace": trace, "n_tool_calls": n_tool_calls, "error": "max_iterations", "_registry": registry}
