"""Summarizer subagent.

Takes the Main Agent's draft reply + emitted files metadata + the three
ExploreFinding digests + decisions ledger and produces a polished, concise
plain-language reply for the user. No tools, fresh context, single Anthropic
call.

Falls back to the Main Agent's draft on any error (empty response, parse
failure, etc.) so a flaky summarizer never breaks the pipeline.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Optional, Sequence

from agent_v2.explore import ExploreFinding


@dataclass
class SummaryResult:
    reply: str
    raw_response: Optional[Any] = None
    used_fallback: bool = False


SYSTEM_PROMPT = """You are the Summarizer for an AI molecular-biology agent. The Main Agent has finished an assembly workflow. Your job: produce a concise, friendly plain-language reply that the user reads as the final answer.

Cover:
- The assembly method used (open with it in bold).
- The single most important design choice + the alternative considered.
- The four output files (assembled.gb, parts_order.csv, protocol.csv, workflow_trace.txt) — one short line per file naming what's in it.
- 2-3 recommended wet-lab next steps.

Constraints:
- Markdown OK.
- <=350 words.
- No raw DNA. No tool-call commentary. No JSON.
- If the Main Agent's draft says "max iterations reached" or otherwise indicates failure, say so honestly and recommend a smaller-scope retry."""


def _format_findings(findings: Sequence[ExploreFinding]) -> str:
    if not findings:
        return "(no findings)"
    out = []
    for f in findings:
        summary = (f.summary_md or "").strip().split("\n")[0][:200]
        out.append(f"- {f.role}: {summary}")
    return "\n".join(out)


def _format_files(files: Optional[list[dict[str, Any]]]) -> str:
    if not files:
        return "(no files)"
    return ", ".join(f.get("fileName", "?") for f in files)


def _format_decisions(decisions: Optional[list[dict[str, Any]]]) -> str:
    if not decisions:
        return "(none recorded)"
    out = []
    for d in decisions:
        choice = d.get("choice", "?")
        alt = d.get("alternative")
        reason = d.get("reason", "")
        if alt:
            out.append(f"- chose {choice} over {alt}: {reason}")
        else:
            out.append(f"- {choice}: {reason}")
    return "\n".join(out)


async def run_summarizer(
    user_message: str,
    main_reply: str,
    files: Optional[list[dict[str, Any]]],
    findings: Sequence[ExploreFinding],
    decisions: Optional[list[dict[str, Any]]] = None,
    *,
    client: Any = None,
    model: Optional[str] = None,
) -> SummaryResult:
    if client is None:
        import anthropic
        client = anthropic.Anthropic()
    model = model or os.getenv("AGENT_MODEL", "claude-sonnet-4-6")

    user_block = (
        f"User prompt:\n{user_message}\n\n"
        f"Main agent draft (rewrite this in your own voice if too verbose):\n"
        f"{main_reply or '(empty)'}\n\n"
        f"Emitted files: {_format_files(files)}\n\n"
        f"Decisions ledger:\n{_format_decisions(decisions)}\n\n"
        f"Findings (one-liners):\n{_format_findings(findings)}\n\n"
        "Produce the final user-facing reply now."
    )

    try:
        response = client.messages.create(
            model=model,
            max_tokens=2048,
            system=SYSTEM_PROMPT,
            messages=[{"role": "user", "content": user_block}],
        )
        text = "".join(getattr(b, "text", "") for b in response.content
                       if getattr(b, "type", None) == "text").strip()
        if not text:
            return SummaryResult(
                reply=main_reply or "(no reply produced)",
                raw_response=response, used_fallback=True,
            )
        return SummaryResult(reply=text, raw_response=response, used_fallback=False)
    except Exception as e:
        return SummaryResult(
            reply=main_reply or f"(summarizer error: {type(e).__name__})",
            raw_response={"error": str(e)}, used_fallback=True,
        )
