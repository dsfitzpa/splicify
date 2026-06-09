"""CRISPR-design clarification responder.

When the user's prompt is missing the information the pipeline needs to
design a sensible guide — no target named, no PAM-resolvable feature,
no edit description for a pegRNA request — we should NOT silently pick a
default. We should ask the user for the missing piece.

Wired in `_run_crispr_pipeline`: after the two Explore subagents return,
the orchestrator collects every `missing_info` string the subagents
flagged. If the list is non-empty, this responder builds the envelope
(no Plan agent, no Main agent, no LLM cost) and the user gets a focused
follow-up question.

Mirrors `agent_v2/rejection.py` in shape: envelope-only, no tool calls,
ok=True so the frontend renders the reply normally.
"""
from __future__ import annotations

from typing import Any, Iterable


def respond(missing_info: Iterable[str], *, user_message: str = "",
             shorthand: str = "") -> dict[str, Any]:
    """Build a clarification chat envelope listing every missing piece.

    The reply is a short markdown response with one bullet per missing
    item — no commentary outside the bullets. The intent is that the
    user reads it, fills in the gap, and resubmits.
    """
    items = [s for s in (missing_info or []) if s]
    if not items:
        items = ["The CRISPR-design request is missing information, but the "
                  "agent could not identify what specifically is missing."]

    bullets = "\n".join(f"- {s}" for s in items)
    user_echo = (shorthand or user_message or "").strip()
    user_block = f"You asked: _{user_echo}_\n\n" if user_echo else ""

    reply = (
        "## I need a bit more information to design your guides\n\n"
        + user_block
        + "Before I can run the CRISPR-design pipeline, please answer the "
        "following:\n\n"
        + bullets
        + "\n\nResubmit with the missing details and I'll run the full "
        "design (sgRNA / pegRNA + NGS + Sanger primers + cloning oligos + "
        "the protocol). If you want to broaden the scope (e.g. design "
        "guides across the entire CDS instead of one residue), say so "
        "explicitly and I'll honour it."
    )

    return {
        "ok": True,
        "reply": reply,
        "files": None,
        "viz": None,
        "agent_trace": [
            {"iteration": 0, "tool": "crispr_clarification",
             "args_summary": f"n_missing={len(items)}",
             "result_keys": ["reply"]}
        ],
        "n_tool_calls": 0,
        "error": None,
        "workflow": "crispr",
        "missing_info": list(items),
    }
