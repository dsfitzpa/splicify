"""CRISPR guide-design agent stub.

Until a full guide-design agent is built, return a fixed 'coming soon' reply
that redirects to existing v1 tools (CRISPR guide-design modal,
sgRNA Golden Gate chat intent).
"""
from __future__ import annotations

from typing import Any


COMING_SOON_REPLY = (
    "**CRISPR guide design is coming soon as a dedicated agent.**\n\n"
    "In the meantime, the existing AI Plasmid Design tools cover most "
    "guide workflows:\n\n"
    "- **CRISPR guide design modal** in the main app — Doench 2014 + heuristic "
    "scoring, drag a `.gb` in and pick a CDS to target.\n"
    "- **sgRNA Golden Gate cloning** — type a Golden Gate guide-cloning prompt "
    "in the regular chat (e.g. \"clone three guides for EMX1 into pX330\") and "
    "the existing handler will design + assemble.\n\n"
    "I'll surface CRISPR guide design as a first-class agent with tool-use "
    "orchestration once the plasmid-cloning agent is stable in production."
)


def respond(user_message: str = "") -> dict[str, Any]:
    """Return the /agent_v2/chat-shaped envelope for a CRISPR guide intent."""
    return {
        "ok": True,
        "reply": COMING_SOON_REPLY,
        "files": None,
        "viz": None,
        "agent_trace": [
            {"iteration": 0, "tool": "stub_crispr",
             "args_summary": "(none)", "result_keys": ["reply"]}
        ],
        "n_tool_calls": 0,
        "error": None,
    }
