"""Off-topic rejection: load the soft-tone template and wrap it as a chat envelope.

The triage agent populates `shorthand` (≤8 words echoing what the user asked)
so the rejection at least mirrors the user's prompt back to them.
"""
from __future__ import annotations

import pathlib
from typing import Any


_TEMPLATE_PATH = pathlib.Path(__file__).parent / "prompts" / "rejection_template.md"


def respond(shorthand: str = "", reason: str | None = None) -> dict[str, Any]:
    template = _TEMPLATE_PATH.read_text()
    reply = template.replace("{shorthand}", shorthand or "(no summary available)")
    return {
        "ok": True,
        "reply": reply,
        "files": None,
        "viz": None,
        "agent_trace": [
            {"iteration": 0, "tool": "rejection_template",
             "args_summary": f"reason={reason or 'off-topic'}",
             "result_keys": ["reply"]}
        ],
        "n_tool_calls": 0,
        "error": None,
    }
