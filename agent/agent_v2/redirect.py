"""Workflow redirect sentinel.

Any pipeline (cloning / qa / crispr) can return a `WorkflowRedirect`
instead of a finished envelope to signal "the user's request would be
better served by a different workflow." The orchestrator catches the
sentinel and re-runs with the new intent (capped to avoid loops).

This avoids the costly failure mode where triage misclassifies a
request and the wrong full pipeline runs to completion before the
user discovers the answer doesn't match what they asked.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


VALID_REDIRECT_TARGETS = {"PLASMID_CLONING", "PLASMID_QA", "CRISPR_GUIDE", "REJECT"}


@dataclass
class WorkflowRedirect:
    """Sentinel envelope. The orchestrator inspects every pipeline
    result and reroutes when this is present."""
    target_intent: str                       # one of VALID_REDIRECT_TARGETS
    reason: str
    notes: str = ""                          # surfaced in the final reply for transparency
    findings: dict[str, Any] = field(default_factory=dict)   # optional context to seed the next pipeline
    source_intent: str = ""                  # filled by the orchestrator

    def to_envelope(self) -> dict[str, Any]:
        """Wrap as a chat-envelope-shaped dict so callers can pattern-match
        on the redirect_to key."""
        return {
            "ok": True,
            "redirect_to": self.target_intent,
            "redirect_reason": self.reason,
            "redirect_notes": self.notes,
            "redirect_findings": dict(self.findings),
            "source_intent": self.source_intent,
            "files": None, "viz": None, "agent_trace": [],
            "n_tool_calls": 0, "error": None,
        }


def is_redirect_envelope(envelope: Any) -> bool:
    return isinstance(envelope, dict) and bool(envelope.get("redirect_to"))


def merge_redirect_into_reply(envelope: dict[str, Any]) -> dict[str, Any]:
    """When a redirect's downstream pipeline finishes, prepend a small
    notice to the reply so the user sees the workflow switch."""
    if not envelope.get("redirected_from"):
        return envelope
    notice = (
        f"_(Note: started in {envelope['redirected_from']} workflow; switched to "
        f"{envelope.get('intent')} because: {envelope.get('redirect_reason', '')})_\n\n"
    )
    envelope = dict(envelope)
    envelope["reply"] = notice + (envelope.get("reply") or "")
    return envelope
