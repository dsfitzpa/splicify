"""Triage agent — one-shot Sonnet 4.6 classifier.

Returns a TriageResult that tells the orchestrator which path to take:
  PLASMID_CLONING -> main agent loop
  CRISPR_GUIDE    -> stub_crispr.respond
  REJECT          -> rejection.respond(shorthand, reason)

Also emits a <=8-word shorthand summary that the frontend renders in an
"actively writing" bubble before the heavy work starts, and an is_new_topic
flag that lets the orchestrator decide whether to keep or reset the session.

The Anthropic client is injectable so tests can mock it without touching the
network. In production, leave `client=None` to use the SDK default
(reads ANTHROPIC_API_KEY from the environment).
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Any, Literal, Optional

import anthropic

Intent = Literal["PLASMID_CLONING", "PLASMID_QA", "CRISPR_GUIDE", "REJECT"]


@dataclass
class TriageResult:
    intent: Intent
    shorthand: str
    is_new_topic: bool
    rejection_reason: Optional[str] = None
    raw_response: Optional[dict] = None


SYSTEM_PROMPT = """You are the triage classifier for an AI molecular-biology agent. Classify the user's request into exactly one of three intents and emit structured JSON. No prose outside the JSON.

INTENTS:
- PLASMID_CLONING - designing or assembling a plasmid (Gibson, Gateway, restriction cloning, Golden Gate, site-directed mutagenesis, sgRNA Golden Gate, target-from-inventory routing). Anything that PRODUCES or MODIFIES a plasmid map. Keywords: "build", "assemble", "design a vector to...", "clone X into Y", "make a plasmid that expresses".
- PLASMID_QA - questions ABOUT a plasmid OR a request to FIND a named/described plasmid (in the uploaded inventory or from public repositories like Addgene). Use this for:
    * coordinate lookups ("Where is Cas9?", "Coordinates of the U6 promoter")
    * amino-acid lookups ("What is the 55th aa in Cas9?")
    * regulatory queries ("What promoter and polyA drives Cas9?")
    * application classification ("What is this plasmid used for?")
    * feature listing, inventory comparison
    * **find/search/identify/suggest a plasmid that has X** — e.g. "Find a pegRNA cloning plasmid with a tevopreQ1 motif", "Suggest a Cas9 lentiviral backbone with puromycin selection", "What is the canonical Addgene plasmid for sgRNA Golden Gate cloning?". These are read-only lookups; the QA pipeline will search Addgene and auto-annotate the candidate, no cloning workflow is needed.
  Even when the user says the word "cloning" inside a noun phrase like "pegRNA cloning plasmid" — that's QA (they want to FIND a backbone), not CLONING (which would mean "design me a custom plasmid"). The distinction is verb: "find / suggest / which / what is" → QA; "build / design / assemble" → CLONING.
- CRISPR_GUIDE - designing CRISPR guides (sgRNA selection, off-target scoring, PAM choice). Note: sgRNA *cloning* into a Golden Gate vector is PLASMID_CLONING, not CRISPR_GUIDE - the line is "design which guide" vs "build the guide-bearing plasmid".
- REJECT - anything outside molecular biology / cloning / guide design. General Claude questions, code, casual chat, unrelated science.

OUTPUT: a single JSON object with exactly these fields:
{
  "intent": "PLASMID_CLONING" | "PLASMID_QA" | "CRISPR_GUIDE" | "REJECT",
  "shorthand": "<<=8 words echoing the user's request - appears in 'actively writing' bubble>",
  "is_new_topic": true | false,
  "rejection_reason": null | "<one short sentence, only when intent=REJECT>"
}

is_new_topic should be true if the prompt does not appear to follow up on a previous turn in the conversation. If you have no history at all, default to true.

Output JSON only. No commentary, no markdown fences."""


def triage(
    user_message: str,
    has_attachments: bool = False,
    *,
    client: Optional[Any] = None,
    model: Optional[str] = None,
) -> TriageResult:
    """One-shot triage; never raises on parse failure (returns REJECT with reason)."""
    client = client if client is not None else anthropic.Anthropic()
    model = model or os.getenv("AGENT_MODEL", "claude-sonnet-4-6")
    user_block = (
        f"User message: {user_message}\n"
        f"Has attachments: {has_attachments}\n\n"
        "Return only the JSON object."
    )

    response = client.messages.create(
        model=model,
        max_tokens=300,
        system=SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_block}],
    )

    text = "".join(
        getattr(b, "text", "") for b in response.content
        if getattr(b, "type", None) == "text"
    )
    try:
        data = json.loads(_strip_fences(text))
    except (json.JSONDecodeError, ValueError):
        return TriageResult(
            intent="REJECT",
            shorthand="(unparseable)",
            is_new_topic=True,
            rejection_reason=f"triage parse failure: {text[:120]!r}",
            raw_response={"text": text},
        )

    intent = data.get("intent")
    if intent not in ("PLASMID_CLONING", "PLASMID_QA", "CRISPR_GUIDE", "REJECT"):
        return TriageResult(
            intent="REJECT",
            shorthand=str(data.get("shorthand", ""))[:80],
            is_new_topic=bool(data.get("is_new_topic", True)),
            rejection_reason=f"unknown intent {intent!r}",
            raw_response=data,
        )

    return TriageResult(
        intent=intent,
        shorthand=str(data.get("shorthand", ""))[:80],
        is_new_topic=bool(data.get("is_new_topic", True)),
        rejection_reason=data.get("rejection_reason"),
        raw_response=data,
    )


def _strip_fences(text: str) -> str:
    text = text.strip()
    if text.startswith("```"):
        first_newline = text.find("\n")
        if first_newline != -1:
            text = text[first_newline + 1:]
        if text.endswith("```"):
            text = text[:-3]
    return text.strip()
