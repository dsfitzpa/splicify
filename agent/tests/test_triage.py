"""Tests for the triage classifier (Anthropic client mocked)."""
import json
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from agent_v2 import triage as triage_mod


class _Block:
    type = "text"

    def __init__(self, text):
        self.text = text


class _Resp:
    def __init__(self, text):
        self.content = [_Block(text)]


class _Messages:
    def __init__(self, text_to_return):
        self._text = text_to_return
        self.calls = []

    def create(self, **kwargs):
        self.calls.append(kwargs)
        return _Resp(self._text)


class _Client:
    def __init__(self, text_to_return):
        self.messages = _Messages(text_to_return)


def _payload(intent, shorthand, is_new=True, reason=None):
    return json.dumps({
        "intent": intent,
        "shorthand": shorthand,
        "is_new_topic": is_new,
        "rejection_reason": reason,
    })


def test_triage_classifies_plasmid_cloning():
    res = triage_mod.triage(
        "Build me a Gibson assembly with hPGK-GFP into pUC19",
        client=_Client(_payload("PLASMID_CLONING", "Gibson: hPGK-GFP into pUC19")),
    )
    assert res.intent == "PLASMID_CLONING"
    assert "Gibson" in res.shorthand
    assert res.is_new_topic is True
    assert res.rejection_reason is None


def test_triage_classifies_crispr():
    res = triage_mod.triage(
        "design a CRISPR guide for EMX1",
        client=_Client(_payload("CRISPR_GUIDE", "guide design for EMX1")),
    )
    assert res.intent == "CRISPR_GUIDE"


def test_triage_rejects_off_topic():
    res = triage_mod.triage(
        "What's your favorite color?",
        client=_Client(_payload("REJECT", "favorite color question",
                                reason="outside molecular biology scope")),
    )
    assert res.intent == "REJECT"
    assert res.rejection_reason == "outside molecular biology scope"


def test_triage_strips_markdown_fences():
    fenced = "```json\n" + _payload("PLASMID_CLONING", "test") + "\n```"
    res = triage_mod.triage("...", client=_Client(fenced))
    assert res.intent == "PLASMID_CLONING"


def test_triage_unparseable_falls_back_to_reject():
    res = triage_mod.triage("...", client=_Client("not json at all"))
    assert res.intent == "REJECT"
    assert "parse failure" in (res.rejection_reason or "")


def test_triage_unknown_intent_falls_back_to_reject():
    bad = json.dumps({"intent": "WEIRD", "shorthand": "x", "is_new_topic": True})
    res = triage_mod.triage("...", client=_Client(bad))
    assert res.intent == "REJECT"
    assert "unknown intent" in (res.rejection_reason or "")


def test_triage_passes_user_message_and_attachments_flag():
    client = _Client(_payload("PLASMID_CLONING", "test"))
    triage_mod.triage("clone gene X", has_attachments=True, client=client)
    call = client.messages.calls[0]
    assert call["model"] == triage_mod.os.getenv("AGENT_MODEL", "claude-sonnet-4-6")
    user_content = call["messages"][0]["content"]
    assert "clone gene X" in user_content
    assert "Has attachments: True" in user_content
