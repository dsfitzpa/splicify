"""Tests for the Summarizer subagent — Anthropic mocked, fallback paths covered."""
import asyncio
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from agent_v2.explore import ExploreFinding
from agent_v2.subagents import summarizer
from agent_v2.subagents.summarizer import run_summarizer


class _Block:
    type = "text"

    def __init__(self, text):
        self.text = text


class _Resp:
    def __init__(self, text):
        self.content = [_Block(text)] if text is not None else []


class _Messages:
    def __init__(self, text=None, raise_exc=None):
        self._text = text
        self._raise = raise_exc
        self.calls = []

    def create(self, **kwargs):
        self.calls.append(kwargs)
        if self._raise is not None:
            raise self._raise
        return _Resp(self._text)


class _Client:
    def __init__(self, text=None, raise_exc=None):
        self.messages = _Messages(text=text, raise_exc=raise_exc)


def _findings():
    return [
        ExploreFinding(role="part_scout", summary_md="Found GFP, AmpR, ori."),
        ExploreFinding(role="target_builder", summary_md="Gibson assembly succeeded."),
        ExploreFinding(role="method_router", summary_md="Gibson 0.85 wins."),
    ]


def test_summarizer_returns_polished_reply():
    polished = ("**Gibson assembly** built.\n\n"
                "I prioritised Gibson over Golden Gate because score 0.85 vs 0.62. "
                "Files: `assembled.gb`, `parts_order.csv`, `protocol.csv`, `workflow_trace.txt`.")
    client = _Client(text=polished)
    result = asyncio.run(run_summarizer(
        user_message="Build a Gibson",
        main_reply="(verbose draft from main agent)",
        files=[{"fileName": "assembled.gb"}],
        findings=_findings(),
    ))
    # default client is anthropic.Anthropic — but we passed our own:
    result = asyncio.run(run_summarizer(
        "Build a Gibson",
        main_reply="(verbose draft from main agent)",
        files=[{"fileName": "assembled.gb"}],
        findings=_findings(),
        client=client,
    ))
    assert result.reply.startswith("**Gibson")
    assert result.used_fallback is False


def test_summarizer_empty_response_falls_back_to_draft():
    client = _Client(text="")
    result = asyncio.run(run_summarizer(
        "Build", main_reply="DRAFT_FALLBACK",
        files=None, findings=[], client=client,
    ))
    assert result.reply == "DRAFT_FALLBACK"
    assert result.used_fallback is True


def test_summarizer_anthropic_error_falls_back():
    client = _Client(raise_exc=RuntimeError("network down"))
    result = asyncio.run(run_summarizer(
        "Build", main_reply="DRAFT", files=None, findings=[], client=client,
    ))
    assert result.reply == "DRAFT"
    assert result.used_fallback is True


def test_summarizer_user_block_contains_inputs():
    client = _Client(text="ok")
    decisions = [{"choice": "Gibson", "alternative": "Golden Gate", "reason": "score 0.85"}]
    asyncio.run(run_summarizer(
        "user prompt here",
        main_reply="main draft here",
        files=[{"fileName": "assembled.gb"}, {"fileName": "parts_order.csv"}],
        findings=_findings(),
        decisions=decisions,
        client=client,
    ))
    sent = client.messages.calls[0]["messages"][0]["content"]
    assert "user prompt here" in sent
    assert "main draft here" in sent
    assert "assembled.gb, parts_order.csv" in sent
    assert "Gibson over Golden Gate" in sent
    assert "Found GFP" in sent  # part_scout finding


def test_summarizer_no_tools_sent_to_claude():
    client = _Client(text="ok")
    asyncio.run(run_summarizer(
        "x", main_reply="y", files=None, findings=[], client=client,
    ))
    assert "tools" not in client.messages.calls[0]


def test_summarizer_default_model_is_sonnet_4_6():
    client = _Client(text="ok")
    asyncio.run(run_summarizer(
        "x", main_reply="y", files=None, findings=[], client=client,
    ))
    assert client.messages.calls[0]["model"] == "claude-sonnet-4-6"


def test_summarizer_handles_none_files_and_decisions():
    client = _Client(text="ok reply")
    result = asyncio.run(run_summarizer(
        "x", main_reply="y", files=None, findings=[], decisions=None, client=client,
    ))
    assert result.reply == "ok reply"
    sent = client.messages.calls[0]["messages"][0]["content"]
    assert "(no files)" in sent
    assert "(none recorded)" in sent


def test_summarizer_empty_response_with_empty_draft():
    client = _Client(text="")
    result = asyncio.run(run_summarizer(
        "x", main_reply="", files=None, findings=[], client=client,
    ))
    assert result.reply == "(no reply produced)"
    assert result.used_fallback is True
