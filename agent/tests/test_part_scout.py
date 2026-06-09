"""Tests for PartScout — Anthropic client + dispatch_tool both mocked."""
import asyncio
import json
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.explore import ExploreFinding
from agent_v2.subagents import part_scout
from splicify_api.agent.agent_tools import AttachmentRegistry


class _Block:
    def __init__(self, type_, **kwargs):
        self.type = type_
        for k, v in kwargs.items():
            setattr(self, k, v)


class _Resp:
    def __init__(self, content, stop_reason="end_turn"):
        self.content = content
        self.stop_reason = stop_reason


class _SeqMessages:
    def __init__(self, responses):
        self._responses = list(responses)
        self.calls = []

    def create(self, **kwargs):
        self.calls.append(kwargs)
        if not self._responses:
            raise RuntimeError("no more mocked responses")
        return self._responses.pop(0)


class _SeqClient:
    def __init__(self, responses):
        self.messages = _SeqMessages(responses)


class _DispatchFake:
    def __init__(self, results_by_name):
        self.results_by_name = results_by_name
        self.calls = []

    async def __call__(self, name, args, registry):
        self.calls.append({"name": name, "args": args})
        return self.results_by_name.get(name, {"error": f"no fake for {name}"})


def _final_summary_response(summary_md, key_facts):
    text = json.dumps({"summary_md": summary_md, "key_facts": key_facts})
    return _Resp([_Block("text", text=text)], stop_reason="end_turn")


def test_part_scout_immediate_summary():
    reg = AttachmentRegistry()
    reg.register_product("test_plasmid", "ACGT" * 200)
    client = _SeqClient([
        _final_summary_response(
            "Found pUC19 backbone + GFP CDS.",
            {"resolved_parts": ["GFP"], "kb_hits": 1,
             "annotated_attachments": ["att_product_1"]},
        ),
    ])
    dispatch = _DispatchFake({})
    finding = asyncio.run(part_scout.run_part_scout(
        "what's on this plasmid?", reg, client=client, dispatch_fn=dispatch,
    ))
    assert isinstance(finding, ExploreFinding)
    assert finding.role == "part_scout"
    assert "pUC19" in finding.summary_md
    assert finding.key_facts["resolved_parts"] == ["GFP"]
    assert finding.references == ["att_product_1"]
    assert dispatch.calls == []
    assert client.messages.calls[0]["model"] == "claude-sonnet-4-6"
    # Tools were sent to Claude
    assert len(client.messages.calls[0]["tools"]) == 3


def test_part_scout_runs_tool_then_summarizes():
    reg = AttachmentRegistry()
    reg.register_product("plasmid_a", "ACGT" * 200)

    tool_use = _Block("tool_use", id="tu_1", name="annotate_attachment",
                      input={"attachment_id": "att_product_1"})
    client = _SeqClient([
        _Resp([tool_use], stop_reason="tool_use"),
        _final_summary_response(
            "Annotated; found GFP, AmpR, ori.",
            {"resolved_parts": ["GFP", "AmpR", "ori"], "kb_hits": 0,
             "annotated_attachments": ["att_product_1"]},
        ),
    ])
    dispatch = _DispatchFake({
        "annotate_attachment": {"features": [{"name": "GFP"}, {"name": "AmpR"}, {"name": "ori"}]},
    })

    finding = asyncio.run(part_scout.run_part_scout(
        "annotate this plasmid", reg, client=client, dispatch_fn=dispatch,
    ))

    assert finding.role == "part_scout"
    assert "GFP" in finding.summary_md
    assert dispatch.calls[0]["name"] == "annotate_attachment"
    assert finding.trace[0]["tool"] == "annotate_attachment"
    assert "features" in finding.trace[0]["result_keys"]
    # Two Anthropic calls: tool_use, then end_turn
    assert len(client.messages.calls) == 2


def test_part_scout_max_iters():
    reg = AttachmentRegistry()
    reg.register_product("plasmid_x", "ACGT" * 100)

    tool_use = _Block("tool_use", id="tu_loop", name="annotate_attachment",
                      input={"attachment_id": "att_product_1"})
    responses = [_Resp([tool_use], stop_reason="tool_use") for _ in range(9)]
    client = _SeqClient(responses)
    dispatch = _DispatchFake({"annotate_attachment": {"features": []}})

    finding = asyncio.run(part_scout.run_part_scout(
        "loop forever", reg, client=client, dispatch_fn=dispatch, max_iters=3,
    ))
    assert "max iterations" in finding.summary_md.lower()
    assert len(finding.trace) == 3


def test_part_scout_falls_back_when_final_text_not_json():
    reg = AttachmentRegistry()
    client = _SeqClient([_Resp([_Block("text", text="just plain prose")])])
    finding = asyncio.run(part_scout.run_part_scout(
        "tell me", reg, client=client, dispatch_fn=_DispatchFake({}),
    ))
    assert finding.summary_md == "just plain prose"
    assert finding.key_facts == {}


def test_part_scout_uses_only_three_tools():
    assert len(part_scout.PART_SCOUT_TOOLS) == 3
    names = {t["name"] for t in part_scout.PART_SCOUT_TOOLS}
    assert names == {"annotate_attachment", "lookup_kb_part", "analyze_design_intent"}
