"""Interpreter agent tests — mock both the Anthropic client and the
dispatch function so we can drive the ReAct loop deterministically."""
from __future__ import annotations

import asyncio
import json
from types import SimpleNamespace

import agent_v2  # noqa: F401
from agent_v2.interpreter.agent import run_interpreter
from agent_v2.interpreter.plasmid_registry import PlasmidRegistry


class _MockClient:
    """Returns a queue of canned responses. Each call to messages.create
    pops the next response."""
    def __init__(self, responses):
        self._responses = list(responses)
        self.calls = []
        self.messages = SimpleNamespace(create=self._create)

    def _create(self, **kwargs):
        self.calls.append(kwargs)
        return self._responses.pop(0)


def _resp_text(text):
    return SimpleNamespace(
        content=[SimpleNamespace(type="text", text=text)],
        stop_reason="end_turn",
    )


def _resp_tool(tool_id, name, inp):
    return SimpleNamespace(
        content=[SimpleNamespace(type="tool_use", id=tool_id, name=name, input=inp)],
        stop_reason="tool_use",
    )


def _resp_text_then_end(text):
    return SimpleNamespace(
        content=[SimpleNamespace(type="text", text=text)],
        stop_reason="end_turn",
    )


def test_empty_registry_runs_agent_anyway():
    """Empty registry no longer short-circuits — the agent still runs
    so find_external_part can answer questions about a named plasmid
    even when nothing has been uploaded. Provide a mock client that
    returns end_turn immediately."""
    reg = PlasmidRegistry()
    client = _MockClient([_resp_text("Drop a .gb to inspect — or name a specific plasmid.")])
    res = asyncio.run(run_interpreter("hi", reg, client=client,
                                       dispatch_fn=lambda *a, **kw: {}))
    assert "drop" in res.answer.lower() or "plasmid" in res.answer.lower()
    assert res.n_tool_calls == 0


def test_tool_use_then_final_text():
    reg = PlasmidRegistry()
    reg.register("p1", {
        "sequence": "A" * 100,
        "annotations": [],
        "modules": [{"module_type": "guide_expression_cassette", "rule_id": "POL3-GG-01",
                      "name": "U6 cassette", "start": 10, "end": 200, "strand": 1,
                      "submodules": [], "golden_gate": {"enzyme": "BsmBI"}}],
        "interactions": [], "hierarchical_annotations": [], "cloning_features": [],
    }, name="lenti")

    captured_args: dict = {}

    def fake_dispatch(name, args, registry):
        captured_args["name"] = name
        captured_args["args"] = args
        return registry.fan_out("find_modules", query=args.get("query", ""))

    client = _MockClient([
        _resp_tool("t1", "find_modules", {"query": "guide_expression"}),
        _resp_text("The guide expression cassette is at 10–200 (POL3-GG-01, BsmBI)."),
    ])

    res = asyncio.run(run_interpreter(
        "Where is the gRNA cloning cassette?",
        reg, client=client, dispatch_fn=fake_dispatch,
    ))
    assert "10" in res.answer and "200" in res.answer
    assert captured_args["name"] == "find_modules"
    assert res.n_tool_calls == 1
    assert len(client.calls) == 2  # 2 Anthropic calls (tool, then final text)
    assert any(c["plasmid_id"] == "p1" for c in res.citations)


def test_no_results_recorded_in_trace():
    reg = PlasmidRegistry()
    reg.register("p1", {
        "sequence": "A" * 100, "annotations": [], "modules": [],
        "interactions": [], "hierarchical_annotations": [], "cloning_features": [],
    })

    def fake_dispatch(name, args, registry):
        return registry.fan_out("find_modules", query=args.get("query", ""))

    client = _MockClient([
        _resp_tool("t1", "find_modules", {"query": "nothing_here"}),
        _resp_text("None of the 1 uploaded plasmid contains anything matching that query."),
    ])
    res = asyncio.run(run_interpreter(
        "Where is the XYZ cassette?", reg, client=client, dispatch_fn=fake_dispatch,
    ))
    assert "none" in res.answer.lower()
    assert res.trace[0]["n_results"] == 0


def test_immediate_end_turn_uses_text_path():
    reg = PlasmidRegistry()
    reg.register("p1", {
        "sequence": "A" * 100, "annotations": [], "modules": [],
        "interactions": [], "hierarchical_annotations": [], "cloning_features": [],
    })
    client = _MockClient([_resp_text("Direct answer with no tool calls.")])
    res = asyncio.run(run_interpreter("trivial", reg, client=client,
                                       dispatch_fn=lambda *a, **kw: {}))
    assert "Direct answer" in res.answer
    assert res.n_tool_calls == 0
