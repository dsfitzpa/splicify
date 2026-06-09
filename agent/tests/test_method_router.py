"""Tests for MethodRouter - Anthropic + dispatch_tool both mocked."""
import asyncio
import json
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.explore import ExploreFinding
from agent_v2.subagents import method_router
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


def _final(summary_md, key_facts):
    text = json.dumps({"summary_md": summary_md, "key_facts": key_facts})
    return _Resp([_Block("text", text=text)], stop_reason="end_turn")


def test_method_router_immediate_summary():
    reg = AttachmentRegistry()
    reg.register_product("target", "ACGT" * 200)
    client = _SeqClient([_final(
        "Single uploaded plasmid - no method to choose.",
        {"recommended_method": "synthesis_fallback", "recommended_score": 0.5,
         "runner_up": None, "feasibility_summary": []},
    )])
    finding = asyncio.run(method_router.run_method_router(
        "what should I do?", reg, client=client, dispatch_fn=_DispatchFake({}),
    ))
    assert isinstance(finding, ExploreFinding)
    assert finding.role == "method_router"
    assert finding.key_facts["recommended_method"] == "synthesis_fallback"
    assert len(client.messages.calls[0]["tools"]) == 3


def test_method_router_route_workflow_then_summary():
    reg = AttachmentRegistry()
    reg.register_product("target", "ACGT" * 200)
    reg.register_product("inv1", "GGCC" * 100)

    route_use = _Block("tool_use", id="tu_route", name="route_workflow",
                       input={"target_attachment_id": "att_product_1",
                              "inventory_attachment_ids": ["att_product_2"]})
    client = _SeqClient([
        _Resp([route_use], stop_reason="tool_use"),
        _final(
            "Gibson wins (score 0.85), Golden Gate runner-up (0.62).",
            {"recommended_method": "gibson", "recommended_score": 0.85,
             "runner_up": "golden_gate",
             "feasibility_summary": [
                 {"method": "gibson", "feasible": True, "score": 0.85},
                 {"method": "golden_gate", "feasible": True, "score": 0.62},
             ]},
        ),
    ])
    dispatch = _DispatchFake({
        "route_workflow": {
            "winner": "gibson",
            "reports": [
                {"method": "gibson", "feasible": True, "score": 0.85},
                {"method": "golden_gate", "feasible": True, "score": 0.62},
            ],
        }
    })
    finding = asyncio.run(method_router.run_method_router(
        "best way to clone gene X into pUC19?", reg, client=client, dispatch_fn=dispatch,
    ))
    assert finding.key_facts["recommended_method"] == "gibson"
    assert finding.key_facts["runner_up"] == "golden_gate"
    assert dispatch.calls[0]["name"] == "route_workflow"
    assert finding.trace[0]["tool"] == "route_workflow"
    assert "winner" in finding.trace[0]["result_keys"]


def test_method_router_digest_path_for_restriction_check():
    reg = AttachmentRegistry()
    reg.register_product("target", "ACGT" * 200)

    digest_use = _Block("tool_use", id="tu_d", name="digest_plasmid",
                        input={"attachment_id": "att_product_1",
                               "enzymes": ["EcoRI", "BamHI"]})
    client = _SeqClient([
        _Resp([digest_use], stop_reason="tool_use"),
        _final(
            "EcoRI/BamHI cuts twice each - restriction not unique-cutter feasible.",
            {"recommended_method": "gibson", "recommended_score": 0.7,
             "runner_up": "restriction",
             "feasibility_summary": [
                 {"method": "restriction", "feasible": False, "score": 0.3},
             ]},
        ),
    ])
    dispatch = _DispatchFake({
        "digest_plasmid": {"cut_positions": {"EcoRI": [120, 950], "BamHI": [430, 1200]}},
    })
    finding = asyncio.run(method_router.run_method_router(
        "EcoRI/BamHI cloning?", reg, client=client, dispatch_fn=dispatch,
    ))
    assert finding.key_facts["recommended_method"] == "gibson"
    assert dispatch.calls[0]["name"] == "digest_plasmid"


def test_method_router_max_iters_bailout():
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 100)
    route_use = _Block("tool_use", id="tu", name="route_workflow",
                       input={"target_attachment_id": "att_product_1"})
    client = _SeqClient([_Resp([route_use], stop_reason="tool_use")] * 9)
    dispatch = _DispatchFake({"route_workflow": {"winner": "gibson", "reports": []}})
    finding = asyncio.run(method_router.run_method_router(
        "loop", reg, client=client, dispatch_fn=dispatch, max_iters=2,
    ))
    assert "max iterations" in finding.summary_md.lower()
    assert len(finding.trace) == 2


def test_method_router_uses_only_three_tools():
    assert len(method_router.METHOD_ROUTER_TOOLS) == 3
    names = {t["name"] for t in method_router.METHOD_ROUTER_TOOLS}
    assert names == {"route_workflow", "verify_assembly", "digest_plasmid"}
