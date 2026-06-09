"""Tests for TargetBuilder - Anthropic + dispatch_tool both mocked."""
import asyncio
import json
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.explore import ExploreFinding
from agent_v2.subagents import target_builder
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


def test_target_builder_immediate_summary():
    reg = AttachmentRegistry()
    reg.register_product("backbone", "ACGT" * 200)
    client = _SeqClient([_final(
        "No assembly required - already a complete plasmid.",
        {"assembly_method": "none", "product_attachment_id": None,
         "verifier_passed": None, "verifier_warnings": []},
    )])
    finding = asyncio.run(target_builder.run_target_builder(
        "what is this?", reg, client=client, dispatch_fn=_DispatchFake({}),
    ))
    assert isinstance(finding, ExploreFinding)
    assert finding.role == "target_builder"
    assert finding.key_facts["assembly_method"] == "none"
    # Tools sent to Claude: should be exactly 4
    assert len(client.messages.calls[0]["tools"]) == 5


def test_target_builder_simulate_assembly_then_summarize():
    reg = AttachmentRegistry()
    reg.register_product("vec", "ACGT" * 200)
    reg.register_product("ins", "GGGG" * 100)

    sim_use = _Block("tool_use", id="tu_sim", name="simulate_assembly",
                     input={"instruction": "Gibson assembly",
                            "target_attachment_id": "att_product_1"})
    verify_use = _Block("tool_use", id="tu_ver", name="verify_assembly",
                        input={"attachment_id": "att_product_3"})
    client = _SeqClient([
        _Resp([sim_use], stop_reason="tool_use"),
        _Resp([verify_use], stop_reason="tool_use"),
        _final(
            "Gibson succeeded; verifier passed.",
            {"assembly_method": "gibson", "product_attachment_id": "att_product_3",
             "verifier_passed": True, "verifier_warnings": []},
        ),
    ])
    dispatch = _DispatchFake({
        "simulate_assembly": {"product_attachment_id": "att_product_3", "ok": True},
        "verify_assembly": {"passed": True, "warnings": []},
    })
    finding = asyncio.run(target_builder.run_target_builder(
        "Gibson assembly please", reg, client=client, dispatch_fn=dispatch,
    ))
    assert finding.key_facts["assembly_method"] == "gibson"
    assert finding.key_facts["verifier_passed"] is True
    assert [c["name"] for c in dispatch.calls] == ["simulate_assembly", "verify_assembly"]
    assert len(finding.trace) == 2


def test_target_builder_golden_gate_path():
    reg = AttachmentRegistry()
    reg.register_product("p1", "ACGT" * 100)
    reg.register_product("p2", "GGCC" * 100)

    gg_use = _Block("tool_use", id="tu_gg", name="golden_gate_assemble",
                    input={"attachment_ids": ["att_product_1", "att_product_2"], "enzyme": "BsaI"})
    client = _SeqClient([
        _Resp([gg_use], stop_reason="tool_use"),
        _final(
            "Golden Gate (BsaI) succeeded.",
            {"assembly_method": "golden_gate", "product_attachment_id": "att_product_3",
             "verifier_passed": None, "verifier_warnings": []},
        ),
    ])
    dispatch = _DispatchFake({
        "golden_gate_assemble": {"product_attachment_id": "att_product_3", "enzyme": "BsaI"},
    })
    finding = asyncio.run(target_builder.run_target_builder(
        "build the Golden Gate", reg, client=client, dispatch_fn=dispatch,
    ))
    assert finding.key_facts["assembly_method"] == "golden_gate"
    assert dispatch.calls[0]["name"] == "golden_gate_assemble"


def test_target_builder_verifier_fail_recorded():
    reg = AttachmentRegistry()
    reg.register_product("vec", "ACGT" * 100)

    verify_use = _Block("tool_use", id="tu_v", name="verify_assembly",
                        input={"attachment_id": "att_product_1"})
    client = _SeqClient([
        _Resp([verify_use], stop_reason="tool_use"),
        _final(
            "Verifier flagged 2 warnings.",
            {"assembly_method": "none", "product_attachment_id": "att_product_1",
             "verifier_passed": False,
             "verifier_warnings": ["GFP lacks promoter", "AmpR strand mismatch"]},
        ),
    ])
    dispatch = _DispatchFake({"verify_assembly": {"passed": False, "warnings": [{"feature": "GFP"}]}})
    finding = asyncio.run(target_builder.run_target_builder(
        "verify this", reg, client=client, dispatch_fn=dispatch,
    ))
    assert finding.key_facts["verifier_passed"] is False
    assert len(finding.key_facts["verifier_warnings"]) == 2


def test_target_builder_uses_only_five_tools():
    assert len(target_builder.TARGET_BUILDER_TOOLS) == 5
    names = {t["name"] for t in target_builder.TARGET_BUILDER_TOOLS}
    assert names == {"annotate_attachment", "simulate_assembly",
                     "golden_gate_assemble", "verify_assembly",
                     "resolve_feature_position"}
