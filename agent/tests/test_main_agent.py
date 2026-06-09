"""Tests for the Main Agent — Anthropic + dispatch_tool both mocked."""
import asyncio
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.explore import ExploreFinding
from agent_v2 import main_agent
from agent_v2.main_agent import run_main_agent, MainAgentResult, _mark_plan_item_done
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


def _findings():
    return [
        ExploreFinding(role="part_scout", summary_md="Found GFP."),
        ExploreFinding(role="target_builder", summary_md="Gibson succeeded."),
        ExploreFinding(role="method_router", summary_md="Gibson 0.85 wins."),
    ]


SAMPLE_PLAN = """## Plan
- [ ] 1. annotate_attachment(att_product_1)
- [ ] 2. simulate_assembly(instruction="Gibson")
- [ ] 3. verify_assembly(att_product_2)
- [ ] 4. emit_assembled_gb(att_product_2)
- [ ] 5. emit_parts_order(att_product_2)
- [ ] 6. emit_protocol(att_product_2)
- [ ] 7. emit_workflow_trace(att_product_2)
"""


def test_mark_plan_item_done_flips_first_match():
    md = "- [ ] annotate_attachment(att_1)\n- [ ] verify_assembly(att_2)\n"
    out, ok = _mark_plan_item_done(md, "annotate_attachment")
    assert ok is True
    assert "- [x] annotate_attachment(att_1)" in out
    assert "- [ ] verify_assembly(att_2)" in out


def test_mark_plan_item_done_no_match():
    md = "- [ ] something_else\n"
    out, ok = _mark_plan_item_done(md, "annotate_attachment")
    assert ok is False
    assert out == md


def test_main_agent_immediate_final_text():
    reg = AttachmentRegistry()
    client = _SeqClient([_Resp([_Block("text", text="Done. See attached files.")])])
    result = asyncio.run(run_main_agent(
        "build a Gibson", _findings(), reg,
        client=client, dispatch_fn=_DispatchFake({}),
        tools=[],  # don't pull v1 schemas in this unit test
    ))
    assert isinstance(result, MainAgentResult)
    assert result.final_text == "Done. See attached files."
    assert result.n_tool_calls == 0
    # default model used:
    assert client.messages.calls[0]["model"] == "claude-sonnet-4-6"


def test_main_agent_tool_call_then_final_text_marks_plan():
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 100)
    tool_use = _Block("tool_use", id="tu_1", name="annotate_attachment",
                      input={"attachment_id": "att_product_1"})
    client = _SeqClient([
        _Resp([tool_use], stop_reason="tool_use"),
        _Resp([_Block("text", text="All steps done.")], stop_reason="end_turn"),
    ])
    dispatch = _DispatchFake({"annotate_attachment": {"features": []}})

    result = asyncio.run(run_main_agent(
        "annotate this", _findings(), reg,
        plan_md=SAMPLE_PLAN, client=client, dispatch_fn=dispatch, tools=[],
    ))

    assert result.n_tool_calls == 1
    assert result.trace[0]["tool"] == "annotate_attachment"
    # plan crossing-off:
    assert "- [x] 1. annotate_attachment" in result.plan_md_final
    # other items remain unchecked:
    assert "- [ ] 3. verify_assembly" in result.plan_md_final


def test_main_agent_writes_plan_path(tmp_path):
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 100)
    tool_use = _Block("tool_use", id="tu_1", name="verify_assembly",
                      input={"attachment_id": "att_product_1"})
    client = _SeqClient([
        _Resp([tool_use], stop_reason="tool_use"),
        _Resp([_Block("text", text="Verified.")], stop_reason="end_turn"),
    ])
    dispatch = _DispatchFake({"verify_assembly": {"passed": True}})
    plan_path = tmp_path / "sess_abc" / "turn_1" / "plan.md"

    asyncio.run(run_main_agent(
        "verify", _findings(), reg,
        plan_md=SAMPLE_PLAN, plan_path=plan_path,
        client=client, dispatch_fn=dispatch, tools=[],
    ))

    assert plan_path.exists()
    body = plan_path.read_text()
    assert "- [x] 3. verify_assembly" in body


def test_main_agent_max_iters_bailout():
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 100)
    looper = _Block("tool_use", id="tu_loop", name="annotate_attachment",
                    input={"attachment_id": "att_product_1"})
    client = _SeqClient([_Resp([looper], stop_reason="tool_use")] * 9)
    dispatch = _DispatchFake({"annotate_attachment": {"features": []}})
    result = asyncio.run(run_main_agent(
        "loop", _findings(), reg,
        plan_md=SAMPLE_PLAN, client=client, dispatch_fn=dispatch,
        tools=[], max_iters=3,
    ))
    assert "max iterations" in result.final_text.lower()
    assert result.n_tool_calls == 3


def test_main_agent_default_tools_full_roster_when_unset():
    reg = AttachmentRegistry()
    client = _SeqClient([_Resp([_Block("text", text="ok")])])
    asyncio.run(run_main_agent(
        "noop", _findings(), reg,
        client=client, dispatch_fn=_DispatchFake({}),
        # tools omitted -> uses v1's AIPLASMIDDESIGN_TOOLS
    ))
    sent_tools = client.messages.calls[0]["tools"]
    names = {t["name"] for t in sent_tools}
    # at least the 11 v1 tools should be present
    expected_subset = {
        "annotate_attachment", "simulate_assembly", "golden_gate_assemble",
        "digest_plasmid", "find_primer_binding_sites", "score_sanger_primer",
        "lookup_kb_part", "route_workflow", "verify_assembly",
        "compare_to_choice",
    }
    assert expected_subset.issubset(names)


def test_main_agent_user_message_includes_plan_and_findings():
    reg = AttachmentRegistry()
    client = _SeqClient([_Resp([_Block("text", text="ok")])])
    asyncio.run(run_main_agent(
        "build a Gibson", _findings(), reg,
        plan_md=SAMPLE_PLAN, client=client,
        dispatch_fn=_DispatchFake({}), tools=[],
    ))
    sent_raw = client.messages.calls[0]["messages"][0]["content"]
    sent = sent_raw[0]["text"] if isinstance(sent_raw, list) else sent_raw
    assert "build a Gibson" in sent
    assert "annotate_attachment(att_product_1)" in sent  # from plan
    assert "Found GFP" in sent  # part_scout finding
    assert "Gibson succeeded" in sent  # target_builder finding
    assert "Gibson 0.85 wins" in sent  # method_router finding
