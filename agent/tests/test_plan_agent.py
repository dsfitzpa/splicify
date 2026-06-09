"""Tests for the Plan agent — Anthropic client mocked, file IO via tmp_path."""
import asyncio
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from agent_v2.explore import ExploreFinding
from agent_v2.subagents import plan_agent


class _Block:
    type = "text"

    def __init__(self, text):
        self.text = text


class _Resp:
    def __init__(self, text):
        self.content = [_Block(text)]


class _Messages:
    def __init__(self, text):
        self._text = text
        self.calls = []

    def create(self, **kwargs):
        self.calls.append(kwargs)
        return _Resp(self._text)


class _Client:
    def __init__(self, text):
        self.messages = _Messages(text)


SAMPLE_PLAN = """## Plan
- [ ] 1. annotate_attachment(att_product_1)
- [ ] 2. simulate_assembly(instruction="Gibson assembly", target_attachment_id="att_product_1")
- [ ] 3. verify_assembly(att_product_2)
- [ ] 4. emit_assembled_gb(att_product_2)
- [ ] 5. emit_parts_order(att_product_2)
- [ ] 6. emit_protocol(att_product_2)
- [ ] 7. emit_workflow_trace(att_product_2)
"""


def _findings():
    return [
        ExploreFinding(role="part_scout", summary_md="Found GFP, AmpR, ori.",
                       key_facts={"resolved_parts": ["GFP", "AmpR", "ori"]},
                       references=["att_product_1"]),
        ExploreFinding(role="target_builder", summary_md="Gibson succeeded.",
                       key_facts={"assembly_method": "gibson",
                                  "product_attachment_id": "att_product_2",
                                  "verifier_passed": True}),
        ExploreFinding(role="method_router", summary_md="Gibson 0.85 wins.",
                       key_facts={"recommended_method": "gibson", "recommended_score": 0.85}),
    ]


def test_plan_agent_returns_markdown_and_step_count():
    client = _Client(SAMPLE_PLAN)
    result = asyncio.run(plan_agent.run_plan_agent(
        "Build me a Gibson", _findings(), client=client,
    ))
    assert isinstance(result, plan_agent.PlanResult)
    assert "## Plan" in result.plan_md
    assert result.n_steps == 7
    assert result.plan_path is None


def test_plan_agent_writes_file_when_path_given(tmp_path):
    client = _Client(SAMPLE_PLAN)
    target = tmp_path / "sess_xyz" / "turn_1" / "plan.md"
    result = asyncio.run(plan_agent.run_plan_agent(
        "build", _findings(), client=client, plan_path=target,
    ))
    assert target.exists()
    assert result.plan_path == target
    assert "emit_assembled_gb" in target.read_text()


def test_plan_agent_empty_response_falls_back():
    client = _Client("")
    result = asyncio.run(plan_agent.run_plan_agent("...", _findings(), client=client))
    assert "no content" in result.plan_md.lower()
    assert result.n_steps == 1  # the placeholder line is itself a checklist item


def test_plan_agent_user_message_includes_findings():
    client = _Client(SAMPLE_PLAN)
    asyncio.run(plan_agent.run_plan_agent(
        "Build a Gibson", _findings(), client=client,
    ))
    sent = client.messages.calls[0]["messages"][0]["content"]
    assert "Build a Gibson" in sent
    assert "part_scout" in sent
    assert "target_builder" in sent
    assert "method_router" in sent
    assert "Found GFP" in sent
    assert "Gibson succeeded" in sent


def test_plan_agent_no_tools_sent_to_claude():
    client = _Client(SAMPLE_PLAN)
    asyncio.run(plan_agent.run_plan_agent("...", _findings(), client=client))
    call = client.messages.calls[0]
    assert "tools" not in call


def test_plan_path_for_helper(tmp_path):
    p = plan_agent.plan_path_for("sess_abc", "turn_3", base=str(tmp_path))
    assert p == tmp_path / "sess_abc" / "turn_3" / "plan.md"


def test_plan_agent_uses_sonnet_4_6_by_default():
    client = _Client(SAMPLE_PLAN)
    asyncio.run(plan_agent.run_plan_agent("...", _findings(), client=client))
    assert client.messages.calls[0]["model"] == "claude-sonnet-4-6"
