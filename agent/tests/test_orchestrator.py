"""Tests for the orchestrator — every Anthropic client mocked."""
import asyncio
import json
import os
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
os.environ.setdefault("REDIS_URL", "redis://localhost:6379/1")

import pytest

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2 import memory
from agent_v2.orchestrator import OrchestratorDeps, run_orchestrator
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
    def __init__(self, results=None):
        self.results = results or {}
        self.calls = []

    async def __call__(self, name, args, registry):
        self.calls.append({"name": name, "args": args})
        return self.results.get(name, {"ok": True})


def _triage_payload(intent, shorthand="(test)", reason=None):
    return _Resp([_Block("text", text=json.dumps({
        "intent": intent, "shorthand": shorthand,
        "is_new_topic": True, "rejection_reason": reason,
    }))])


def _explore_summary(summary, key_facts=None):
    return _Resp([_Block("text", text=json.dumps({
        "summary_md": summary, "key_facts": key_facts or {},
    }))])


def _plain_text(text):
    return _Resp([_Block("text", text=text)])


@pytest.fixture
def fresh_session():
    sid = memory.mint_session_id()
    yield sid
    memory.reset(sid)


def test_orchestrator_reject_path(fresh_session):
    deps = OrchestratorDeps(
        triage_client=_SeqClient([_triage_payload(
            "REJECT", shorthand="favorite color question",
            reason="outside molecular biology scope",
        )]),
    )
    result = asyncio.run(run_orchestrator(
        "what's your favorite color?", AttachmentRegistry(),
        session_id=fresh_session, deps=deps,
    ))
    assert result["intent"] == "REJECT"
    assert "favorite color question" in result["reply"]
    assert "claude.ai" in result["reply"].lower()
    assert result["files"] is None
    assert result["session_id"] == fresh_session
    state = memory.load(fresh_session)
    assert state is not None
    assert state.last_user_message == "what's your favorite color?"


def test_orchestrator_crispr_path(fresh_session):
    """CRISPR_GUIDE intent routes through _run_crispr_pipeline (was stub_crispr pre-iter-37).

    Mocks all five Anthropic clients touched by the slim pipeline:
    triage, target_locator, guide_strategist, plan, main, summarizer.
    """
    reg = AttachmentRegistry()
    reg.register_product("test", "ACGT" * 200)

    triage_c = _SeqClient([_triage_payload("CRISPR_GUIDE", shorthand="guide for EMX1")])
    target_locator_c = _SeqClient([
        _explore_summary("Resolved E1 of EMX1.", {"resolved_targets": [{"plasmid_position": 100}]}),
    ])
    guide_strategist_c = _SeqClient([
        _explore_summary("Cas9 NGG / doench2014.",
                          {"strategy": "sgRNA", "n_targets": 1, "pam": "NGG"}),
    ])
    plan_c = _SeqClient([_plain_text("## Plan\n- [ ] 1. design_guides on att_product_1")])
    main_c = _SeqClient([_plain_text("Designed 5 candidate sgRNAs; top score 78.")])
    summarizer_c = _SeqClient([_plain_text("**Cas9 sgRNA design** ... top guide spacer ATGCATGC...")])

    deps = OrchestratorDeps(
        triage_client=triage_c,
        target_locator_client=target_locator_c,
        guide_strategist_client=guide_strategist_c,
        plan_client=plan_c,
        main_client=main_c,
        summarizer_client=summarizer_c,
        dispatch_fn=_DispatchFake(),
        tools=[],
    )
    result = asyncio.run(run_orchestrator(
        "design a CRISPR guide for EMX1", reg,
        session_id=fresh_session, deps=deps,
    ))
    assert result["intent"] == "CRISPR_GUIDE"
    assert result["ok"] is True
    assert result["workflow"] == "crispr"
    assert "main_agent_draft" in result
    assert len(result["findings"]) == 2
    finding_roles = {f["role"] for f in result["findings"]}
    assert finding_roles == {"target_locator", "guide_strategist"}


def test_orchestrator_plasmid_full_flow(fresh_session):
    reg = AttachmentRegistry()
    reg.register_product("test", "ACGT" * 200)

    triage_c = _SeqClient([_triage_payload("PLASMID_CLONING", shorthand="Gibson into pUC19")])
    part_c = _SeqClient([_explore_summary("Found GFP.", {"resolved_parts": ["GFP"]})])
    target_c = _SeqClient([_explore_summary("Gibson succeeded.",
                                             {"assembly_method": "gibson",
                                              "verifier_passed": True})])
    method_c = _SeqClient([_explore_summary("Gibson 0.85 wins.",
                                             {"recommended_method": "gibson",
                                              "recommended_score": 0.85})])
    plan_c = _SeqClient([_plain_text(
        "## Plan\n- [ ] 1. annotate_attachment(att_product_1)\n"
        "- [ ] 2. emit_assembled_gb(att_product_1)\n"
    )])
    main_c = _SeqClient([_plain_text(
        "Built a Gibson assembly. See assembled.gb."
    )])

    deps = OrchestratorDeps(
        triage_client=triage_c,
        part_scout_client=part_c,
        target_builder_client=target_c,
        method_router_client=method_c,
        plan_client=plan_c,
        main_client=main_c,
        dispatch_fn=_DispatchFake(),
        tools=[],
    )
    result = asyncio.run(run_orchestrator(
        "Build me a Gibson into pUC19", reg,
        session_id=fresh_session, deps=deps,
    ))
    assert result["intent"] == "PLASMID_CLONING"
    assert "Gibson" in result["reply"]
    assert result["error"] is None
    assert "## Plan" in result["plan_md"]
    assert len(result["findings"]) == 3
    roles = [f["role"] for f in result["findings"]]
    assert roles == ["part_scout", "target_builder", "method_router"]
    # The shorthand made it onto the envelope
    assert result["shorthand"] == "Gibson into pUC19"


def test_orchestrator_mints_session_id_when_missing():
    deps = OrchestratorDeps(
        triage_client=_SeqClient([_triage_payload("REJECT", reason="x")]),
        skip_memory=True,
    )
    result = asyncio.run(run_orchestrator(
        "off-topic", AttachmentRegistry(), deps=deps,
    ))
    assert result["session_id"].startswith("sess_")


def test_orchestrator_skip_memory_does_not_persist(fresh_session):
    deps = OrchestratorDeps(
        triage_client=_SeqClient([_triage_payload("REJECT", reason="x")]),
        skip_memory=True,
    )
    asyncio.run(run_orchestrator(
        "off-topic", AttachmentRegistry(),
        session_id=fresh_session, deps=deps,
    ))
    assert memory.load(fresh_session) is None


def test_orchestrator_writes_plan_md_when_plan_dir_set(tmp_path, fresh_session):
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 100)

    deps = OrchestratorDeps(
        triage_client=_SeqClient([_triage_payload("PLASMID_CLONING")]),
        part_scout_client=_SeqClient([_explore_summary("ok")]),
        target_builder_client=_SeqClient([_explore_summary("ok")]),
        method_router_client=_SeqClient([_explore_summary("ok")]),
        plan_client=_SeqClient([_plain_text("## Plan\n- [ ] 1. emit_workflow_trace(x)\n")]),
        main_client=_SeqClient([_plain_text("done")]),
        dispatch_fn=_DispatchFake(),
        tools=[],
        output_dir=str(tmp_path),
    )
    result = asyncio.run(run_orchestrator(
        "build", reg, session_id=fresh_session, deps=deps,
    ))
    # plan.md was written under the configured base
    written = list(tmp_path.glob(f"{fresh_session}/turn_*/plan.md"))
    assert len(written) == 1
    body = written[0].read_text()
    assert "## Plan" in body
