"""Tests for output_dir wiring + envelope.files collection."""
import asyncio
import json
import os
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
os.environ.setdefault("REDIS_URL", "redis://localhost:6379/1")

import pytest

from agent_v2 import main_agent, memory
from agent_v2.orchestrator import OrchestratorDeps, run_orchestrator
from splicify_api.agent.agent_tools import AttachmentRegistry


# ---- shared mock infra (kept local to file) ----

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


# ---- main_agent: trace captures file envelope ----

def test_main_agent_trace_captures_file_envelope():
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 100)
    use = _Block("tool_use", id="tu", name="emit_assembled_gb",
                 input={"attachment_id": "att_product_1"})
    client = _SeqClient([
        _Resp([use], stop_reason="tool_use"),
        _Resp([_Block("text", text="done")], stop_reason="end_turn"),
    ])

    class _FakeDispatch:
        async def __call__(self, name, args, registry):
            return {"ok": True,
                    "file": {"fileName": "assembled.gb", "dataBase64": "AAA="},
                    "length_bp": 100}

    result = asyncio.run(main_agent.run_main_agent(
        "x", [], reg, plan_md="",
        client=client, dispatch_fn=_FakeDispatch(), tools=[],
    ))
    assert result.trace[0]["tool"] == "emit_assembled_gb"
    assert "file" in result.trace[0]
    assert result.trace[0]["file"]["fileName"] == "assembled.gb"


# ---- orchestrator: end-to-end output_dir wiring + envelope.files ----

def test_orchestrator_collects_emitter_files_and_writes_to_disk(tmp_path, fresh_session):
    """Real dispatch_with_emitters; main_agent calls emit_assembled_gb once."""
    reg = AttachmentRegistry()
    aid = reg.register_product("p", "ACGT" * 200)

    triage_c = _SeqClient([_triage_payload("PLASMID_CLONING")])
    part_c = _SeqClient([_explore_summary("ok")])
    target_c = _SeqClient([_explore_summary("ok")])
    method_c = _SeqClient([_explore_summary("ok")])
    plan_c = _SeqClient([_plain_text("## Plan\n- [ ] 1. emit_assembled_gb(att)\n")])
    emit_use = _Block("tool_use", id="tu_e", name="emit_assembled_gb",
                      input={"attachment_id": aid})
    main_c = _SeqClient([
        _Resp([emit_use], stop_reason="tool_use"),
        _Resp([_Block("text", text="Done.")], stop_reason="end_turn"),
    ])

    deps = OrchestratorDeps(
        triage_client=triage_c,
        part_scout_client=part_c,
        target_builder_client=target_c,
        method_router_client=method_c,
        plan_client=plan_c,
        main_client=main_c,
        # dispatch_fn=None: real dispatch_with_emitters used (with bound output_dir)
        # tools=None: full roster used
        output_dir=str(tmp_path),
    )

    result = asyncio.run(run_orchestrator(
        "build", reg, session_id=fresh_session, deps=deps,
    ))

    # On disk:
    written = list(tmp_path.glob(f"{fresh_session}/turn_*/assembled.gb"))
    assert len(written) == 1, f"expected 1 file, found: {list(tmp_path.rglob('*'))}"

    # Envelope:
    assert result["files"] is not None
    file_names = [f["fileName"] for f in result["files"]]
    assert "assembled.gb" in file_names
    # turn_id surfaces:
    assert result["turn_id"].startswith("turn_")


def test_orchestrator_no_output_dir_means_files_none(fresh_session):
    """When output_dir is unset, no files written and envelope.files stays None."""
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 200)

    deps = OrchestratorDeps(
        triage_client=_SeqClient([_triage_payload("PLASMID_CLONING")]),
        part_scout_client=_SeqClient([_explore_summary("ok")]),
        target_builder_client=_SeqClient([_explore_summary("ok")]),
        method_router_client=_SeqClient([_explore_summary("ok")]),
        plan_client=_SeqClient([_plain_text("## Plan\n- [ ] 1. (nothing)\n")]),
        main_client=_SeqClient([_plain_text("done immediately, no tools called")]),
        # No output_dir -> default dispatch_with_emitters with output_dir=None
        # No tool calls in main_agent -> no files possible
    )
    result = asyncio.run(run_orchestrator(
        "build", reg, session_id=fresh_session, deps=deps,
    ))
    # Server-side safety net: workflow_trace.txt is always emitted, AND
    # assembled.gb is auto-emitted when the Main agent skipped it and the
    # registry has any attachment to synthesise a construct from. The
    # fixture registers product 'p' above, so both files appear here.
    assert isinstance(result["files"], list)
    file_names = sorted(f["fileName"] for f in result["files"])
    assert file_names == ["assembled.gb", "workflow_trace.txt"]


def test_orchestrator_default_tools_is_full_roster(fresh_session, tmp_path):
    """When deps.tools is None, the full roster (v1 + 4 emitters) reaches Claude."""
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 100)
    main_c = _SeqClient([_plain_text("done")])

    deps = OrchestratorDeps(
        triage_client=_SeqClient([_triage_payload("PLASMID_CLONING")]),
        part_scout_client=_SeqClient([_explore_summary("ok")]),
        target_builder_client=_SeqClient([_explore_summary("ok")]),
        method_router_client=_SeqClient([_explore_summary("ok")]),
        plan_client=_SeqClient([_plain_text("## Plan\n- [ ] 1.\n")]),
        main_client=main_c,
        output_dir=str(tmp_path),
    )
    asyncio.run(run_orchestrator(
        "build", reg, session_id=fresh_session, deps=deps,
    ))
    sent_tools = main_c.messages.calls[0]["tools"]
    names = {t["name"] for t in sent_tools}
    # at least all 4 emitters + a sample of v1 tools
    # emit_workflow_trace is NOT in the LLM roster — it's auto-emitted
    # server-side now.
    for n in ("emit_assembled_gb", "emit_parts_order", "emit_protocol",
              "annotate_attachment", "verify_assembly"):
        assert n in names
    assert "emit_workflow_trace" not in names


# ---- summarizer wiring ----

def test_orchestrator_uses_summarizer_to_polish_reply(fresh_session):
    reg = AttachmentRegistry()
    reg.register_product("p", "ACGT" * 100)

    deps = OrchestratorDeps(
        triage_client=_SeqClient([_triage_payload("PLASMID_CLONING")]),
        part_scout_client=_SeqClient([_explore_summary("ok")]),
        target_builder_client=_SeqClient([_explore_summary("ok")]),
        method_router_client=_SeqClient([_explore_summary("ok")]),
        plan_client=_SeqClient([_plain_text("## Plan\n- [ ] 1.\n")]),
        main_client=_SeqClient([_plain_text("verbose draft from main")]),
        summarizer_client=_SeqClient([_plain_text("**Polished by summarizer.**")]),
        skip_memory=True,
    )
    result = asyncio.run(run_orchestrator(
        "build", reg, session_id=fresh_session, deps=deps,
    ))
    assert result["reply"] == "**Polished by summarizer.**"
    # The original main-agent draft is preserved separately:
    assert result["main_agent_draft"] == "verbose draft from main"
