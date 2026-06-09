"""Tests for the CRISPR pipeline + its two new Explore subagents.

All Anthropic clients are mocked via _SeqClient (which mimics the
.messages.create interface). v1 tool dispatch is faked via _DispatchFake.
"""
import asyncio
import json
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
os.environ.setdefault("REDIS_URL", "redis://localhost:6379/1")

import pytest

import agent_v2  # noqa: F401
from agent_v2.crispr_pipeline import (
    CRISPR_MAIN_SYSTEM_PROMPT,
    CRISPR_PLAN_SYSTEM_PROMPT,
    _run_crispr_pipeline,
)
from agent_v2.orchestrator import OrchestratorDeps
from agent_v2.subagents.guide_strategist import run_guide_strategist
from agent_v2.subagents.target_locator import run_target_locator
from splicify_api.agent.agent_tools import AttachmentRegistry

def _system_text(sys_arg):
    """Normalise an Anthropic `system` payload to a plain string. The SDK
    accepts both legacy str and the prompt-caching list shape
    [{"type": "text", "text": ..., "cache_control": ...}]; tests should
    compare against the concatenated text either way."""
    if isinstance(sys_arg, str):
        return sys_arg
    if isinstance(sys_arg, list):
        return "".join(b.get("text", "") for b in sys_arg if isinstance(b, dict))
    return ""



# ---------------------------------------------------------------------------
# Test doubles (same shape as test_orchestrator's)
# ---------------------------------------------------------------------------
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


def _explore_summary_resp(summary, key_facts=None):
    return _Resp([_Block("text", text=json.dumps({
        "summary_md": summary, "key_facts": key_facts or {},
    }))])


def _text_resp(text):
    return _Resp([_Block("text", text=text)])


def _registry():
    reg = AttachmentRegistry()
    reg.register_product("test_plasmid", "ACGT" * 1000, circular=True)
    return reg


# ---------------------------------------------------------------------------
# TargetLocator subagent
# ---------------------------------------------------------------------------
def test_target_locator_immediate_summary():
    """No tool calls; subagent emits the JSON digest in one turn."""
    reg = _registry()
    client = _SeqClient([
        _explore_summary_resp(
            "Resolved D10A in Cas9 at plasmid position 5138.",
            {"resolved_targets": [{"feature_name": "Cas9", "plasmid_position": 5138,
                                    "codon": "GAT", "amino_acid": "D"}]},
        ),
    ])
    finding = asyncio.run(run_target_locator(
        "edit D10A in Cas9 on the uploaded plasmid", reg,
        client=client, dispatch_fn=_DispatchFake(),
    ))
    assert finding.role == "target_locator"
    assert "5138" in finding.summary_md
    assert finding.key_facts["resolved_targets"][0]["codon"] == "GAT"
    # System prompt must mention the two-tool roster.
    sent = client.messages.calls[0]
    assert "annotate_attachment" in _system_text(sent["system"])
    assert "resolve_feature_position" in _system_text(sent["system"])


# ---------------------------------------------------------------------------
# GuideStrategist subagent
# ---------------------------------------------------------------------------
def test_guide_strategist_immediate_summary():
    reg = _registry()
    client = _SeqClient([
        _explore_summary_resp(
            "Cas9 sgRNA, NGG, doench2014, n_targets=1, readouts: NGS + Sanger.",
            {"strategy": "sgRNA", "n_targets": 1, "nuclease": "SpCas9",
             "pam": "NGG", "pam_position": "3prime", "guide_length": 20,
             "scoring_method": "doench2014",
             "readouts": ["illumina_ngs", "sanger"]},
        ),
    ])
    finding = asyncio.run(run_guide_strategist(
        "use Cas9 to knock in a stop codon at residue 10 of KEAP1", reg,
        client=client, dispatch_fn=_DispatchFake(),
    ))
    assert finding.role == "guide_strategist"
    assert finding.key_facts["strategy"] == "sgRNA"
    assert finding.key_facts["pam"] == "NGG"
    sent = client.messages.calls[0]
    # System prompt enumerates the strategy decisions.
    sys_text = _system_text(sent["system"])
    assert "Cas9" in sys_text and "pegRNA" in sys_text


# ---------------------------------------------------------------------------
# Full _run_crispr_pipeline
# ---------------------------------------------------------------------------
def test_crispr_pipeline_end_to_end_no_files():
    """Pipeline runs through both subagents + plan + main + summarizer.

    Main agent returns final text with no tool_use blocks, so `files` is None.
    """
    reg = _registry()
    target_c = _SeqClient([
        _explore_summary_resp("Resolved.", {"resolved_targets": [{"plasmid_position": 1000}]}),
    ])
    strat_c = _SeqClient([
        _explore_summary_resp("Cas9 NGG.", {"strategy": "sgRNA", "n_targets": 1, "pam": "NGG"}),
    ])
    plan_c = _SeqClient([_text_resp("## Plan\n- [ ] 1. design_guides")])
    main_c = _SeqClient([_text_resp("Top sgRNA spacer ATGCATGCATGCATGCATGC, score 78.")])
    summarizer_c = _SeqClient([_text_resp("**Cas9 sgRNA** -- top guide ATGC...")])

    class _Deps:
        triage_client = None
        target_locator_client = target_c
        guide_strategist_client = strat_c
        plan_client = plan_c
        main_client = main_c
        summarizer_client = summarizer_c
        dispatch_fn = _DispatchFake()
        tools = []
        skip_memory = True
        output_dir = None

    envelope = asyncio.run(_run_crispr_pipeline(
        user_message="design a CRISPR guide for EMX1",
        registry=reg,
        session_id="sess_test",
        state=None,
        deps=_Deps(),
    ))
    assert envelope["ok"] is True
    assert envelope["workflow"] == "crispr"
    assert envelope["reply"].startswith("**Cas9 sgRNA**")
    assert envelope["main_agent_draft"].startswith("Top sgRNA")
    # auto_emit_workflow_trace now produces workflow_trace.txt server-side
    # even when the main agent's trace recorded no other emitter calls.
    assert isinstance(envelope["files"], list)
    assert len(envelope["files"]) == 1
    assert envelope["files"][0]["fileName"] == "workflow_trace.txt"
    assert len(envelope["findings"]) == 2
    roles = {f["role"] for f in envelope["findings"]}
    assert roles == {"target_locator", "guide_strategist"}


def test_crispr_pipeline_writes_plan_md_when_output_dir_set(tmp_path):
    reg = _registry()
    target_c = _SeqClient([_explore_summary_resp("Resolved.", {})])
    strat_c = _SeqClient([_explore_summary_resp("Cas9.", {"strategy": "sgRNA"})])
    plan_c = _SeqClient([_text_resp("## Plan\n- [ ] 1. design_guides\n- [ ] 2. design_primers")])
    main_c = _SeqClient([_text_resp("done")])
    summarizer_c = _SeqClient([_text_resp("Final reply")])

    class _Deps:
        triage_client = None
        target_locator_client = target_c
        guide_strategist_client = strat_c
        plan_client = plan_c
        main_client = main_c
        summarizer_client = summarizer_c
        dispatch_fn = _DispatchFake()
        tools = []
        skip_memory = True
        output_dir = str(tmp_path)

    envelope = asyncio.run(_run_crispr_pipeline(
        user_message="design a guide",
        registry=reg, session_id="sess_a", state=None, deps=_Deps(),
    ))
    turn_id = envelope["turn_id"]
    plan_file = tmp_path / "sess_a" / turn_id / "plan.md"
    assert plan_file.exists()
    body = plan_file.read_text()
    # Main agent's plan-crossing-off doesn't fire on text-only response, so
    # plan.md should still carry the original 2 unchecked items.
    assert body.count("- [ ]") == 2


def test_crispr_pipeline_sends_crispr_system_prompt_to_plan_agent():
    """The plan agent must receive CRISPR_PLAN_SYSTEM_PROMPT, not the plasmid one."""
    reg = _registry()
    target_c = _SeqClient([_explore_summary_resp("ok", {})])
    strat_c = _SeqClient([_explore_summary_resp("ok", {})])
    plan_c = _SeqClient([_text_resp("## Plan")])
    main_c = _SeqClient([_text_resp("done")])
    summarizer_c = _SeqClient([_text_resp("reply")])

    class _Deps:
        triage_client = None
        target_locator_client = target_c
        guide_strategist_client = strat_c
        plan_client = plan_c
        main_client = main_c
        summarizer_client = summarizer_c
        dispatch_fn = _DispatchFake()
        tools = []
        skip_memory = True
        output_dir = None

    asyncio.run(_run_crispr_pipeline(
        user_message="design a guide",
        registry=reg, session_id="sess_b", state=None, deps=_Deps(),
    ))
    plan_call = plan_c.messages.calls[0]
    assert _system_text(plan_call["system"]) == CRISPR_PLAN_SYSTEM_PROMPT
    # Plan agent must NOT see the plasmid emit_assembled_gb requirement.
    assert "assembly method from MethodRouter" not in _system_text(plan_call["system"])
    # Main agent must see the CRISPR main prompt.
    main_call = main_c.messages.calls[0]
    assert _system_text(main_call["system"]) == CRISPR_MAIN_SYSTEM_PROMPT


def test_crispr_pipeline_collects_files_from_main_trace(monkeypatch):
    """When the main agent's trace records a file envelope, the pipeline surfaces it."""
    reg = _registry()

    async def fake_main(*args, **kwargs):
        # Return a MainAgentResult-shaped object with a file in the trace.
        from agent_v2.main_agent import MainAgentResult
        return MainAgentResult(
            final_text="done",
            trace=[
                {"iteration": 0, "tool": "emit_parts_order", "args_summary": "...",
                 "result_keys": ["file", "n_parts"],
                 "file": {"fileName": "parts_order.csv", "dataBase64": "Zm9v"}},
            ],
            plan_md_final="## Plan",
            n_tool_calls=1,
        )

    import agent_v2.crispr_pipeline as cp
    monkeypatch.setattr(cp, "run_main_agent", fake_main)

    target_c = _SeqClient([_explore_summary_resp("ok", {})])
    strat_c = _SeqClient([_explore_summary_resp("ok", {})])
    plan_c = _SeqClient([_text_resp("## Plan")])
    summarizer_c = _SeqClient([_text_resp("reply")])

    class _Deps:
        triage_client = None
        target_locator_client = target_c
        guide_strategist_client = strat_c
        plan_client = plan_c
        main_client = None  # unused due to fake_main
        summarizer_client = summarizer_c
        dispatch_fn = _DispatchFake()
        tools = []
        skip_memory = True
        output_dir = None

    envelope = asyncio.run(_run_crispr_pipeline(
        user_message="design",
        registry=reg, session_id="sess_c", state=None, deps=_Deps(),
    ))
    assert envelope["files"] is not None
    # 2 files: the parts_order.csv recorded in the main agent's trace +
    # the auto-emitted workflow_trace.txt the orchestrator now produces.
    fileNames = {f["fileName"] for f in envelope["files"]}
    assert fileNames == {"parts_order.csv", "workflow_trace.txt"}


def test_orchestrator_deps_has_new_crispr_fields():
    """OrchestratorDeps gained target_locator_client + guide_strategist_client."""
    deps = OrchestratorDeps()
    assert deps.target_locator_client is None
    assert deps.guide_strategist_client is None
