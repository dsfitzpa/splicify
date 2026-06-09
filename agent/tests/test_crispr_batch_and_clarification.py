"""Tests for the iter 41 features:
  - Missing-info short-circuit in _run_crispr_pipeline.
  - Multi-target batch dispatch (3 targets -> 3 parallel design_guides calls).
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
from agent_v2 import crispr_clarification
from agent_v2.crispr_pipeline import _run_crispr_pipeline
from splicify_api.agent.agent_tools import AttachmentRegistry


# ---------------------------------------------------------------------------
# Test doubles (mirror test_crispr_pipeline + test_orchestrator)
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


def _explore_summary(summary, key_facts=None):
    return _Resp([_Block("text", text=json.dumps({
        "summary_md": summary, "key_facts": key_facts or {},
    }))])


def _text_resp(text):
    return _Resp([_Block("text", text=text)])


def _registry():
    reg = AttachmentRegistry()
    reg.register_product("test_plasmid", "ACGT" * 1000, circular=True)
    return reg


class _Deps:
    """Minimal injectable deps for the pipeline."""
    triage_client = None
    target_locator_client = None
    guide_strategist_client = None
    plan_client = None
    main_client = None
    summarizer_client = None
    dispatch_fn = _DispatchFake()
    tools = []
    skip_memory = True
    output_dir = None


# ---------------------------------------------------------------------------
# crispr_clarification.respond — standalone unit
# ---------------------------------------------------------------------------
def test_clarification_envelope_shape():
    env = crispr_clarification.respond(
        ["No editing target specified.", "pegRNA design needs alt sequence."],
        user_message="design a guide",
        shorthand="vague guide design",
    )
    assert env["ok"] is True
    assert env["files"] is None
    assert env["viz"] is None
    assert env["workflow"] == "crispr"
    assert env["n_tool_calls"] == 0
    assert env["missing_info"] == ["No editing target specified.", "pegRNA design needs alt sequence."]
    # Bullet list in reply.
    assert "- No editing target specified." in env["reply"]
    assert "- pegRNA design needs alt sequence." in env["reply"]
    # Echoes the user prompt for context.
    assert "vague guide design" in env["reply"]
    # Trace records one local tool call.
    assert env["agent_trace"][0]["tool"] == "crispr_clarification"


def test_clarification_empty_missing_info_uses_fallback():
    env = crispr_clarification.respond([], user_message="something")
    assert env["ok"] is True
    assert env["missing_info"] == ["The CRISPR-design request is missing information, but the agent could not identify what specifically is missing."]
    assert "could not identify" in env["reply"]


# ---------------------------------------------------------------------------
# Pipeline short-circuit when an Explore subagent flags missing_info
# ---------------------------------------------------------------------------
def test_pipeline_short_circuits_on_missing_target():
    """TargetLocator flags 'no target specified' -> pipeline returns
    clarification envelope without running Plan / Main / Summarizer."""
    target_c = _SeqClient([
        _explore_summary(
            "User did not name a residue or feature.",
            {
                "target_attachment_id": "att_product_1",
                "resolved_targets": [],
                "annotated_attachments": ["att_product_1"],
                "missing_info": ["No editing target specified - tell me which residue / base / feature / gene region you want to cut."],
            },
        ),
    ])
    strat_c = _SeqClient([
        _explore_summary("Cas9 NGG.", {"strategy": "sgRNA", "n_targets": 0, "missing_info": []}),
    ])
    # No plan/main/summarizer clients should be invoked — the short-circuit
    # must skip them. Set them to None and the pipeline must NOT touch them.
    deps = _Deps()
    deps.target_locator_client = target_c
    deps.guide_strategist_client = strat_c
    deps.plan_client = None
    deps.main_client = None
    deps.summarizer_client = None

    envelope = asyncio.run(_run_crispr_pipeline(
        user_message="design a guide for me",
        registry=_registry(), session_id="sess_x", state=None, deps=deps,
    ))
    assert envelope["ok"] is True
    assert envelope["workflow"] == "crispr"
    assert envelope["files"] is None
    assert "No editing target specified" in envelope["reply"]
    assert envelope["missing_info"] == ["No editing target specified - tell me which residue / base / feature / gene region you want to cut."]
    # Findings still surface so the frontend can show what the agent learned.
    assert {f["role"] for f in envelope["findings"]} == {"target_locator", "guide_strategist"}


def test_pipeline_short_circuits_on_missing_pegrna_edit():
    """GuideStrategist sets strategy=pegRNA but flags missing_info for the alt."""
    target_c = _SeqClient([_explore_summary("Resolved.", {
        "resolved_targets": [{"feature_name": "KEAP1", "plasmid_position": 1234}],
        "missing_info": [],
    })])
    strat_c = _SeqClient([_explore_summary(
        "pegRNA strategy chosen, but no edit alt given.",
        {"strategy": "pegRNA", "n_targets": 1,
         "missing_info": ["pegRNA design needs the desired edit - tell me the new amino acid (e.g. D10A) or the substitution / insertion / deletion at the target base."]},
    )])
    deps = _Deps()
    deps.target_locator_client = target_c
    deps.guide_strategist_client = strat_c

    envelope = asyncio.run(_run_crispr_pipeline(
        user_message="make a precise edit at residue 33 of KEAP1",
        registry=_registry(), session_id="sess_y", state=None, deps=deps,
    ))
    assert envelope["ok"] is True
    assert "pegRNA design needs the desired edit" in envelope["reply"]
    assert envelope["missing_info"] == [
        "pegRNA design needs the desired edit - tell me the new amino acid (e.g. D10A) or the substitution / insertion / deletion at the target base."
    ]


def test_pipeline_short_circuits_deduplicates_missing_info():
    """Both subagents flag the same missing piece -> reply lists it once."""
    target_c = _SeqClient([_explore_summary("ok", {
        "missing_info": ["No editing target specified."],
    })])
    strat_c = _SeqClient([_explore_summary("ok", {
        "missing_info": ["No editing target specified.", "Pick a nuclease - SpCas9 or AsCas12a."],
    })])
    deps = _Deps()
    deps.target_locator_client = target_c
    deps.guide_strategist_client = strat_c

    envelope = asyncio.run(_run_crispr_pipeline(
        user_message="design something",
        registry=_registry(), session_id="sess_dup", state=None, deps=deps,
    ))
    # Both items present, but the duplicate is collapsed.
    assert envelope["missing_info"] == [
        "No editing target specified.",
        "Pick a nuclease - SpCas9 or AsCas12a.",
    ]


def test_pipeline_proceeds_when_missing_info_is_empty():
    """Pipeline does NOT short-circuit when both subagents say everything is fine."""
    target_c = _SeqClient([_explore_summary("Resolved 1 target.", {
        "resolved_targets": [{"plasmid_position": 100}],
        "missing_info": [],
    })])
    strat_c = _SeqClient([_explore_summary("Cas9 NGG.", {
        "strategy": "sgRNA", "n_targets": 1, "missing_info": [],
    })])
    plan_c = _SeqClient([_text_resp("## Plan\n- [ ] 1. design_guides")])
    main_c = _SeqClient([_text_resp("designed.")])
    summarizer_c = _SeqClient([_text_resp("**reply**")])

    deps = _Deps()
    deps.target_locator_client = target_c
    deps.guide_strategist_client = strat_c
    deps.plan_client = plan_c
    deps.main_client = main_c
    deps.summarizer_client = summarizer_c

    envelope = asyncio.run(_run_crispr_pipeline(
        user_message="design a Cas9 sgRNA targeting residue 10 of GeneX",
        registry=_registry(), session_id="sess_ok", state=None, deps=deps,
    ))
    # Full pipeline ran — plan + main + summarizer all consumed their mocked responses.
    assert envelope["ok"] is True
    assert envelope["reply"] == "**reply**"
    assert "main_agent_draft" in envelope
    # missing_info NOT present on a normal envelope.
    assert "missing_info" not in envelope


# ---------------------------------------------------------------------------
# Multi-target batch — N=3 targets fan out in one Main agent response
# ---------------------------------------------------------------------------
def test_multi_target_main_agent_fans_out_three_design_calls_in_parallel():
    """When the Main agent returns 3 design_guides tool_use blocks in one
    response, the harness must dispatch all 3 in parallel via asyncio.gather
    (not sequentially)."""
    target_c = _SeqClient([_explore_summary("Resolved 3 targets.", {
        "resolved_targets": [
            {"feature_name": "KEAP1", "kind": "aa_residue", "offset": 10, "plasmid_position": 100},
            {"feature_name": "KEAP1", "kind": "aa_residue", "offset": 33, "plasmid_position": 200},
            {"feature_name": "KEAP1", "kind": "aa_residue", "offset": 88, "plasmid_position": 300},
        ],
        "missing_info": [],
    })])
    strat_c = _SeqClient([_explore_summary("Cas9 NGG, 3 targets.", {
        "strategy": "sgRNA", "n_targets": 3, "missing_info": [],
    })])
    plan_c = _SeqClient([_text_resp(
        "## Plan\n- [ ] 1. design_guides batched across 3 targets\n- [ ] 2. emit_guides_csv"
    )])

    # Main agent returns 3 tool_use blocks for design_guides in ONE response,
    # then a final text response (no more tools). The harness dispatches all
    # 3 in parallel via asyncio.gather.
    tool_use_resp = _Resp(
        [
            _Block("tool_use", id="t1", name="design_guides",
                    input={"attachment_id": "att_product_1", "region_start": 70,  "region_end": 130}),
            _Block("tool_use", id="t2", name="design_guides",
                    input={"attachment_id": "att_product_1", "region_start": 170, "region_end": 230}),
            _Block("tool_use", id="t3", name="design_guides",
                    input={"attachment_id": "att_product_1", "region_start": 270, "region_end": 330}),
        ],
        stop_reason="tool_use",
    )
    final_resp = _Resp([_Block("text", text="done")])
    main_c = _SeqClient([tool_use_resp, final_resp])
    summarizer_c = _SeqClient([_text_resp("**3-target reply**")])

    dispatch = _DispatchFake(results={"design_guides": {"ok": True, "guides": [], "summary": {}}})
    deps = _Deps()
    deps.target_locator_client = target_c
    deps.guide_strategist_client = strat_c
    deps.plan_client = plan_c
    deps.main_client = main_c
    deps.summarizer_client = summarizer_c
    deps.dispatch_fn = dispatch

    envelope = asyncio.run(_run_crispr_pipeline(
        user_message="design Cas9 guides against residues 10, 33, and 88 of KEAP1",
        registry=_registry(), session_id="sess_batch", state=None, deps=deps,
    ))
    assert envelope["ok"] is True
    # All 3 design_guides calls dispatched.
    design_calls = [c for c in dispatch.calls if c["name"] == "design_guides"]
    assert len(design_calls) == 3
    # And they came with the 3 distinct regions.
    regions = sorted((c["args"]["region_start"], c["args"]["region_end"]) for c in design_calls)
    assert regions == [(70, 130), (170, 230), (270, 330)]
    # Reply came from the summarizer.
    assert envelope["reply"] == "**3-target reply**"
