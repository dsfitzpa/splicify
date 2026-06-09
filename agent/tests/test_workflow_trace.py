"""Tests for emit_workflow_trace — section assembly + optional blocks."""
import asyncio
import base64
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from agent_v2.outputs.workflow_trace import emit_workflow_trace
from agent_v2 import tools as v2_tools
from splicify_api.agent.agent_tools import AttachmentRegistry


def _decode(result):
    return base64.b64decode(result["file"]["dataBase64"]).decode()


def test_emit_workflow_trace_minimal_args():
    result = asyncio.run(emit_workflow_trace({}, AttachmentRegistry()))
    assert result["ok"] is True
    text = _decode(result)
    assert "=== agent_v2 workflow trace ===" in text
    # required sections present even with empty input:
    assert "--- AGENT TRACE ---" in text
    assert "--- PLAN.MD ---" in text
    assert "--- DECISIONS LEDGER ---" in text
    # optional sections absent
    assert "--- DESIGN VERIFICATION ---" not in text
    assert "--- EXPLORE FINDINGS ---" not in text


def test_emit_workflow_trace_full_args(tmp_path):
    args = {
        "session_id": "sess_abc",
        "turn_id": "turn_1_1715000000",
        "assembly_method": "gibson",
        "product_attachment_id": "att_product_2",
        "timestamp": "2026-05-07T12:34:56Z",
        "agent_trace": [
            {"iteration": 0, "tool": "annotate_attachment",
             "args_summary": "attachment_id=att_product_1",
             "result_keys": ["features", "modules"]},
            {"iteration": 1, "tool": "simulate_assembly",
             "args_summary": "instruction=Gibson", "result_keys": ["ok"]},
        ],
        "plan_md": "## Plan\n- [x] 1. annotate_attachment(...)\n- [x] 2. simulate_assembly(...)\n",
        "decisions": [
            {"choice": "Gibson", "alternative": "Golden Gate",
             "reason": "Score 0.85 vs 0.62 + PCR-friendly"},
        ],
        "verifier": {
            "passed": True,
            "warnings": [],
        },
        "findings": [
            {"role": "part_scout", "summary_md": "Found GFP, AmpR.",
             "key_facts": {"resolved_parts": ["GFP", "AmpR"]}},
        ],
    }
    result = asyncio.run(emit_workflow_trace(args, AttachmentRegistry(),
                                              output_dir=str(tmp_path)))
    text = _decode(result)
    assert "session: sess_abc" in text
    assert "turn:    turn_1_1715000000" in text
    assert "method:  gibson" in text
    assert "product: att_product_2" in text
    assert "iter 0: annotate_attachment(attachment_id=att_product_1) -> features, modules" in text
    assert "iter 1: simulate_assembly" in text
    assert "- [x] 1. annotate_attachment" in text
    assert "Chose Gibson; runner-up Golden Gate. Reason: Score 0.85" in text
    assert "--- DESIGN VERIFICATION ---" in text
    assert "Passed." in text
    assert "--- EXPLORE FINDINGS ---" in text
    assert "[part_scout] Found GFP, AmpR." in text
    assert "key_facts:" in text

    # on-disk write
    out_file = tmp_path / "workflow_trace.txt"
    assert out_file.exists()
    assert "Chose Gibson" in out_file.read_text()


def test_emit_workflow_trace_verifier_failed_with_warnings():
    args = {
        "verifier": {
            "passed": False,
            "warnings": [
                {"feature_name": "GFP", "remediation": "pair_with_cds"},
                "Generic string warning",
            ],
        },
    }
    text = _decode(asyncio.run(emit_workflow_trace(args, AttachmentRegistry())))
    assert "FAILED." in text
    assert "GFP: pair_with_cds" in text
    assert "Generic string warning" in text


def test_emit_workflow_trace_empty_trace_falls_back():
    text = _decode(asyncio.run(emit_workflow_trace({"agent_trace": []},
                                                    AttachmentRegistry())))
    assert "(no tool calls)" in text


def test_emit_workflow_trace_empty_decisions_falls_back():
    text = _decode(asyncio.run(emit_workflow_trace({}, AttachmentRegistry())))
    assert "(no decisions recorded)" in text


def test_emit_workflow_trace_counts():
    args = {
        "agent_trace": [{"iteration": 0, "tool": "x"}, {"iteration": 1, "tool": "y"}],
        "decisions": [{"choice": "a"}],
    }
    result = asyncio.run(emit_workflow_trace(args, AttachmentRegistry()))
    assert result["n_trace_entries"] == 2
    assert result["n_decisions"] == 1
    assert result["n_chars"] > 100


def test_emit_workflow_trace_NOT_in_llm_roster():
    """emit_workflow_trace is emitted SERVER-SIDE by the orchestrator now,
    not by the LLM — so it must NOT appear in the Main agent's tool list.
    Dispatcher access is preserved (covered by test_dispatch_chain_routes_emit_workflow_trace)."""
    roster = v2_tools.make_full_tool_roster()
    names = {t["name"] for t in roster}
    assert "emit_workflow_trace" not in names
    # The other emitters are still in the LLM roster.
    for emitter in ("emit_assembled_gb", "emit_parts_order", "emit_protocol"):
        assert emitter in names


def test_dispatch_chain_routes_emit_workflow_trace(tmp_path):
    result = asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_workflow_trace", {"assembly_method": "gibson"},
        AttachmentRegistry(), output_dir=str(tmp_path),
    ))
    assert result["ok"] is True
    assert (tmp_path / "workflow_trace.txt").exists()


def test_all_four_emitters_round_trip(tmp_path):
    """End-to-end: call all 4 emitters via the dispatch chain into one dir."""
    reg = AttachmentRegistry()
    aid = reg.register_product("test", "ACGT" * 200)

    # 1. assembled.gb
    asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_assembled_gb", {"attachment_id": aid}, reg, output_dir=str(tmp_path),
    ))
    # 2. parts_order.csv
    asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_parts_order",
        {"parts": [{"part_id": "p1", "name": "test", "length_bp": 800,
                    "origin": "user_supplied"}]},
        reg, output_dir=str(tmp_path),
    ))
    # 3. protocol.csv
    asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_protocol", {"assembly_method": "gibson"}, reg, output_dir=str(tmp_path),
    ))
    # 4. workflow_trace.txt
    asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_workflow_trace", {"assembly_method": "gibson"}, reg, output_dir=str(tmp_path),
    ))

    written = sorted(p.name for p in tmp_path.iterdir())
    assert written == ["assembled.gb", "parts_order.csv", "protocol.csv", "workflow_trace.txt"]
