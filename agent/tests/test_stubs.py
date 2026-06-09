"""Tests for the CRISPR stub + off-topic rejection."""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from agent_v2 import stub_crispr, rejection


def test_crispr_stub_envelope_shape():
    out = stub_crispr.respond("design a guide for EMX1")
    assert out["ok"] is True
    assert "coming soon" in out["reply"].lower()
    assert "sgrna golden gate" in out["reply"].lower()
    assert out["files"] is None
    assert out["viz"] is None
    assert out["n_tool_calls"] == 0
    assert out["error"] is None
    assert len(out["agent_trace"]) == 1
    assert out["agent_trace"][0]["tool"] == "stub_crispr"


def test_crispr_stub_default_arg():
    out = stub_crispr.respond()
    assert out["ok"] is True


def test_rejection_envelope_substitutes_shorthand():
    out = rejection.respond(shorthand="favorite color question")
    assert out["ok"] is True
    assert "favorite color question" in out["reply"]
    assert "claude.ai" in out["reply"].lower()
    assert "plasmid cloning" in out["reply"].lower()
    assert out["files"] is None
    assert out["viz"] is None
    assert out["agent_trace"][0]["tool"] == "rejection_template"


def test_rejection_default_no_shorthand():
    out = rejection.respond(shorthand="")
    assert out["ok"] is True
    assert "no summary available" in out["reply"].lower()


def test_rejection_records_reason_in_trace():
    out = rejection.respond(shorthand="hello", reason="off-topic-sports")
    trace = out["agent_trace"][0]
    assert "off-topic-sports" in trace["args_summary"]
