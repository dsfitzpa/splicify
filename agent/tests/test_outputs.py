"""Tests for output emitters + the two-tier dispatch chain."""
import asyncio
import base64
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.outputs.assembled_gb import emit_assembled_gb
from agent_v2 import tools as v2_tools
from splicify_api.agent.agent_tools import AttachmentRegistry


def test_emit_assembled_gb_returns_file_envelope():
    reg = AttachmentRegistry()
    aid = reg.register_product("test_plasmid", "ACGT" * 200, circular=True)

    result = asyncio.run(emit_assembled_gb({"attachment_id": aid}, reg))
    assert result["ok"] is True
    assert result["file"]["fileName"] == "assembled.gb"
    decoded = base64.b64decode(result["file"]["dataBase64"]).decode("utf-8")
    assert decoded.startswith("LOCUS")
    assert "ORIGIN" in decoded
    assert result["length_bp"] == 800
    assert result["topology"] == "circular"
    assert result["written_path"] is None


def test_emit_assembled_gb_writes_to_disk(tmp_path):
    reg = AttachmentRegistry()
    aid = reg.register_product("test_plasmid", "ACGT" * 100)
    result = asyncio.run(emit_assembled_gb(
        {"attachment_id": aid}, reg, output_dir=str(tmp_path),
    ))
    assert result["written_path"] == str(tmp_path / "assembled.gb")
    body = (tmp_path / "assembled.gb").read_text()
    assert body.startswith("LOCUS")


def test_emit_assembled_gb_unknown_attachment_id():
    reg = AttachmentRegistry()
    result = asyncio.run(emit_assembled_gb({"attachment_id": "att_nope"}, reg))
    assert result["ok"] is False
    assert "unknown" in result["error"].lower()


def test_emit_assembled_gb_missing_arg():
    reg = AttachmentRegistry()
    result = asyncio.run(emit_assembled_gb({}, reg))
    assert result["ok"] is False


def test_emitter_tool_schema_well_formed():
    schema = v2_tools.EMIT_ASSEMBLED_GB_TOOL
    assert schema["name"] == "emit_assembled_gb"
    assert "description" in schema
    assert schema["input_schema"]["type"] == "object"
    assert "attachment_id" in schema["input_schema"]["properties"]
    assert schema["input_schema"]["required"] == ["attachment_id"]


def test_make_full_tool_roster_includes_v1_and_emitters():
    roster = v2_tools.make_full_tool_roster()
    names = {t["name"] for t in roster}
    assert "emit_assembled_gb" in names  # v2
    assert "annotate_attachment" in names  # v1
    assert "simulate_assembly" in names    # v1


def test_dispatch_with_emitters_routes_emitter_locally(tmp_path):
    reg = AttachmentRegistry()
    aid = reg.register_product("test", "ACGT" * 100)
    result = asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_assembled_gb", {"attachment_id": aid}, reg,
        output_dir=str(tmp_path),
    ))
    assert result["ok"] is True
    assert (tmp_path / "assembled.gb").exists()


def test_dispatch_with_emitters_falls_through_to_v1():
    """Names not in EMITTER_HANDLERS go to v1's dispatch_tool.

    v1's dispatch returns {"error": ...} for unknown tools (rather than raising),
    so we exercise that path with a clearly-bogus name.
    """
    reg = AttachmentRegistry()
    result = asyncio.run(v2_tools.dispatch_with_emitters(
        "definitely_not_a_real_tool_xyz", {}, reg,
    ))
    assert "error" in result
    assert "unknown tool" in result["error"].lower()
