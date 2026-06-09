"""Tests for emit_protocol — per-method templates + CSV format + overrides."""
import asyncio
import base64
import csv
import io
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from agent_v2.outputs.protocol_csv import emit_protocol, PROTOCOLS
from agent_v2 import tools as v2_tools
from splicify_api.agent.agent_tools import AttachmentRegistry


def _csv_rows(result):
    text = base64.b64decode(result["file"]["dataBase64"]).decode()
    return list(csv.reader(io.StringIO(text)))


def test_emit_protocol_gibson_template_categories():
    result = asyncio.run(emit_protocol({"assembly_method": "gibson"}, AttachmentRegistry()))
    assert result["ok"] is True
    assert result["method"] == "gibson"
    cats = result["categories"]
    for c in ("PCR", "ligation", "transform", "verify", "sanger"):
        assert c in cats


def test_emit_protocol_gateway_has_clonase_step():
    result = asyncio.run(emit_protocol({"assembly_method": "gateway"}, AttachmentRegistry()))
    rows = _csv_rows(result)
    body = "\n".join(",".join(r) for r in rows)
    assert "Clonase" in body or "clonase" in body


def test_emit_protocol_restriction_template():
    result = asyncio.run(emit_protocol({"assembly_method": "restriction"}, AttachmentRegistry()))
    cats = result["categories"]
    assert "digest" in cats
    assert "ligation" in cats


def test_emit_protocol_sdm_has_dpni_step():
    result = asyncio.run(emit_protocol({"assembly_method": "sdm"}, AttachmentRegistry()))
    rows = _csv_rows(result)
    body = "\n".join(",".join(r) for r in rows)
    assert "DpnI" in body


def test_emit_protocol_sgrna_gg_has_anneal_and_pnk():
    result = asyncio.run(emit_protocol({"assembly_method": "sgrna_gg"}, AttachmentRegistry()))
    rows = _csv_rows(result)
    body = "\n".join(",".join(r) for r in rows)
    assert "anneal" in body.lower() or "Annealed" in body
    assert "PNK" in body


def test_emit_protocol_golden_gate_template():
    result = asyncio.run(emit_protocol({"assembly_method": "golden_gate"}, AttachmentRegistry()))
    rows = _csv_rows(result)
    body = "\n".join(",".join(r) for r in rows)
    assert "Type-IIs" in body or "BsaI" in body


def test_emit_protocol_unknown_method_falls_back_to_generic():
    result = asyncio.run(emit_protocol({"assembly_method": "made_up_method"}, AttachmentRegistry()))
    assert result["ok"] is True
    # generic + transform tail = at least the 5 transform-block steps + 1 generic
    assert result["n_steps"] >= 6


def test_emit_protocol_csv_header_and_step_numbering():
    result = asyncio.run(emit_protocol({"assembly_method": "gibson"}, AttachmentRegistry()))
    rows = _csv_rows(result)
    assert rows[0] == ["step_num", "category", "inputs", "output", "instrument",
                       "time_min", "temp_C", "reagents", "notes"]
    # step numbers are 1..N consecutive
    nums = [int(r[0]) for r in rows[1:]]
    assert nums == list(range(1, len(rows)))


def test_emit_protocol_total_time_is_summed():
    result = asyncio.run(emit_protocol({"assembly_method": "gibson"}, AttachmentRegistry()))
    expected = sum(int(s.get("time_min") or 0) for s in PROTOCOLS["gibson"])
    expected += sum(int(s.get("time_min") or 0) for s in PROTOCOLS["gibson"][3:5])  # tail in main template
    # Actually the result is the sum across all generated steps. Recompute from rows:
    rows = _csv_rows(result)
    actual_total = sum(int(r[5] or 0) for r in rows[1:])
    assert result["total_time_min"] == actual_total


def test_emit_protocol_writes_to_disk(tmp_path):
    result = asyncio.run(emit_protocol(
        {"assembly_method": "gibson"}, AttachmentRegistry(),
        output_dir=str(tmp_path),
    ))
    body = (tmp_path / "protocol.csv").read_text()
    assert body.startswith("step_num,category")


def test_emit_protocol_custom_step_overrides():
    custom = [{"step_num": 1, "time_min": 999, "notes": "extended cycling"}]
    result = asyncio.run(emit_protocol(
        {"assembly_method": "gibson", "custom_steps": custom},
        AttachmentRegistry(),
    ))
    rows = _csv_rows(result)
    # step 1 (PCR) should now show 999 min and the custom note
    step1 = rows[1]
    assert step1[0] == "1"
    assert step1[5] == "999"
    assert "extended cycling" in step1[8]


def test_emit_protocol_custom_step_appended_when_no_step_num():
    custom = [{"category": "verify",
               "inputs": "Final plasmid + 100 ng input",
               "output": "Whole-plasmid sequencing read",
               "instrument": "(external)",
               "time_min": 1440, "temp_C": "n/a",
               "reagents": "Plasmidsaurus tube",
               "notes": "Confirm whole-plasmid before starting wet lab"}]
    result = asyncio.run(emit_protocol(
        {"assembly_method": "gibson", "custom_steps": custom},
        AttachmentRegistry(),
    ))
    rows = _csv_rows(result)
    # last row contains the appended step
    last = rows[-1]
    assert "Plasmidsaurus" in last[7]


def test_emit_protocol_tool_schema_in_full_roster():
    roster = v2_tools.make_full_tool_roster()
    names = {t["name"] for t in roster}
    assert "emit_protocol" in names


def test_dispatch_chain_routes_emit_protocol(tmp_path):
    result = asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_protocol", {"assembly_method": "gibson"}, AttachmentRegistry(),
        output_dir=str(tmp_path),
    ))
    assert result["ok"] is True
    assert (tmp_path / "protocol.csv").exists()
