"""Tests for emit_parts_order — source classification, CSV format, on-disk write."""
import asyncio
import base64
import csv
import io
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from agent_v2.outputs.parts_order_csv import emit_parts_order, _classify_source
from agent_v2 import tools as v2_tools
from splicify_api.agent.agent_tools import AttachmentRegistry


def test_classify_source_addgene_wins_over_length():
    # Even a 5-bp Addgene record routes to order_addgene
    assert _classify_source({"addgene_id": "12345", "length_bp": 5}) == "order_addgene"


def test_classify_source_inventory():
    assert _classify_source({"origin": "inventory", "length_bp": 2686}) == "inventory"
    assert _classify_source({"origin": "user_supplied", "length_bp": 2686}) == "inventory"


def test_classify_source_oligo_under_40_bp():
    assert _classify_source({"length_bp": 22}) == "order_oligo"
    assert _classify_source({"origin": "designed_oligo", "length_bp": 60}) == "order_oligo"


def test_classify_source_synthesis_default():
    assert _classify_source({"length_bp": 720}) == "order_synthesis"
    assert _classify_source({"origin": "knowledge_base", "length_bp": 720}) == "order_synthesis"


def test_emit_parts_order_csv_header_and_rows():
    parts = [
        {"part_id": "p1", "name": "EGFP", "length_bp": 720, "origin": "knowledge_base"},
        {"part_id": "p2", "name": "T7-fwd", "length_bp": 22, "origin": "designed_oligo"},
        {"part_id": "p3", "name": "pUC19", "length_bp": 2686, "origin": "inventory"},
        {"part_id": "p4", "name": "px330", "length_bp": 8500, "addgene_id": "42230"},
    ]
    result = asyncio.run(emit_parts_order({"parts": parts}, AttachmentRegistry()))
    assert result["ok"] is True
    assert result["n_parts"] == 4
    assert result["source_counts"] == {
        "order_synthesis": 1, "order_oligo": 1, "inventory": 1, "order_addgene": 1,
    }

    csv_text = base64.b64decode(result["file"]["dataBase64"]).decode()
    rows = list(csv.reader(io.StringIO(csv_text)))
    assert rows[0] == ["part_id", "name", "source", "length_bp", "vendor",
                       "est_cost_usd", "est_lead_time_days", "notes"]
    assert len(rows) == 5  # header + 4 parts
    # spot-check classifications + default vendors
    by_id = {r[0]: r for r in rows[1:]}
    assert by_id["p1"][2] == "order_synthesis"
    assert by_id["p1"][4] == "Twist gBlock"
    assert by_id["p2"][2] == "order_oligo"
    assert by_id["p2"][4] == "IDT oligo"
    assert by_id["p3"][2] == "inventory"
    assert by_id["p4"][2] == "order_addgene"
    assert "42230" in by_id["p4"][4]


def test_emit_parts_order_writes_to_disk(tmp_path):
    parts = [{"part_id": "p1", "name": "EGFP", "length_bp": 720}]
    result = asyncio.run(emit_parts_order(
        {"parts": parts}, AttachmentRegistry(), output_dir=str(tmp_path),
    ))
    assert result["written_path"] == str(tmp_path / "parts_order.csv")
    body = (tmp_path / "parts_order.csv").read_text()
    assert body.startswith("part_id,name,source")


def test_emit_parts_order_empty_parts_returns_header_only():
    result = asyncio.run(emit_parts_order({"parts": []}, AttachmentRegistry()))
    assert result["ok"] is True
    csv_text = base64.b64decode(result["file"]["dataBase64"]).decode()
    rows = list(csv.reader(io.StringIO(csv_text)))
    assert len(rows) == 1
    assert result["source_counts"] == {}


def test_emit_parts_order_overrides_vendor_cost_lead():
    parts = [{
        "part_id": "p", "name": "custom", "length_bp": 500,
        "vendor_hint": "GenScript", "cost_usd": 250.0, "lead_time_days": 21,
        "notes": "custom protein, host-optimized",
    }]
    result = asyncio.run(emit_parts_order({"parts": parts}, AttachmentRegistry()))
    csv_text = base64.b64decode(result["file"]["dataBase64"]).decode()
    rows = list(csv.reader(io.StringIO(csv_text)))
    row = rows[1]
    assert row[4] == "GenScript"
    assert row[5] == "250.00"
    assert row[6] == "21"
    assert "host-optimized" in row[7]


def test_emit_parts_order_quotes_commas_in_names():
    parts = [{"part_id": "p", "name": "Bsmb-I, version 2", "length_bp": 500}]
    result = asyncio.run(emit_parts_order({"parts": parts}, AttachmentRegistry()))
    csv_text = base64.b64decode(result["file"]["dataBase64"]).decode()
    # csv module handles quoting; round-trip via csv.reader to verify
    rows = list(csv.reader(io.StringIO(csv_text)))
    assert rows[1][1] == "Bsmb-I, version 2"


def test_emit_parts_order_bad_input_returns_error():
    result = asyncio.run(emit_parts_order({"parts": "not a list"}, AttachmentRegistry()))
    assert result["ok"] is False


def test_emit_parts_order_tool_schema_in_full_roster():
    roster = v2_tools.make_full_tool_roster()
    names = {t["name"] for t in roster}
    assert "emit_parts_order" in names
    assert "emit_assembled_gb" in names


def test_dispatch_chain_routes_emit_parts_order(tmp_path):
    parts = [{"part_id": "p1", "name": "x", "length_bp": 100}]
    result = asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_parts_order", {"parts": parts}, AttachmentRegistry(),
        output_dir=str(tmp_path),
    ))
    assert result["ok"] is True
    assert (tmp_path / "parts_order.csv").exists()
