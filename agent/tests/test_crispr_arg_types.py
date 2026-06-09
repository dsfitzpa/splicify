"""Tests for _arg_to_int + every defensive int conversion in crispr_tools.

Sonnet 4.6 occasionally passes the entire response of a prior tool (typically
resolve_feature_position) as the argument to an int-typed input on a later
tool. The pegRNA + multi-target iter-41 live smoke hit this on the
design_pegrnas_tool edit_start argument. _arg_to_int handles the common
bad shapes (dict / list / tuple / string) so the pipeline keeps working
even when Claude mis-formats its tool args.
"""
import asyncio
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest

import agent_v2  # noqa: F401
from agent_v2.crispr_tools import (
    _arg_to_int,
    design_guides_tool,
    design_pegrnas_tool,
    design_primers_tool,
)
from splicify_api.agent.agent_tools import AttachmentRegistry


# ---------------------------------------------------------------------------
# _arg_to_int unit tests
# ---------------------------------------------------------------------------
def test_arg_to_int_passes_through_int():
    assert _arg_to_int(42) == 42
    assert _arg_to_int(0) == 0
    assert _arg_to_int(-3) == -3


def test_arg_to_int_truncates_float():
    assert _arg_to_int(3.7) == 3
    assert _arg_to_int(-1.9) == -1


def test_arg_to_int_handles_bool_as_int():
    assert _arg_to_int(True) == 1
    assert _arg_to_int(False) == 0


def test_arg_to_int_parses_numeric_string():
    assert _arg_to_int("42") == 42
    assert _arg_to_int(" 42 ") == 42
    assert _arg_to_int("3.7") == 3   # float-string path


def test_arg_to_int_returns_default_on_empty_or_garbage_string():
    assert _arg_to_int("", default=-1) == -1
    assert _arg_to_int("xyz", default=99) == 99


def test_arg_to_int_unwraps_tuple_and_list():
    """primer3 returns left_pos as Tuple[int, int] = (start, length)."""
    assert _arg_to_int([10, 20]) == 10
    assert _arg_to_int((10, 20)) == 10
    assert _arg_to_int([100]) == 100


def test_arg_to_int_unwraps_resolve_feature_position_dict():
    """Most common Sonnet 4.6 mistake: pass the whole resolve_* response."""
    fake = {"ok": True, "plasmid_position": 14694, "codon": "CGC",
            "amino_acid": "R", "feature_strand": "-",
            "feature_start": 1724, "feature_end": 18955}
    assert _arg_to_int(fake) == 14694


def test_arg_to_int_falls_back_through_known_keys():
    assert _arg_to_int({"position": 100}) == 100
    assert _arg_to_int({"start": 200}) == 200
    assert _arg_to_int({"value": 300}) == 300


def test_arg_to_int_picks_any_int_value_when_no_known_key():
    """Dict with no known key but a single int-typed value works as last resort."""
    assert _arg_to_int({"obscure_field": 42}) == 42


def test_arg_to_int_returns_default_on_none_or_empty():
    assert _arg_to_int(None, default=5) == 5
    assert _arg_to_int(None) is None
    assert _arg_to_int([], default=5) == 5
    assert _arg_to_int({}, default=5) == 5


def test_arg_to_int_nested_dict_with_position_key():
    """Position deep in a more-complex dict still resolves."""
    payload = {"ok": True, "feature_name": "KEAP1",
                "plasmid_position": 14694, "spans_intron": False}
    assert _arg_to_int(payload) == 14694


# ---------------------------------------------------------------------------
# design_pegrnas_tool — the exact iter 41 live-smoke failure
# ---------------------------------------------------------------------------
def _registry():
    reg = AttachmentRegistry()
    reg.register_product("test", "ACGT" * 5000, circular=False)
    return reg


def test_pegrnas_tool_recovers_from_dict_edit_start(monkeypatch):
    """Reproduces the iter 41 live-smoke failure mode: Claude passed the
    full resolve_feature_position result dict instead of plasmid_position."""
    captured = {}

    def fake_v1(**kw):
        captured.update(kw)
        return {"ok": True, "pegrnas": [], "summary": {"n_returned": 0}}

    import splicify_api.pegrna_designer as pd_mod
    monkeypatch.setattr(pd_mod, "design_pegrnas", fake_v1)

    rfp_response = {"ok": True, "plasmid_position": 1234, "codon": "CGC",
                     "amino_acid": "R", "feature_strand": "-",
                     "feature_start": 1000, "feature_end": 2000}
    out = asyncio.run(design_pegrnas_tool(
        {"attachment_id": "att_product_1",
         "edit_start": rfp_response, "edit_end": rfp_response,
         "alt": "TGC", "edit_type": "substitution"},
        _registry(),
    ))
    assert out["ok"] is True
    # _arg_to_int unwrapped both to the plasmid_position int.
    assert captured["edit_start_1based"] == 1234
    assert captured["edit_end_1based"] == 1234


def test_pegrnas_tool_recovers_from_list_edit_coords(monkeypatch):
    """Tuple/list shape (e.g. [start, length]) also unwraps."""
    captured = {}
    import splicify_api.pegrna_designer as pd_mod
    monkeypatch.setattr(pd_mod, "design_pegrnas",
                          lambda **kw: (captured.update(kw) or
                                        {"ok": True, "pegrnas": [], "summary": {}}))
    out = asyncio.run(design_pegrnas_tool(
        {"attachment_id": "att_product_1",
         "edit_start": [1234, 3], "edit_end": [1236, 3],
         "alt": "TGC"}, _registry(),
    ))
    assert out["ok"] is True
    assert captured["edit_start_1based"] == 1234
    assert captured["edit_end_1based"] == 1236


def test_pegrnas_tool_clean_error_on_missing_edit_coords():
    """No edit_start at all -> the new error message tells the LLM exactly
    what to pass next time."""
    out = asyncio.run(design_pegrnas_tool(
        {"attachment_id": "att_product_1", "alt": "TGC"}, _registry(),
    ))
    assert out["ok"] is False
    assert "plasmid_position" in out["error"]


# ---------------------------------------------------------------------------
# design_guides_tool — dict-region defence
# ---------------------------------------------------------------------------
def test_guides_tool_recovers_from_dict_region(monkeypatch):
    captured = {}

    def fake_v1(**kw):
        captured.update(kw)
        return {"ok": True, "guides": [], "summary": {}}

    import splicify_api.guide_designer as gd
    monkeypatch.setattr(gd, "design_guides", fake_v1)

    fake = {"plasmid_position": 100, "codon": "ATG"}
    out = asyncio.run(design_guides_tool(
        {"attachment_id": "att_product_1",
         "region_start": fake, "region_end": 200},
        _registry(),
    ))
    assert out["ok"] is True
    assert captured["region_start"] == 100
    assert captured["region_end"] == 200


# ---------------------------------------------------------------------------
# design_primers_tool — primer3 (start, length) tuple AND LLM dict args
# ---------------------------------------------------------------------------
def test_primers_tool_unwraps_primer3_tuple_left_pos(monkeypatch):
    """primer3 left_pos = (start_in_template, length)."""
    def fake_v1(req):
        return {"ok": True, "left_pos": [10, 20], "right_pos": [290, 20],
                 "application": req.application}
    import splicify_api.pcr as pcr_mod
    monkeypatch.setattr(pcr_mod, "design_primers", fake_v1)

    out = asyncio.run(design_primers_tool(
        {"attachment_id": "att_product_1",
         "region_start": 1000, "region_end": 1500,
         "application": "sanger"}, _registry(),
    ))
    assert out["ok"] is True
    # left_pos 10 in template -> plasmid 1000 + 10 - 1 = 1009.
    assert out["left_pos_plasmid"] == 1009
    assert out["right_pos_plasmid"] == 1289


def test_primers_tool_recovers_from_dict_region(monkeypatch):
    def fake_v1(req):
        return {"ok": True, "left_pos": [0, 20], "right_pos": [10, 20]}
    import splicify_api.pcr as pcr_mod
    monkeypatch.setattr(pcr_mod, "design_primers", fake_v1)

    fake_anchor = {"plasmid_position": 1000}
    out = asyncio.run(design_primers_tool(
        {"attachment_id": "att_product_1",
         "region_start": fake_anchor, "region_end": 1500},
        _registry(),
    ))
    assert out["ok"] is True
    assert out["region_start_plasmid"] == 1000


def test_primers_tool_recovers_from_dict_excluded(monkeypatch):
    captured = {}
    def fake_v1(req):
        captured["excluded_start"] = req.excluded_start
        captured["excluded_length"] = req.excluded_length
        return {"ok": True, "left_pos": [0, 20], "right_pos": [10, 20]}
    import splicify_api.pcr as pcr_mod
    monkeypatch.setattr(pcr_mod, "design_primers", fake_v1)

    anchor_start = {"plasmid_position": 1100}
    anchor_end = {"plasmid_position": 1140}
    out = asyncio.run(design_primers_tool(
        {"attachment_id": "att_product_1",
         "region_start": 1000, "region_end": 1500,
         "excluded_start": anchor_start, "excluded_end": anchor_end},
        _registry(),
    ))
    assert out["ok"] is True
    # 1100 - 1000 = 100 (template-relative start); length 1140-1100+1 = 41.
    assert captured["excluded_start"] == 100
    assert captured["excluded_length"] == 41
