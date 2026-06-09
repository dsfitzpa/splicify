"""Tests for the three CRISPR + primer design tool wrappers in agent_v2.crispr_tools.

v1 designers are monkeypatched so the suite stays fast (~no primer3, no pickled
PE3 XGBoost load) and tests are deterministic.
"""
import asyncio
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest

import agent_v2  # noqa: F401 - triggers path shim
from agent_v2.crispr_tools import (
    design_guides_tool,
    design_pegrnas_tool,
    design_primers_tool,
)
from agent_v2 import tools as v2_tools
from splicify_api.agent.agent_tools import AttachmentRegistry


def _registry_with_seq(seq: str = "ACGT" * 1000):
    reg = AttachmentRegistry()
    reg.register_product("test", seq, circular=True)
    return reg


# ---------------------------------------------------------------------------
# design_guides
# ---------------------------------------------------------------------------
def test_design_guides_unknown_attachment():
    reg = AttachmentRegistry()
    out = asyncio.run(design_guides_tool(
        {"attachment_id": "att_nope", "region_start": 1, "region_end": 100}, reg,
    ))
    assert out["ok"] is False
    assert "unknown" in out["error"].lower()


def test_design_guides_invalid_region():
    reg = _registry_with_seq("ACGT" * 100)
    out = asyncio.run(design_guides_tool(
        {"attachment_id": "att_product_1", "region_start": 100, "region_end": 50}, reg,
    ))
    assert out["ok"] is False
    assert "invalid region" in out["error"].lower()


def test_design_guides_happy_path(monkeypatch):
    reg = _registry_with_seq("ACGT" * 250)  # 1000 bp

    captured = {}

    def fake_v1(*, sequence, region_start, region_end, **kw):
        captured["sequence_len"] = len(sequence)
        captured["region_start"] = region_start
        captured["region_end"] = region_end
        captured.update(kw)
        return {
            "ok": True,
            "guides": [
                {
                    "name": "guide_001_TGTTCA",
                    "spacer": "TGTTCAACGAAATGAAAATC",
                    "pam": "TGG",
                    "start": 57, "end": 77, "direction": 1,
                    "score": 62.0,
                    "score_components": {"intercept": 0.5},
                    "score_method": "doench2014",
                    "context_30mer": "GCG...",
                    "gc_fraction": 0.45,
                    "max_homopolymer": 3,
                    "n_offtargets": 0,
                },
            ],
            "summary": {"n_candidates": 1, "n_returned": 1},
        }

    import splicify_api.guide_designer as gd
    monkeypatch.setattr(gd, "design_guides", fake_v1)

    out = asyncio.run(design_guides_tool(
        {
            "attachment_id": "att_product_1",
            "region_start": 35, "region_end": 101,
            "pam": "NGG", "guide_length": 20, "max_guides": 5,
        }, reg,
    ))
    assert out["ok"] is True
    assert out["attachment_id"] == "att_product_1"
    assert len(out["guides"]) == 1
    g = out["guides"][0]
    # Digest preserves the useful fields.
    assert g["spacer"] == "TGTTCAACGAAATGAAAATC"
    assert g["pam"] == "TGG"
    assert g["score"] == 62.0
    assert g["score_method"] == "doench2014"
    # Heavy fields stripped.
    assert "context_30mer" not in g
    assert "score_components" not in g
    # v1 was called with the right args.
    assert captured["region_start"] == 35
    assert captured["region_end"] == 101
    assert captured["pam"] == "NGG"
    assert captured["sequence_len"] == 1000


# ---------------------------------------------------------------------------
# design_pegrnas
# ---------------------------------------------------------------------------
def test_design_pegrnas_unknown_attachment():
    reg = AttachmentRegistry()
    out = asyncio.run(design_pegrnas_tool(
        {"attachment_id": "att_nope", "edit_start": 100, "edit_end": 102, "alt": "TAA"}, reg,
    ))
    assert out["ok"] is False
    assert "unknown" in out["error"].lower()


def test_design_pegrnas_missing_edit_coords():
    reg = _registry_with_seq()
    out = asyncio.run(design_pegrnas_tool(
        {"attachment_id": "att_product_1"}, reg,
    ))
    assert out["ok"] is False
    assert "edit_start" in out["error"]


def test_design_pegrnas_happy_path(monkeypatch):
    reg = _registry_with_seq()
    captured = {}

    def fake_v1(*, sequence, edit_start_1based, edit_end_1based, alt, edit_type, n_results, use_pe3, params=None):
        captured.update({
            "edit_start_1based": edit_start_1based,
            "edit_end_1based": edit_end_1based,
            "alt": alt, "edit_type": edit_type,
            "n_results": n_results, "use_pe3": use_pe3,
            "sequence_len": len(sequence),
        })
        return {
            "ok": True,
            "pegrnas": [
                {
                    "rank": 1,
                    "name": "pegRNA_1_TGTTCA",
                    "predicted_efficiency": 0.412,
                    "spacer": "TGTTCAACGAAATGAAAATC", "pam": "TGG",
                    "spacer_start": 1207, "spacer_end": 1227, "direction": 1,
                    "cas9_score": 62.0,
                    "rtt": "GCGATTAATCG", "rtt_length": 11, "rtt_gc": 0.45,
                    "pbs": "TTCAACGAAATGAAA", "pbs_length": 13, "pbs_gc": 0.40,
                    "scaffold": "GTTTTAGAGCT...",
                    "full_pegrna": "ACGT" * 30, "full_pegrna_length": 120,
                    "is_dpam": False, "is_pe3b": True,
                    "ngrna": {"spacer": "AAAGCT", "start": 1310, "end": 1330, "strand": "-", "cas9_score": 54.3},
                    "edit_type": "substitution",
                    "edit_ref": "AAA", "edit_alt": "TAA",
                    "edit_start_1based": 1234, "edit_end_1based": 1236,
                    "score_components": {"cas9_score": 62.0},
                },
            ],
            "summary": {"n_returned": 1, "use_pe3": True},
        }

    import splicify_api.pegrna_designer as pd_mod
    monkeypatch.setattr(pd_mod, "design_pegrnas", fake_v1)

    out = asyncio.run(design_pegrnas_tool(
        {
            "attachment_id": "att_product_1",
            "edit_start": 1234, "edit_end": 1236,
            "alt": "TAA", "edit_type": "substitution",
            "n_results": 3, "use_pe3": True,
        }, reg,
    ))
    assert out["ok"] is True
    assert len(out["pegrnas"]) == 1
    p = out["pegrnas"][0]
    assert p["spacer"] == "TGTTCAACGAAATGAAAATC"
    assert p["full_pegrna_length"] == 120
    assert p["is_pe3b"] is True
    assert p["ngrna"]["spacer"] == "AAAGCT"
    # Heavy fields stripped from digest.
    assert "score_components" not in p
    assert captured["edit_start_1based"] == 1234
    assert captured["alt"] == "TAA"
    assert captured["use_pe3"] is True


# ---------------------------------------------------------------------------
# design_primers
# ---------------------------------------------------------------------------
def test_design_primers_unknown_attachment():
    reg = AttachmentRegistry()
    out = asyncio.run(design_primers_tool(
        {"attachment_id": "att_nope", "region_start": 1, "region_end": 500}, reg,
    ))
    assert out["ok"] is False
    assert "unknown" in out["error"].lower()


def test_design_primers_invalid_region():
    reg = _registry_with_seq("ACGT" * 100)  # 400 bp
    out = asyncio.run(design_primers_tool(
        {"attachment_id": "att_product_1", "region_start": 500, "region_end": 600}, reg,
    ))
    assert out["ok"] is False
    assert "invalid region" in out["error"].lower()


def test_design_primers_excluded_region_outside():
    reg = _registry_with_seq("ACGT" * 1000)  # 4000 bp
    out = asyncio.run(design_primers_tool(
        {
            "attachment_id": "att_product_1",
            "region_start": 100, "region_end": 600,
            "excluded_start": 800, "excluded_end": 850,
        }, reg,
    ))
    assert out["ok"] is False
    assert "outside" in out["error"].lower()


def test_design_primers_happy_path_coord_translation(monkeypatch):
    reg = _registry_with_seq("ACGT" * 1000)
    captured = {}

    def fake_v1(req):
        captured["fragments_in_len"] = len(req.fragments_in)
        captured["excluded_start"] = req.excluded_start
        captured["excluded_length"] = req.excluded_length
        captured["application"] = req.application
        captured["product_size_min"] = req.product_size_min
        captured["product_size_max"] = req.product_size_max
        return {
            "ok": True,
            "pair_index": 0, "pair_penalty": 0.3,
            "product_size": 280,
            "left_pos": [10, 20], "right_pos": [290, 20],
            "left_primer": "ACGTACGTACGTACGTACGT",
            "right_primer": "TGCATGCATGCATGCATGCA",
            "left_annealing": "ACGTACGTACGTACGTACGT",
            "right_annealing": "TGCATGCATGCATGCATGCA",
            "left_adapter": "", "right_adapter": "",
            "sanger_scores": [{"score": 72.0, "rating": "good"}, {"score": 68.0, "rating": "good"}],
            "selection_method": "sanger_aware",
            "selection_rationale": "Picked candidate #0 ...",
            "num_candidates_considered": 9,
            "application": "sanger",
        }

    import splicify_api.pcr as pcr_mod
    monkeypatch.setattr(pcr_mod, "design_primers", fake_v1)

    out = asyncio.run(design_primers_tool(
        {
            "attachment_id": "att_product_1",
            "region_start": 1000, "region_end": 1500,
            "excluded_start": 1230, "excluded_end": 1270,
            "application": "sanger",
        }, reg,
    ))
    assert out["ok"] is True
    # Region slicing: 1500 - 1000 + 1 = 501 bp template.
    assert captured["fragments_in_len"] == 501
    # Excluded coords translated to template-relative (1230 - 1000 = 230).
    assert captured["excluded_start"] == 230
    # Length: 1270 - 1230 + 1 = 41.
    assert captured["excluded_length"] == 41
    # Sanger default product size band.
    assert captured["product_size_min"] == 250
    assert captured["product_size_max"] == 500
    # Returned primer positions translated back to plasmid coords:
    # left_pos=10 in template -> 1000 + 10 - 1 = 1009 on plasmid.
    assert out["left_pos_plasmid"] == 1009
    assert out["right_pos_plasmid"] == 1289
    assert out["region_start_plasmid"] == 1000
    assert out["region_end_plasmid"] == 1500
    assert out["attachment_id"] == "att_product_1"


def test_design_primers_application_defaults(monkeypatch):
    """illumina + fragment hit their own default product-size bands."""
    reg = _registry_with_seq("ACGT" * 1000)
    seen = []

    def fake_v1(req):
        seen.append({
            "application": req.application,
            "min": req.product_size_min, "max": req.product_size_max,
        })
        return {"ok": True, "left_pos": [0, 20], "right_pos": [10, 20], "application": req.application}

    import splicify_api.pcr as pcr_mod
    monkeypatch.setattr(pcr_mod, "design_primers", fake_v1)

    asyncio.run(design_primers_tool(
        {"attachment_id": "att_product_1", "region_start": 1, "region_end": 600, "application": "illumina"}, reg,
    ))
    asyncio.run(design_primers_tool(
        {"attachment_id": "att_product_1", "region_start": 1, "region_end": 600, "application": "fragment"}, reg,
    ))
    assert seen[0] == {"application": "illumina", "min": 150, "max": 290}
    assert seen[1] == {"application": "fragment", "min": 100, "max": 300}


# ---------------------------------------------------------------------------
# Dispatch chain + roster
# ---------------------------------------------------------------------------
def test_dispatch_routes_design_guides_locally(monkeypatch):
    """dispatch_with_emitters routes design_guides to our local handler."""
    reg = _registry_with_seq("ACGT" * 250)

    def fake_v1(**kw):
        return {"ok": True, "guides": [], "summary": {"n_returned": 0}}

    import splicify_api.guide_designer as gd
    monkeypatch.setattr(gd, "design_guides", fake_v1)

    out = asyncio.run(v2_tools.dispatch_with_emitters(
        "design_guides",
        {"attachment_id": "att_product_1", "region_start": 10, "region_end": 200},
        reg,
    ))
    assert out["ok"] is True
    assert out["attachment_id"] == "att_product_1"


def test_full_tool_roster_includes_crispr_tools():
    roster = v2_tools.make_full_tool_roster()
    names = {t["name"] for t in roster}
    assert "design_guides" in names
    assert "design_pegrnas" in names
    assert "design_primers" in names
    # Existing entries still present.
    assert "emit_assembled_gb" in names
    assert "resolve_feature_position" in names
