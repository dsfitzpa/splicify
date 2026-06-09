"""Inventory fan-out + plasmid resolution tests."""
from __future__ import annotations

import agent_v2  # noqa: F401
from agent_v2.interpreter.plasmid_registry import PlasmidRegistry


def _env(name: str, **extra):
    base = {
        "sequence": "A" * 100,
        "annotations": [],
        "modules": [],
        "hierarchical_annotations": [],
        "interactions": [],
        "cloning_features": [],
    }
    base.update(extra)
    base["_name"] = name
    return base


def test_register_and_fanout_summary():
    reg = PlasmidRegistry()
    reg.register("p1", _env("alpha"), name="alpha.gb")
    reg.register("p2", _env("beta"), name="beta.gb")
    res = reg.fan_out("summary")
    assert res["n_plasmids_searched"] == 2
    assert len(res["results"]) == 2
    assert {r["plasmid_id"] for r in res["results"]} == {"p1", "p2"}


def test_resolve_plasmid_exact_name():
    reg = PlasmidRegistry()
    reg.register("p1", _env("alpha"), name="lentiCRISPR_v2")
    reg.register("p2", _env("beta"), name="pBSK")
    res = reg.resolve_plasmid("lentiCRISPR_v2")
    assert res["ok"] and res["method"] == "exact"
    assert res["matches"][0]["plasmid_id"] == "p1"


def test_resolve_plasmid_substring():
    reg = PlasmidRegistry()
    reg.register("p1", _env("alpha"), name="lentiCRISPR_v2_unannotated")
    res = reg.resolve_plasmid("lentiCRISPR")
    assert res["ok"] and res["method"] == "substring"


def test_resolve_plasmid_feature_overlap():
    """'the one with EGFP and PuroR' should match a plasmid whose
    annotations contain both EGFP and PuroR by distinguishing tokens."""
    reg = PlasmidRegistry()
    reg.register("p1", {
        "sequence": "A" * 100,
        "annotations": [
            {"name": "EGFP", "kb_data": {"gene_name": "EGFP"}, "start": 0, "end": 10},
            {"name": "PuroR", "kb_data": {"gene_name": "PuroR"}, "start": 20, "end": 30},
        ],
        "modules": [],
        "interactions": [],
        "hierarchical_annotations": [],
        "cloning_features": [],
    }, name="construct_a")
    reg.register("p2", {
        "sequence": "A" * 100,
        "annotations": [
            {"name": "Cas9", "kb_data": {"gene_name": "Cas9"}, "start": 0, "end": 10},
        ],
        "modules": [],
        "interactions": [],
        "hierarchical_annotations": [],
        "cloning_features": [],
    }, name="construct_b")
    res = reg.resolve_plasmid("the one with EGFP and PuroR")
    assert res["ok"]
    assert res["method"] == "feature_overlap"
    assert res["matches"][0]["plasmid_id"] == "p1"


def test_resolve_plasmid_no_match_returns_ok_false():
    reg = PlasmidRegistry()
    reg.register("p1", _env("alpha"), name="alpha.gb")
    res = reg.resolve_plasmid("zzzzz_doesnt_exist")
    assert not res["ok"]
    assert "registry_names" in res


def test_fan_out_records_empty_plasmids():
    reg = PlasmidRegistry()
    reg.register("p1", {
        "sequence": "A" * 100,
        "annotations": [{"name": "EGFP", "start": 0, "end": 10, "kb_data": {}}],
        "modules": [], "interactions": [], "hierarchical_annotations": [], "cloning_features": [],
    }, name="a")
    reg.register("p2", _env("b"), name="b")
    res = reg.fan_out("find_features", query="EGFP")
    assert res["ok"]
    assert len(res["results"]) == 1
    assert res["results"][0]["plasmid_id"] == "p1"
    assert "p2" in res["no_results_in"]


def test_fan_out_scoped_to_single_plasmid():
    reg = PlasmidRegistry()
    reg.register("p1", _env("a"), name="a")
    reg.register("p2", _env("b"), name="b")
    res = reg.fan_out("summary", plasmid_id="p1")
    assert res["scope"] == "single"
    assert len(res["results"]) == 1
    assert res["results"][0]["plasmid_id"] == "p1"


def test_fan_out_unknown_plasmid_id():
    reg = PlasmidRegistry()
    reg.register("p1", _env("a"))
    res = reg.fan_out("summary", plasmid_id="p404")
    assert not res["ok"]
    assert "p404" in res["reason"]
