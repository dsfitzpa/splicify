"""Tests for resolve_feature_position — annotate_cached mocked, coordinate math verified."""
import asyncio
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest

import agent_v2  # noqa: F401
from agent_v2.feature_resolver import resolve_feature_position, _rc, _slice_circular
from splicify_api.agent.agent_tools import AttachmentRegistry


@pytest.fixture
def patched_annotation(monkeypatch):
    """Mock annotate_cached to return a small canned annotation."""
    async def fake(seq, *, circular=True, depth="full"):
        return {
            "annotations": [
                {"name": "Cas9", "start": 5111, "end": 9215, "strand": 1},
                {"name": "AmpR", "start": 13882, "end": 14741, "strand": -1},
                {"name": "U6 promoter", "start": 2607, "end": 2849, "strand": 1},
            ],
        }
    import splicify_api.annotation_cache as ac
    monkeypatch.setattr(ac, "annotate_cached", fake)


def _make_registry_with_cas9_at_5111():
    """Build a synthetic plasmid where Cas9 starts at position 5111 with codon GAT (Asp)."""
    seq = "N" * 5111 + "GAT" + "N" * (14873 - 5114)  # 14873 bp circular
    reg = AttachmentRegistry()
    reg.register_product("lenti", seq, circular=True)
    return reg


def test_resolve_aa_residue_plus_strand(patched_annotation):
    reg = _make_registry_with_cas9_at_5111()
    result = asyncio.run(resolve_feature_position(
        {"attachment_id": "att_product_1", "feature_name": "Cas9",
         "kind": "aa_residue", "offset": 1},
        reg,
    ))
    assert result["ok"] is True
    assert result["feature_name"] == "Cas9"
    assert result["feature_strand"] == "+"
    assert result["plasmid_position"] == 5111
    assert result["codon"] == "GAT"
    assert result["amino_acid"] == "D"


def test_resolve_aa_residue_d10a_codon(patched_annotation):
    """Residue 10 of Cas9 on + strand: position 5111 + (10-1)*3 = 5138."""
    seq = bytearray(b"N" * 14873)
    seq[5138:5141] = b"GAT"  # D10 codon
    reg = AttachmentRegistry()
    reg.register_product("lenti", seq.decode(), circular=True)
    result = asyncio.run(resolve_feature_position(
        {"attachment_id": "att_product_1", "feature_name": "Cas9",
         "kind": "aa_residue", "offset": 10},
        reg,
    ))
    assert result["plasmid_position"] == 5138
    assert result["codon"] == "GAT"
    assert result["amino_acid"] == "D"


def test_resolve_aa_residue_minus_strand_reverse_complemented(patched_annotation):
    """AmpR is on - strand, end=14741. Residue 1 codon should be at end-3=14738 raw,
    then reverse-complemented to give the sense codon."""
    seq = bytearray(b"N" * 14873)
    # Sense codon at residue 1 of AmpR should be e.g. ATG (Met start codon).
    # On the genome (+ strand), that means the bases at 14738-14740 are the reverse-complement: CAT.
    seq[14738:14741] = b"CAT"  # rev-comp of ATG
    reg = AttachmentRegistry()
    reg.register_product("lenti", seq.decode(), circular=True)
    result = asyncio.run(resolve_feature_position(
        {"attachment_id": "att_product_1", "feature_name": "AmpR",
         "kind": "aa_residue", "offset": 1},
        reg,
    ))
    assert result["ok"] is True
    assert result["feature_strand"] == "-"
    assert result["codon"] == "ATG"
    assert result["amino_acid"] == "M"


def test_resolve_feature_start_and_end(patched_annotation):
    reg = _make_registry_with_cas9_at_5111()
    r1 = asyncio.run(resolve_feature_position(
        {"attachment_id": "att_product_1", "feature_name": "Cas9",
         "kind": "feature_start", "offset": 0}, reg,
    ))
    assert r1["plasmid_position"] == 5111
    r2 = asyncio.run(resolve_feature_position(
        {"attachment_id": "att_product_1", "feature_name": "Cas9",
         "kind": "feature_end", "offset": 0}, reg,
    ))
    assert r2["plasmid_position"] == 9215


def test_resolve_bp_offset_plus_strand(patched_annotation):
    reg = _make_registry_with_cas9_at_5111()
    result = asyncio.run(resolve_feature_position(
        {"attachment_id": "att_product_1", "feature_name": "U6 promoter",
         "kind": "bp_offset", "offset": 100},
        reg,
    ))
    assert result["plasmid_position"] == 2707  # 2607 + 100


def test_resolve_unknown_attachment(patched_annotation):
    reg = AttachmentRegistry()
    result = asyncio.run(resolve_feature_position(
        {"attachment_id": "att_nope", "feature_name": "Cas9", "kind": "aa_residue", "offset": 10},
        reg,
    ))
    assert result["ok"] is False
    assert "unknown" in result["error"].lower()


def test_resolve_unknown_feature_returns_candidates(patched_annotation):
    reg = _make_registry_with_cas9_at_5111()
    result = asyncio.run(resolve_feature_position(
        {"attachment_id": "att_product_1", "feature_name": "Frobnicator", "kind": "feature_start", "offset": 0},
        reg,
    ))
    assert result["ok"] is False
    assert "no feature matching" in result["error"]
    assert "available_features" in result
    assert "Cas9" in result["available_features"]


def test_rc_helper():
    assert _rc("ATGC") == "GCAT"
    assert _rc("GAT") == "ATC"


def test_resolve_in_make_full_tool_roster():
    from agent_v2 import tools as v2_tools
    roster = v2_tools.make_full_tool_roster()
    names = {t["name"] for t in roster}
    assert "resolve_feature_position" in names


def test_dispatch_chain_routes_resolve_feature_position(patched_annotation):
    reg = _make_registry_with_cas9_at_5111()
    from agent_v2 import tools as v2_tools
    result = asyncio.run(v2_tools.dispatch_with_emitters(
        "resolve_feature_position",
        {"attachment_id": "att_product_1", "feature_name": "Cas9",
         "kind": "aa_residue", "offset": 1},
        reg,
    ))
    assert result["ok"] is True
    assert result["amino_acid"] == "D"
