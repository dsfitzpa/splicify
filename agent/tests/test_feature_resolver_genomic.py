"""Genomic-aware feature_resolver tests.

Covers the new branch in resolve_feature_position that walks joined exon
coords for aa_residue resolution against a KEAP1.gb attachment.
"""
import asyncio
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest

import agent_v2  # noqa: F401
from agent_v2 import attachment_kinds
from agent_v2.feature_resolver import resolve_feature_position, _walk_intervals_to_codon
from agent_v2.file_kind import classify_genbank
from splicify_api.agent.agent_tools import (
    AttachmentRegistry, extract_seq_from_genbank,
)


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


@pytest.fixture(autouse=True)
def _reset_kinds():
    attachment_kinds.clear_all_kinds()
    yield
    attachment_kinds.clear_all_kinds()


def _registered_keap1():
    gb = (FIXTURES / "KEAP1.gb").read_text()
    reg = AttachmentRegistry()
    seq = extract_seq_from_genbank(gb)
    aid = reg.register_product("KEAP1", seq, circular=False)
    fk = classify_genbank(gb)
    attachment_kinds.stash_kind(aid, fk, gb_text=gb)
    return reg, aid


def test_resolve_33rd_aa_of_KEAP1_returns_tyrosine():
    reg, aid = _registered_keap1()
    out = asyncio.run(resolve_feature_position(
        {"attachment_id": aid, "feature_name": "KEAP1",
         "kind": "aa_residue", "offset": 33}, reg,
    ))
    assert out["ok"] is True
    assert out["amino_acid"] == "Y"
    assert out["feature_strand"] == "-"
    assert out["n_exons"] == 5
    assert out["spliced_cds_length_bp"] == 5 * len(out["codon"]) or out["spliced_cds_length_bp"] > 600 * 3
    # Position should fall inside the gene's genomic span.
    assert out["feature_start"] <= out["plasmid_position"] <= out["feature_end"]


def test_resolve_1st_aa_of_KEAP1_returns_M_start_codon():
    reg, aid = _registered_keap1()
    out = asyncio.run(resolve_feature_position(
        {"attachment_id": aid, "feature_name": "KEAP1",
         "kind": "aa_residue", "offset": 1}, reg,
    ))
    assert out["ok"] is True
    assert out["amino_acid"] == "M"
    assert out["codon"] in ("ATG",)  # start codon


def test_resolve_out_of_range_residue():
    reg, aid = _registered_keap1()
    out = asyncio.run(resolve_feature_position(
        {"attachment_id": aid, "feature_name": "KEAP1",
         "kind": "aa_residue", "offset": 9999}, reg,
    ))
    assert out["ok"] is False
    assert "out of range" in out["error"].lower()


def test_resolve_unknown_gene_returns_candidates():
    reg, aid = _registered_keap1()
    out = asyncio.run(resolve_feature_position(
        {"attachment_id": aid, "feature_name": "MYC",
         "kind": "aa_residue", "offset": 10}, reg,
    ))
    assert out["ok"] is False
    assert "no CDS matching" in out["error"]
    assert "KEAP1" in out.get("available_features", [])


def test_resolve_by_protein_id_works():
    """Look up by protein_id (NP_036421.2) instead of gene name."""
    reg, aid = _registered_keap1()
    out = asyncio.run(resolve_feature_position(
        {"attachment_id": aid, "feature_name": "NP_036421.2",
         "kind": "aa_residue", "offset": 33}, reg,
    ))
    assert out["ok"] is True
    assert out["amino_acid"] == "Y"
    assert out["protein_id"] == "NP_036421.2"


def test_resolve_plasmid_kind_uses_legacy_path():
    """When kind is plasmid (no stashed kind), feature_resolver falls back
    to the original v1 annotate_cached path."""
    reg = AttachmentRegistry()
    reg.register_product("test", "ACGT" * 200, circular=True)
    # No stash_kind() call => default plasmid path. annotate_cached will be
    # invoked and (without mocking) will fail; we just need the branch to
    # NOT call the genomic resolver.
    out = asyncio.run(resolve_feature_position(
        {"attachment_id": "att_product_1", "feature_name": "ZZZ",
         "kind": "aa_residue", "offset": 1}, reg,
    ))
    # Either annotation fails or no-feature-match — but we shouldn't see
    # the genomic error string.
    assert "no genomic annotation cached" not in (out.get("error") or "")


def test_walk_intervals_simple_forward_strand():
    """3 bases, all in one segment, + strand."""
    seq = "AAACCCTTT"  # 9 bp
    intervals = [(0, 9)]
    pos, codon, spans, _lo, _hi = _walk_intervals_to_codon(intervals, +1, seq, cumulative_bp_offset=3)
    assert pos == 3
    assert codon == "CCC"
    assert spans is False


def test_walk_intervals_spans_exon_boundary():
    """Codon crosses an exon boundary: should set spans_intron=True."""
    # Exon1 = bases 0..2 (3 bp); Exon2 = bases 10..12 (3 bp). Spliced CDS = 6 bp.
    # Codon offset 1..3 (1-indexed offsets 2..4) crosses the boundary.
    seq = "AAA" + "x" * 7 + "TTT"   # 13 bp
    intervals = [(0, 3), (10, 13)]
    pos, codon, spans, _lo, _hi = _walk_intervals_to_codon(intervals, +1, seq, cumulative_bp_offset=1)
    # Bases at offset 1,2,3 in the spliced CDS = seq[1] + seq[2] + seq[10] = 'A','A','T'
    assert codon == "AAT"
    assert spans is True


def test_walk_intervals_minus_strand_revcomp():
    """Minus strand: codon is reverse-complemented from genome to sense."""
    # 9-base segment on - strand. Bases at genome positions 6,7,8 are
    # the FIRST codon when reading + strand bases 6,7,8, but the codon
    # is read as the reverse-complement on the - strand.
    seq = "AAATTTGGG"  # 9 bp
    intervals = [(0, 9)]
    pos, codon, spans, _lo, _hi = _walk_intervals_to_codon(intervals, -1, seq, cumulative_bp_offset=0)
    # First codon on - strand starts from base 8 (end of interval),
    # reading 8 -> 7 -> 6: bases G,G,G -> revcomp = CCC.
    assert codon == "CCC"
    assert pos == 8
    assert spans is False
