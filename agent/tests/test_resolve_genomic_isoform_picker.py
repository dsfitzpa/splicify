"""Tests for the canonical-isoform picker in _resolve_genomic.

CGAS in the human genomic record has TWO CDS isoforms:
  NP_612450     / 522 aa (canonical)
  NP_001397840 / 497 aa (alternative splice)

Before this fix, `feature_name='CGAS'` returned cds_matches[0] —
whichever appeared first in the GenBank features — and residues 511 /
527 / 530 silently fell out of range on the shorter isoform. After:
the LONGEST translation wins by default; alternative_isoforms surfaces
the others; the LLM can override by passing feature_name=<protein_id>.
"""
import asyncio
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest

import agent_v2  # noqa
from agent_v2 import attachment_kinds
from agent_v2.feature_resolver import _resolve_genomic
from agent_v2.genomic_annotator import GenomicAnnotation, GenomicFeature
from splicify_api.agent.agent_tools import AttachmentRegistry


@pytest.fixture(autouse=True)
def _reset():
    attachment_kinds.clear_all_kinds()
    yield
    attachment_kinds.clear_all_kinds()


def _make_two_isoform_annotation():
    """Synthetic GenomicAnnotation modelling CGAS's two-isoform layout:
    a SHORTER CDS appearing first in the features (NP_001397840, 497 aa),
    and the canonical LONGER CDS second (NP_612450, 522 aa)."""
    short_cds = GenomicFeature(
        type="CDS", gene="CGAS",
        transcript_id="NM_001410911.1", protein_id="NP_001397840.1",
        strand=1, start=2000, end=3500,
        intervals=[(2000, 3500)],
        translation="M" + "A" * 496,   # 497 aa
        qualifiers={}, note=None, label=None,
    )
    long_cds = GenomicFeature(
        type="CDS", gene="CGAS",
        transcript_id="NM_138441.4", protein_id="NP_612450.2",
        strand=1, start=2000, end=3700,
        intervals=[(2000, 3700)],
        translation="M" + "A" * 521,   # 522 aa
        qualifiers={}, note=None, label=None,
    )
    return GenomicAnnotation(
        organism="Homo sapiens", accession="NC_000006.12", chromosome="6",
        length_bp=4000,
        features=[short_cds, long_cds],   # short appears FIRST in features
        transcripts={
            "NM_001410911.1": {"gene": "CGAS"},
            "NM_138441.4": {"gene": "CGAS"},
        },
        genes={"CGAS": {"transcripts": ["NM_001410911.1", "NM_138441.4"]}},
    )


def _register_with_annotation(ann, *, name="CGAS_slice"):
    reg = AttachmentRegistry()
    aid = reg.register_product(name, "A" * ann.length_bp, circular=False)
    from agent_v2.file_kind import FileKind
    fk = FileKind(kind="genomic", topology="linear", organism="Homo sapiens",
                   confidence=0.9, signals=["test"])
    attachment_kinds.stash_kind(aid, fk, gb_text="LOCUS test 4000 bp\nORIGIN\n        1 a\n//\n")
    # Bypass lazy build — inject the prepared annotation directly.
    with attachment_kinds._LOCK:
        attachment_kinds._CACHE[aid]._genomic_annotation = ann
    return reg, aid


# ---------------------------------------------------------------------------
# Core: longest-isoform wins on a bare gene-symbol query
# ---------------------------------------------------------------------------
def test_resolves_to_canonical_longest_isoform_by_default():
    """feature_name='CGAS' must pick the 522-aa isoform, not the 497-aa one
    that happens to appear first in ann.features."""
    ann = _make_two_isoform_annotation()
    reg, aid = _register_with_annotation(ann)
    out = asyncio.run(_resolve_genomic({
        "attachment_id": aid, "feature_name": "CGAS",
        "kind": "aa_residue", "offset": 511,
    }, reg))
    assert out["ok"] is True
    assert out["protein_id"] == "NP_612450.2"
    assert out["transcript_id"] == "NM_138441.4"
    assert out["cds_length_aa"] == 522
    # Residue 511 falls inside the 522-aa isoform.
    assert out["amino_acid"] == "A"


def test_isoform_pick_reason_explains_choice():
    ann = _make_two_isoform_annotation()
    reg, aid = _register_with_annotation(ann)
    out = asyncio.run(_resolve_genomic({
        "attachment_id": aid, "feature_name": "CGAS",
        "kind": "aa_residue", "offset": 1,
    }, reg))
    assert out["ok"] is True
    assert "longest" in out["isoform_pick_reason"].lower()
    assert "2 CDS isoforms" in out["isoform_pick_reason"]


def test_alternative_isoforms_array_surfaces_every_match():
    ann = _make_two_isoform_annotation()
    reg, aid = _register_with_annotation(ann)
    out = asyncio.run(_resolve_genomic({
        "attachment_id": aid, "feature_name": "CGAS",
        "kind": "aa_residue", "offset": 1,
    }, reg))
    alt = out["alternative_isoforms"]
    assert len(alt) == 2
    # Sorted longest-first.
    assert alt[0]["protein_id"] == "NP_612450.2"
    assert alt[0]["length_aa"] == 522
    assert alt[1]["protein_id"] == "NP_001397840.1"
    assert alt[1]["length_aa"] == 497


# ---------------------------------------------------------------------------
# Exact protein_id / transcript_id passed by the LLM wins over longest
# ---------------------------------------------------------------------------
def test_exact_protein_id_match_overrides_longest_preference():
    """When the LLM passes feature_name=<protein_id>, that exact match
    wins even if a longer alternative exists."""
    ann = _make_two_isoform_annotation()
    reg, aid = _register_with_annotation(ann)
    out = asyncio.run(_resolve_genomic({
        "attachment_id": aid, "feature_name": "NP_001397840.1",
        "kind": "aa_residue", "offset": 1,
    }, reg))
    assert out["ok"] is True
    assert out["protein_id"] == "NP_001397840.1"
    assert out["cds_length_aa"] == 497
    assert "exact" in out["isoform_pick_reason"].lower()


def test_exact_transcript_id_match_overrides_longest_preference():
    ann = _make_two_isoform_annotation()
    reg, aid = _register_with_annotation(ann)
    out = asyncio.run(_resolve_genomic({
        "attachment_id": aid, "feature_name": "NM_001410911.1",
        "kind": "aa_residue", "offset": 1,
    }, reg))
    assert out["ok"] is True
    assert out["transcript_id"] == "NM_001410911.1"
    assert out["cds_length_aa"] == 497


# ---------------------------------------------------------------------------
# Out-of-range residue is still rejected — but the message now includes
# cds_length_aa so the LLM can see why.
# ---------------------------------------------------------------------------
def test_out_of_range_residue_returns_range_aware_error():
    ann = _make_two_isoform_annotation()
    reg, aid = _register_with_annotation(ann)
    out = asyncio.run(_resolve_genomic({
        "attachment_id": aid, "feature_name": "NP_001397840.1",
        "kind": "aa_residue", "offset": 511,   # > 497
    }, reg))
    assert out["ok"] is False
    assert "out of range" in out["error"]
    assert "1..497" in out["error"]
    # cds_length_aa is surfaced on the error envelope so the LLM can pivot.
    assert out["cds_length_aa"] == 497
