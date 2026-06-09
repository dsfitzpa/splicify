"""Tests for the kind-aware /annotate-on-upload routing.

Plasmid uploads keep the existing KB-driven pipeline; genomic uploads
go through agent_v2.genomic_annotator and return native GenBank
features only — no Pol2 cassettes, no false-positive cloning features.
"""
import asyncio
import io
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest

import agent_v2  # noqa: F401


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


class _FakeUploadFile:
    """Mimics FastAPI's UploadFile.read() contract."""
    def __init__(self, path: pathlib.Path):
        self.filename = path.name
        self._bytes = path.read_bytes()

    async def read(self) -> bytes:
        return self._bytes


def test_keap1_routes_through_genomic_annotator(monkeypatch):
    """KEAP1.gb classifies as genomic — the plasmid annotate_llm_cached
    must NOT run, and the viz must carry native GenBank features only."""
    plasmid_pipeline_called = {"count": 0}

    async def fake_plasmid_pipeline(seq, *, circular=True, **kw):
        plasmid_pipeline_called["count"] += 1
        return {"annotations": [{"name": "POISONED", "start": 0, "end": 10}]}

    import splicify_api.annotation_cache as ac
    monkeypatch.setattr(ac, "annotate_llm_cached", fake_plasmid_pipeline)

    from agent_v2.router import annotate_on_upload
    res = asyncio.run(annotate_on_upload(_FakeUploadFile(FIXTURES / "KEAP1.gb")))

    assert res["ok"] is True
    assert res["kind"] == "genomic"
    assert res["kind_confidence"] >= 0.7
    assert res["viz"]["type"] == "genomic"
    assert res["viz"]["circular"] is False
    # Plasmid KB pipeline must NOT have run.
    assert plasmid_pipeline_called["count"] == 0
    names = {a["name"] for a in res["viz"]["annotations"]}
    # The poisoned name from the fake plasmid pipeline is absent.
    assert "POISONED" not in names
    # KEAP1 gene appears as a native annotation.
    assert "KEAP1" in names
    # genomic_summary block carries metadata for the frontend.
    gs = res["viz"]["genomic_summary"]
    assert gs["organism"] == "Homo sapiens"
    assert gs["chromosome"] == "19"
    assert gs["n_genes"] >= 1
    assert "KEAP1" in gs["genes"]


def test_keap1_features_have_strand_and_intervals():
    """Each genomic feature row exposes start/end/strand the SeqViz viewer reads."""
    from agent_v2.router import annotate_on_upload
    res = asyncio.run(annotate_on_upload(_FakeUploadFile(FIXTURES / "KEAP1.gb")))
    cds_rows = [a for a in res["viz"]["annotations"] if a["type"] == "CDS"]
    assert cds_rows, "expected at least one CDS row"
    for row in cds_rows:
        assert isinstance(row["start"], int)
        assert isinstance(row["end"], int)
        assert row["direction"] in (1, -1)
        assert row["source"] == "genbank_native"
    # Multi-exon CDS exploded into one row per exon segment.
    keap1_cds_rows = [a for a in cds_rows if "KEAP1" in (a.get("name") or "")]
    assert len(keap1_cds_rows) >= 2


def test_plasmid_keeps_existing_pipeline(monkeypatch):
    """pHAGE_TRE_dCas9_KRAB.gb stays on annotate_llm_cached (plasmid path)."""
    plasmid_pipeline_called = {"count": 0}

    async def fake_plasmid_pipeline(seq, *, circular=True, **kw):
        plasmid_pipeline_called["count"] += 1
        return {"annotations": [{"name": "AmpR", "start": 7712, "end": 8574}],
                 "modules": [], "cloning_features": [],
                 "interactions": [], "hierarchical_annotations": []}

    import splicify_api.annotation_cache as ac
    monkeypatch.setattr(ac, "annotate_llm_cached", fake_plasmid_pipeline)

    from agent_v2.router import annotate_on_upload
    res = asyncio.run(annotate_on_upload(
        _FakeUploadFile(FIXTURES / "pHAGE_TRE_dCas9_KRAB_v034244.gb"),
    ))

    assert res["ok"] is True
    assert res["kind"] == "plasmid"
    assert res["viz"]["type"] == "plasmid"
    assert res["viz"]["circular"] is True
    # Plasmid KB pipeline ran exactly once.
    assert plasmid_pipeline_called["count"] == 1
    # Plasmid path does NOT carry the genomic_summary block.
    assert "genomic_summary" not in res["viz"]


def test_genomic_response_includes_kind_signals():
    """Frontend should be able to surface why the file was classified that way."""
    from agent_v2.router import annotate_on_upload
    res = asyncio.run(annotate_on_upload(_FakeUploadFile(FIXTURES / "KEAP1.gb")))
    sigs = res["kind_signals"]
    assert any("Homo sapiens" in s for s in sigs)
    assert any("RefSeq" in s for s in sigs)
    assert any("join()" in s for s in sigs)
