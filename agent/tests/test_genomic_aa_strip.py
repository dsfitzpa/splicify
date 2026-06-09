"""Tests for the per-exon translation annotations the genomic upload emits."""
import asyncio
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest

import agent_v2  # noqa: F401


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


class _Upload:
    def __init__(self, path: pathlib.Path):
        self.filename = path.name
        self._bytes = path.read_bytes()

    async def read(self) -> bytes:
        return self._bytes


def test_keap1_emits_per_exon_translations():
    """KEAP1's 5-exon CDS yields 5 translation annotations covering all 625 aa."""
    from agent_v2.router import annotate_on_upload
    res = asyncio.run(annotate_on_upload(_Upload(FIXTURES / "KEAP1.gb")))
    assert res["ok"] is True and res["kind"] == "genomic"
    trans = [a for a in res["viz"]["annotations"] if a.get("layer") == "translation"]
    assert len(trans) == 5, f"expected 5 exon translations; got {len(trans)}"
    # AA ranges are contiguous + cover 1..625 with no gaps and no overlap.
    sorted_t = sorted(trans, key=lambda a: a["metadata"]["aa_start_global"])
    starts = [t["metadata"]["aa_start_global"] for t in sorted_t]
    ends = [t["metadata"]["aa_end_global"] for t in sorted_t]
    assert starts[0] == 1
    assert ends[-1] == 625
    for i in range(1, len(sorted_t)):
        # Adjacent exon AA-ranges abut (codon-spanning-intron AAs sit on the upstream exon).
        assert starts[i] == ends[i - 1] or starts[i] == ends[i - 1] + 1


def test_translation_annotation_carries_aa_sequence_metadata():
    from agent_v2.router import annotate_on_upload
    res = asyncio.run(annotate_on_upload(_Upload(FIXTURES / "KEAP1.gb")))
    trans = [a for a in res["viz"]["annotations"] if a.get("layer") == "translation"]
    for t in trans:
        m = t["metadata"]
        assert isinstance(m["aa_sequence"], str)
        # AA sequence length matches the annotated [aa_start, aa_end] window
        # (within 1 aa — codon-spanning-intron rounding).
        expected_len = m["aa_end_global"] - m["aa_start_global"] + 1
        assert abs(len(m["aa_sequence"]) - expected_len) <= 1
        assert m["cds_total_aa"] == 625
        assert m["gene"] == "KEAP1"
        assert m["protein_id"] == "NP_036421.2"
        assert m["orf_detected"] is True
        # Frontend AA-strip rendering keys.
        assert t["color"] == "#673AB7"
        assert t["module_type"] == "translation"
        assert t["source"] == "orf_detection"


def test_first_translation_starts_with_M_for_KEAP1():
    """KEAP1's protein starts MQPDPR — the first exon's AA slice must begin with M."""
    from agent_v2.router import annotate_on_upload
    res = asyncio.run(annotate_on_upload(_Upload(FIXTURES / "KEAP1.gb")))
    trans = sorted(
        [a for a in res["viz"]["annotations"] if a.get("layer") == "translation"],
        key=lambda a: a["metadata"]["exon_idx"],
    )
    first_exon = trans[0]
    assert first_exon["metadata"]["aa_sequence"].startswith("MQPDPR")
    assert first_exon["metadata"]["exon_idx"] == 1
    assert first_exon["metadata"]["n_exons"] == 5


def test_translation_emitted_only_for_CDS():
    """gene / mRNA / misc_feature / regulatory rows do NOT get a translation."""
    from agent_v2.router import annotate_on_upload
    res = asyncio.run(annotate_on_upload(_Upload(FIXTURES / "KEAP1.gb")))
    trans = [a for a in res["viz"]["annotations"] if a.get("layer") == "translation"]
    types_with_trans = {t["type"] for t in trans}
    assert types_with_trans == {"translation"}


def test_minus_strand_translation_carries_direction_minus1():
    """KEAP1 is on the minus strand — translations must inherit direction=-1."""
    from agent_v2.router import annotate_on_upload
    res = asyncio.run(annotate_on_upload(_Upload(FIXTURES / "KEAP1.gb")))
    trans = [a for a in res["viz"]["annotations"] if a.get("layer") == "translation"]
    assert all(t["direction"] == -1 for t in trans)
    assert all(t["strand"] == -1 for t in trans)
