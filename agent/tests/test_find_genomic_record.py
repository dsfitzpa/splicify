"""Tests for find_genomic_record_tool — NCBI Gene -> .gb -> registered as
genomic attachment. v1 helpers (search_ncbi_gene / fetch_ncbi_genbank)
are monkeypatched so the suite stays offline + deterministic."""
import asyncio
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest

import agent_v2  # noqa: F401
from agent_v2 import attachment_kinds, tools as v2_tools
from agent_v2.crispr_tools import find_genomic_record_tool
from splicify_api.agent.agent_tools import AttachmentRegistry


# A tiny synthetic GenBank record that classifies as genomic + parses
# cleanly via Biopython. 200 bp + a single CDS so the genomic_annotator
# can extract a translation.
_FAKE_GB = """LOCUS       NG_FAKE                  200 bp    DNA     linear   PRI 18-MAY-2026
DEFINITION  Homo sapiens fake gene (FAKE), RefSeqGene.
ACCESSION   NG_FAKE
VERSION     NG_FAKE.1
KEYWORDS    RefSeq.
SOURCE      Homo sapiens (human)
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata.
FEATURES             Location/Qualifiers
     source          1..200
                     /organism="Homo sapiens"
                     /mol_type="genomic DNA"
                     /chromosome="1"
     gene            10..198
                     /gene="FAKE"
                     /db_xref="GeneID:9999"
     mRNA            10..198
                     /gene="FAKE"
                     /transcript_id="NM_FAKE.1"
     CDS             10..198
                     /gene="FAKE"
                     /protein_id="NP_FAKE.1"
                     /codon_start=1
                     /translation="MQPDPRPSGAGACCRFLPLQSQCPEGAGDAVMYASTECKAE"
ORIGIN
        1 ggggggggga tgcagccaga tcctaggcca tctggagcag gagcgtgctg ccgctttctg
       61 ccgttgcaga gccagtgccc cgagggcgca ggcgatgctg tgatgtatgc cagcactgag
      121 tgcaaggcag aatgaaacga cgtcgtactt agtttgcatt gtacgactgg tgcttgttat
      181 gccatgcctt ggccaaatt
//
"""


@pytest.fixture(autouse=True)
def _reset_kinds():
    attachment_kinds.clear_all_kinds()
    yield
    attachment_kinds.clear_all_kinds()


def _make_search_hit(db_source="refseqgene"):
    from splicify_api.external_search import NCBIGeneHit
    return NCBIGeneHit(
        gene_id="9999", accession="NG_FAKE.1",
        db_source=db_source, organism="Homo sapiens",
        title="Homo sapiens fake gene (FAKE), RefSeqGene",
    )


def _patch_ncbi(monkeypatch, *, hit=None, gb_text=_FAKE_GB):
    import splicify_api.external_search as es
    async def _fake_search(symbol, organism="Homo sapiens", *, client=None):
        return hit
    async def _fake_fetch(accession, *, client=None, seq_start=None, seq_stop=None):
        return gb_text
    monkeypatch.setattr(es, "search_ncbi_gene", _fake_search)
    monkeypatch.setattr(es, "fetch_ncbi_genbank", _fake_fetch)


# ---------------------------------------------------------------------------
# happy paths
# ---------------------------------------------------------------------------
def test_find_genomic_record_registers_refseqgene(monkeypatch):
    """User asks for CGAS -> tool downloads NG_*, registers as genomic."""
    _patch_ncbi(monkeypatch, hit=_make_search_hit("refseqgene"))
    reg = AttachmentRegistry()

    out = asyncio.run(find_genomic_record_tool(
        {"gene_symbol": "FAKE", "organism": "Homo sapiens"}, reg,
    ))

    assert out["ok"] is True
    assert out["source"] == "ncbi"
    assert out["db_source"] == "refseqgene"
    assert out["accession"] == "NG_FAKE.1"
    assert out["gene_symbol"] == "FAKE"
    assert out["organism"] == "Homo sapiens"
    assert out["ncbi_gene_id"] == "9999"
    assert 195 <= out["length_bp"] <= 205   # synthetic test record is approx 200 bp

    # Registry got the new attachment as a non-circular genomic record.
    aid = out["attachment_id"]
    att = reg.get(aid)
    assert att is not None
    assert att.circular is False
    assert 195 <= len(att.sequence) <= 205

    # attachment_kinds stashed kind + gb_text.
    fk = attachment_kinds.get_kind(aid)
    assert fk is not None
    assert fk.kind == "genomic"
    assert out["kind"] == "genomic"

    # genomic_annotator findings surfaced in the digest.
    assert "FAKE" in out["genes"]
    assert "NM_FAKE.1" in out["transcripts"]
    assert out["n_cds"] >= 1
    assert out["primary_cds_aa"] >= 50


def test_find_genomic_record_mrna_fallback(monkeypatch):
    """When RefSeqGene isn't available the search returns an mRNA hit and
    the response flags db_source='mrna' + warns about no introns."""
    _patch_ncbi(monkeypatch, hit=_make_search_hit("mrna"))
    reg = AttachmentRegistry()
    out = asyncio.run(find_genomic_record_tool({"gene_symbol": "FAKE"}, reg))
    assert out["ok"] is True
    assert out["db_source"] == "mrna"
    assert "mRNA fallback" in out["note"]
    assert "intronic" in out["note"]


# ---------------------------------------------------------------------------
# error paths
# ---------------------------------------------------------------------------
def test_find_genomic_record_no_hit(monkeypatch):
    _patch_ncbi(monkeypatch, hit=None)
    out = asyncio.run(find_genomic_record_tool({"gene_symbol": "UNOBTANIUM"}, AttachmentRegistry()))
    assert out["ok"] is False
    assert "no NCBI gene record" in out["error"]


def test_find_genomic_record_gb_download_fails(monkeypatch):
    _patch_ncbi(monkeypatch, hit=_make_search_hit(), gb_text=None)
    out = asyncio.run(find_genomic_record_tool({"gene_symbol": "FAKE"}, AttachmentRegistry()))
    assert out["ok"] is False
    assert ".gb download failed" in out["error"]
    assert out["accession"] == "NG_FAKE.1"


def test_find_genomic_record_empty_symbol():
    out = asyncio.run(find_genomic_record_tool({"gene_symbol": ""}, AttachmentRegistry()))
    assert out["ok"] is False
    assert "required" in out["error"]


# ---------------------------------------------------------------------------
# Tool roster + dispatch
# ---------------------------------------------------------------------------
def test_find_genomic_record_in_full_tool_roster():
    roster = v2_tools.make_full_tool_roster()
    names = {t["name"] for t in roster}
    assert "find_genomic_record" in names


def test_dispatch_routes_find_genomic_record(monkeypatch):
    _patch_ncbi(monkeypatch, hit=_make_search_hit())
    reg = AttachmentRegistry()
    out = asyncio.run(v2_tools.dispatch_with_emitters(
        "find_genomic_record", {"gene_symbol": "FAKE"}, reg,
    ))
    assert out["ok"] is True
    assert out["source"] == "ncbi"


def test_target_locator_tool_roster_includes_find_genomic_record():
    from agent_v2.subagents.target_locator import TARGET_LOCATOR_TOOLS
    names = {t["name"] for t in TARGET_LOCATOR_TOOLS}
    assert "find_genomic_record" in names
    assert "annotate_attachment" in names
    assert "resolve_feature_position" in names


def test_find_genomic_record_chromosomal_slice_routing(monkeypatch):
    """db_source='chromosomal_slice' wires through correctly: the helper is
    called with seq_start/seq_stop and the response surfaces the slice
    spatial context fields."""
    from splicify_api.external_search import NCBIGeneHit
    import splicify_api.external_search as es

    captured_fetch = {}

    async def _fake_search(symbol, organism="Homo sapiens", *, client=None):
        return NCBIGeneHit(
            gene_id="115004",
            accession="NC_000006.12",
            db_source="chromosomal_slice",
            organism="Homo sapiens",
            title="Homo sapiens chromosome 6, GRCh38.p14",
            seq_start=73421711, seq_stop=73454297, flanking_bp=2000,
            gene_chr_start=73423711, gene_chr_stop=73452297,
        )

    async def _fake_fetch(accession, *, client=None, seq_start=None, seq_stop=None):
        captured_fetch["accession"] = accession
        captured_fetch["seq_start"] = seq_start
        captured_fetch["seq_stop"] = seq_stop
        return _FAKE_GB

    monkeypatch.setattr(es, "search_ncbi_gene", _fake_search)
    monkeypatch.setattr(es, "fetch_ncbi_genbank", _fake_fetch)

    reg = AttachmentRegistry()
    out = asyncio.run(find_genomic_record_tool(
        {"gene_symbol": "CGAS", "organism": "Homo sapiens"}, reg,
    ))
    assert out["ok"] is True
    assert out["db_source"] == "chromosomal_slice"
    assert out["chromosome_accession"] == "NC_000006.12"
    assert out["slice_seq_start"] == 73421711
    assert out["slice_seq_stop"] == 73454297
    assert out["flanking_bp"] == 2000
    assert out["gene_chr_start"] == 73423711
    assert out["gene_chr_stop"] == 73452297
    assert "Chromosomal slice" in out["note"]

    # fetch_ncbi_genbank received the slice coords (not a full-chromosome fetch).
    assert captured_fetch["accession"] == "NC_000006.12"
    assert captured_fetch["seq_start"] == 73421711
    assert captured_fetch["seq_stop"] == 73454297
