"""Classifier tests against real .gb fixtures + minimal synthetic inputs.

KEAP1.gb (NCBI RefSeq genomic) and pHAGE_TRE_dCas9_KRAB.gb (SnapGene
synthetic plasmid) are the two canonical examples — they sit at opposite
ends of every heuristic the classifier examines.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401
from agent_v2.file_kind import classify_genbank


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


def test_classify_plasmid_pHAGE():
    gb = (FIXTURES / "pHAGE_TRE_dCas9_KRAB_v034244.gb").read_text()
    fk = classify_genbank(gb)
    assert fk.kind == "plasmid"
    assert fk.confidence >= 0.7
    assert fk.topology == "circular"
    assert "synthetic" in (fk.organism or "").lower()
    assert fk.division == "SYN"


def test_classify_genomic_KEAP1():
    gb = (FIXTURES / "KEAP1.gb").read_text()
    fk = classify_genbank(gb)
    assert fk.kind == "genomic"
    assert fk.confidence >= 0.7
    assert fk.topology == "linear"
    assert fk.organism == "Homo sapiens"
    assert fk.accession.startswith("NC_000019")
    assert fk.has_join_coords is True
    assert fk.has_refseq_keyword is True
    assert fk.n_genes >= 1 and fk.n_mrna >= 1 and fk.n_cds >= 1


def test_classify_unknown_empty():
    fk = classify_genbank("")
    assert fk.kind == "unknown"
    assert fk.confidence == 0.0


def test_classify_synthetic_minimal_plasmid():
    gb = """LOCUS       my_plasmid              3500 bp    DNA     circular SYN 15-MAY-2026
DEFINITION  Test plasmid.
ACCESSION   .
SOURCE      synthetic DNA construct
  ORGANISM  synthetic DNA construct
FEATURES             Location/Qualifiers
     misc_feature    1..100
                     /label="prom"
ORIGIN
        1 acgtacgtac gtacgtacgt
//
"""
    fk = classify_genbank(gb)
    assert fk.kind == "plasmid"
    assert fk.topology == "circular"
    assert "synthetic" in fk.signals[0].lower() or "synthetic" in (fk.organism or "").lower()


def test_classify_minimal_genomic_refseq_join_signal():
    """Even without ORGANISM, join() coords + linear topology + RefSeq KEYWORDS lean genomic."""
    gb = """LOCUS       NG_TEST                 5000 bp    DNA     linear   CON 15-MAY-2026
DEFINITION  Test contig.
ACCESSION   NG_TEST
KEYWORDS    RefSeq.
SOURCE      Homo sapiens
  ORGANISM  Homo sapiens
FEATURES             Location/Qualifiers
     gene            1..5000
                     /gene="TESTGENE"
     mRNA            join(1..100,500..600,1000..2000)
                     /gene="TESTGENE"
                     /transcript_id="NM_TEST.1"
     CDS             join(50..100,500..600,1000..1500)
                     /gene="TESTGENE"
                     /protein_id="NP_TEST.1"
ORIGIN
        1 acgtacgtac
//
"""
    fk = classify_genbank(gb)
    assert fk.kind == "genomic"
    assert fk.has_join_coords is True
    assert fk.n_genes >= 1
