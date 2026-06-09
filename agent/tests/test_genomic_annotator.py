"""Genomic annotator tests on KEAP1.gb + targeted synthetic inputs."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401
from agent_v2.genomic_annotator import annotate_genomic_gb, _looks_like_cds


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


def test_KEAP1_basic_shape():
    gb = (FIXTURES / "KEAP1.gb").read_text()
    ann = annotate_genomic_gb(gb)
    assert ann.organism == "Homo sapiens"
    assert ann.chromosome == "19"
    assert ann.accession.startswith("NC_000019")
    assert ann.length_bp == 20678
    # KEAP1 has 1 gene, 2 mRNA transcript variants, 1 canonical CDS.
    assert "KEAP1" in ann.genes
    assert len(ann.transcripts) == 2
    cds_features = [f for f in ann.features if f.type == "CDS"]
    assert len(cds_features) == 1, f"expected 1 canonical CDS, got {len(cds_features)}"


def test_KEAP1_CDS_translates_to_625aa_starting_MQPDPR():
    gb = (FIXTURES / "KEAP1.gb").read_text()
    ann = annotate_genomic_gb(gb)
    cds = next(f for f in ann.features if f.type == "CDS")
    assert cds.gene == "KEAP1"
    assert cds.strand == -1  # gene is on the minus strand
    assert len(cds.intervals) == 5  # 5 exons in the CDS
    assert cds.protein_id == "NP_036421.2"
    tr = cds.translation or ""
    assert len(tr) == 625, f"KEAP1 protein is 624 residues + stop; got {len(tr)}"
    assert tr.startswith("MQPDPR")
    # 33rd amino acid of KEAP1 is Y (Tyrosine).
    assert tr[32] == "Y", f"expected Y at residue 33; got {tr[32]}"


def test_protein_domain_misc_features_NOT_upgraded():
    """KEAP1's misc_feature annotations (Kelch repeats, Sensor sites) carry
    UniProt xrefs in /note — they must NOT be reclassified as CDS."""
    gb = (FIXTURES / "KEAP1.gb").read_text()
    ann = annotate_genomic_gb(gb)
    upgraded = [f for f in ann.features if f.upgraded_from == "misc_feature"]
    assert upgraded == [], f"no misc_features should upgrade; got {len(upgraded)}"


def test_looks_like_cds_strong_qualifiers():
    """misc_feature with /protein_id, /codon_start, or /translation upgrades."""
    class _F:
        def __init__(self, q): self.qualifiers = q
    assert _looks_like_cds(_F({"protein_id": ["NP_test.1"]}))
    assert _looks_like_cds(_F({"codon_start": ["1"]}))
    assert _looks_like_cds(_F({"translation": ["MKVL..."]}))
    # Note explicitly mentioning standalone "exon" token
    assert _looks_like_cds(_F({"note": ["exon 3 of NM_TEST.1"]}))
    assert _looks_like_cds(_F({"note": ["This is the coding sequence."]}))
    # Just a UniProt xref is NOT enough.
    assert not _looks_like_cds(_F({"db_xref": ["UniProtKB/Swiss-Prot:Q14145.2"]}))
    # Protein-domain note (typical of RefSeq misc_features).
    assert not _looks_like_cds(_F({"note": ["propagated from UniProtKB/Swiss-Prot (Q14145.2); Region: Kelch 1"]}))


def test_dedup_collapses_identical_features():
    """Two identical features with the same coords+gene collapse to one row."""
    gb = """LOCUS       TEST                    1000 bp    DNA     linear   PRI 15-MAY-2026
SOURCE      Homo sapiens
  ORGANISM  Homo sapiens
FEATURES             Location/Qualifiers
     gene            1..900
                     /gene="DUP"
     gene            1..900
                     /gene="DUP"
ORIGIN
        1 acgt
//
"""
    ann = annotate_genomic_gb(gb)
    genes = [f for f in ann.features if f.type == "gene"]
    assert len(genes) == 1, f"two identical /gene features must dedup; got {len(genes)}"
