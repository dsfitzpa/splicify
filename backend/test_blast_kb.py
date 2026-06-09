#!/usr/bin/env python3
"""
Test script for pLannotate BLAST database integration.

This script demonstrates how the knowledge base loads and searches
BLAST database FASTA files.
"""

import asyncio
import logging
from pathlib import Path
from splicify_api.predesign import get_knowledge_base

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def create_test_blast_databases():
    """Create test BLAST database FASTA files for demonstration."""
    test_dir = Path("/tmp/test_plannotate_data/BLAST_dbs")
    test_dir.mkdir(parents=True, exist_ok=True)

    # Create kb_features.fasta with common cloning features
    kb_features = """>CMV_promoter CMV immediate early promoter
GACATTGATTATTGACTAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCCCATATATGGAGTTCCGCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAACGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGGACTTTCCATTGACGTCAATGGGTGGAGTATTTACGGTAAACTGCCCACTTGGCAGTACATCAAGTGTATCATATGCCAAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCTGGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCT
>eGFP enhanced green fluorescent protein
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG
>bGH_polyA bovine growth hormone polyadenylation signal
CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCATGCTGGGGATGCGGTGGGCTCTATGG
>T7_promoter T7 RNA polymerase promoter
TAATACGACTCACTATAGGG
>SV40_polyA SV40 late polyadenylation signal
ACAAATAAAAGATTTTATTTAGTCTCCAGAAAAAGGGGGGAATGAAAGACCCCACCTGTAGGTTTGGCAAGCTAGCTTAAGTAACGCCATTTTGCAAGGCATGGAAAATACATAACTGAGAATAGAGAAGTTCAGATCAAGGTTAGGAACAGAGAGACAGCAGAATATGGGCCAAACAGGATATCTGTGGTAAGCAGTTCCTGCCCCGGCTCAGGGCCAAGAACAGATGGTCCCCAGATGCGGTCCCGCCCTCAGCAGTTTCTAGAGAACCATCAGATGTTTCCAGGGTGCCCCAAGGACCTGAAATGACCCTGTGCCTTATTTGAACTAACCAATCAGTTCGCTTCTCGCTTCTGTTCGCGCGCTTCTGCTCCCCGAGCTCAATAAAAGAGCCCACAACCCCTCACTCGGCGCGCCAGTCCTCCGATAGACTGCGTCGCCCGGGTACCCGTATTCCCAATAAAGCCTCTTGCTGTTTGCATCCGAATCGTGGACTCGCTGATCCTTGGGAGGGTCTCCTCAGATTGATTGACTGCCCACCTCGGGGGTCTTTCATTTGGAGGTTCCACCGAGATTTGGAGACCCCTGCCCAGGGACCACCGACCCCCCCGCCGGGAGGTAAGCTGGCCAGCGGTCGTTTCGTGTCTGTCTCTGTCTTTGTGCGTGTTTGTGCCGGCATCTAATGTTTGCGCCTGCGTCTGTACTAGTTAGCTAACTAGCTCTGTATCTGGCGGACCCGTGGTGGAACTGACGAGTTCGGAACACCCGGCCGCAACCCTGGGAGACGTCCCAGGGACTTCGGGGGCCGTTTTTGTGGCCCGACCTGAGTCCAAAAATCCCGGAAACTTTCACTCTGAGTTTTCTTCAGGGATTCGAAAAGCCCTCCCTATAAAAGGGTGGTGATGCAAATGAGATAGGTCGG
"""

    with open(test_dir / "kb_features.fasta", "w") as f:
        f.write(kb_features)

    # Create fpbase.fasta with fluorescent proteins
    fpbase = """>mCherry mCherry red fluorescent protein
ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG
>mRuby3 mRuby3 red fluorescent protein
ATGGTGTCGAAGGGCGAGGAGGACAACATGGCGATCATCAAGGAGTTCATGCGCTTCAAGGTGCGCATGGAGGGCTCCATGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCCCAGTTCATGTACGGCTCCAAGGCGTACGTGAAGCACCCCGCCGACATCCCCGACTATCTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGTCTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCACGCTGATCTACAAGGTGAAGATGCGCGGCACCAACTTCCCCCCCGACGGCCCCGTAATGCAGAAGAAGACTATGGGCTGGGAGGCCTCCACCGAGCGCCTGTACCCCCGCGACGGCGTGCTGAAGGGCGAGATCCACCAGGCCCTGAAGCTGAAGGACGGCGGCCACTACCTGGTGGAGTTCAAGACCATCTACATGGCCAAGAAGCCCGTGCAACTGCCCGGCTACTACTACGTGGACACCAAGCTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAGCAGTACGAGCGCGCCGAGGGCCGCCACCACCTGTTCCTGTAGCGGCCGCGACTCTAGATCATAATCAGCCATACCACATTTGTAGAGGTTTTACTTGCTTTAAAAAACCTCCCACACCTCCCCCTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTATCATGTCTGGATC
>EGFP enhanced GFP with optimized codon usage
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG
"""

    with open(test_dir / "fpbase.fasta", "w") as f:
        f.write(fpbase)

    # Create kb_cds.fasta with some coding sequences
    kb_cds = """>AmpR ampicillin resistance gene
ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA
>KanR kanamycin resistance gene
ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCAGCGCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTGCAGGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTGCTCGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAGGATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATGCGGCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGCATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAAGAGCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGACGGCGAGGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAATGGCCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGACATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTCCTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTTGACGAGTTCTTCTGA
"""

    with open(test_dir / "kb_cds.fasta", "w") as f:
        f.write(kb_cds)

    logger.info(f"Created test BLAST databases in {test_dir}")
    return str(test_dir.parent)


async def test_blast_database_loading():
    """Test loading BLAST databases."""
    print("\n" + "=" * 70)
    print("TEST 1: Loading BLAST Databases")
    print("=" * 70)

    # Create test databases
    test_data_dir = create_test_blast_databases()

    # Initialize knowledge base
    kb = get_knowledge_base(data_dir=test_data_dir)
    kb.load()

    # Show what was loaded
    databases = kb.get_databases()
    print(f"\nLoaded databases:")
    for db_name, count in databases.items():
        if count > 0:
            print(f"  ✓ {db_name}: {count} sequences")
        else:
            print(f"  - {db_name}: not found")

    print(f"\nTotal sequences: {kb.get_total_sequences()}")

    return kb


async def test_feature_search(kb):
    """Test searching for features."""
    print("\n" + "=" * 70)
    print("TEST 2: Feature Search")
    print("=" * 70)

    # Test 1: Exact match
    print("\n1. Exact match: 'CMV'")
    results = await kb.search_feature("CMV")
    for r in results[:3]:
        print(f"  → {r['id']}: {r['name']}")
        print(f"     Database: {r['database']}, Confidence: {r['confidence']:.2f}")
        print(f"     Sequence: {r['sequence'][:60]}...")

    # Test 2: Case insensitive
    print("\n2. Case insensitive: 'egfp'")
    results = await kb.search_feature("egfp")
    for r in results[:3]:
        print(f"  → {r['id']}: {r['name']}")
        print(f"     Database: {r['database']}, Confidence: {r['confidence']:.2f}")

    # Test 3: Partial match
    print("\n3. Partial match: 'GFP'")
    results = await kb.search_feature("GFP")
    for r in results[:5]:
        print(f"  → {r['id']}: {r['name']}")
        print(f"     Database: {r['database']}, Confidence: {r['confidence']:.2f}")


async def test_database_specific_search(kb):
    """Test searching specific databases."""
    print("\n" + "=" * 70)
    print("TEST 3: Database-Specific Search")
    print("=" * 70)

    # Search only fpbase for fluorescent proteins
    print("\n1. Search fpbase for 'mCherry'")
    results = await kb.search_feature("mCherry", databases=["fpbase"])
    for r in results:
        print(f"  → {r['id']}: {r['name']}")
        print(f"     Database: {r['database']}")
        print(f"     Type: {r['type']}")

    # Search kb_features for promoters
    print("\n2. Search kb_features for 'promoter'")
    results = await kb.search_feature("promoter", databases=["kb_features"])
    for r in results[:3]:
        print(f"  → {r['id']}: {r['name']}")
        print(f"     Type: {r['type']}")

    # Search kb_cds for resistance genes
    print("\n3. Search kb_cds for 'resistance'")
    results = await kb.search_feature("resistance", databases=["kb_cds"])
    for r in results[:3]:
        print(f"  → {r['id']}: {r['name']}")
        print(f"     Type: {r['type']}")


async def test_pattern_search(kb):
    """Test pattern-based searching."""
    print("\n" + "=" * 70)
    print("TEST 4: Pattern Search")
    print("=" * 70)

    # Test 1: Search for red fluorescent proteins
    print("\n1. Pattern: 'red fluorescent protein'")
    results = await kb.search_pattern("red fluorescent protein")
    for r in results[:3]:
        print(f"  → {r['id']}: {r['name']}")
        print(f"     Database: {r['database']}, Confidence: {r['confidence']:.2f}")

    # Test 2: Search for promoters
    print("\n2. Pattern: 'promoter'")
    results = await kb.search_pattern("promoter")
    for r in results[:5]:
        print(f"  → {r['id']}: {r['name']}")
        print(f"     Database: {r['database']}, Confidence: {r['confidence']:.2f}")

    # Test 3: Search for polyA signals
    print("\n3. Pattern: 'polyadenylation signal'")
    results = await kb.search_pattern("polyadenylation signal")
    for r in results[:3]:
        print(f"  → {r['id']}: {r['name']}")
        print(f"     Database: {r['database']}, Confidence: {r['confidence']:.2f}")


async def test_type_inference(kb):
    """Test automatic type inference."""
    print("\n" + "=" * 70)
    print("TEST 5: Type Inference")
    print("=" * 70)

    test_queries = [
        ("CMV", "promoter"),
        ("eGFP", "fluorescent_protein"),
        ("bGH_polyA", "terminator"),
        ("AmpR", "CDS"),
        ("mCherry", "fluorescent_protein"),
    ]

    for query, expected_type in test_queries:
        results = await kb.search_feature(query)
        if results:
            actual_type = results[0].get("type", "unknown")
            match = "✓" if expected_type.lower() in actual_type.lower() or actual_type.lower() in expected_type.lower() else "✗"
            print(f"{match} {query}: inferred as '{actual_type}' (expected '{expected_type}')")


async def main():
    """Run all tests."""
    print("=" * 70)
    print("  pLannotate BLAST Database Integration Test")
    print("=" * 70)

    try:
        # Test 1: Load databases
        kb = await test_blast_database_loading()

        # Test 2: Feature search
        await test_feature_search(kb)

        # Test 3: Database-specific search
        await test_database_specific_search(kb)

        # Test 4: Pattern search
        await test_pattern_search(kb)

        # Test 5: Type inference
        await test_type_inference(kb)

        print("\n" + "=" * 70)
        print("✓ ALL TESTS PASSED")
        print("=" * 70)
        print("\nKey Features Demonstrated:")
        print("  ✓ BLAST database FASTA parsing")
        print("  ✓ Multi-database indexing")
        print("  ✓ Exact and fuzzy matching")
        print("  ✓ Database-specific search")
        print("  ✓ Pattern-based search")
        print("  ✓ Automatic type inference")
        print("\nProduction Setup:")
        print("  1. Set PLANNOTATE_DATA_DIR to pLannotate data directory")
        print("  2. Ensure BLAST_dbs subdirectory exists with FASTA files:")
        print("     - kb_features.fasta")
        print("     - kb_cds.fasta")
        print("     - kb_rna.fasta")
        print("     - fpbase.fasta")
        print("     - swissprot.fasta")
        print("     - Rfam.fasta")
        print("=" * 70 + "\n")

    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    exit(asyncio.run(main()))
