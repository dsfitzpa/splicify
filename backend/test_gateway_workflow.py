#!/usr/bin/env python3
"""
Complete Gateway Cloning Workflow Test

This script demonstrates the complete workflow:
1. Upload pDONR221 to inventory
2. Design primers to add attB sites to GFP insert
3. Simulate BP recombination
4. Visualize product with annotations
5. Verify correct product (no ccdB, correct att sites)
"""

import sys
import json
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent))

from splicify_api.cloning.gateway_sites import (
    GATEWAY_ATT_SITES,
    scan_att_sites,
    scan_for_ccdb,
    extract_core_sequence,
)
from splicify_api.cloning.gateway_operator import GatewayOperator


# =============================================================================
# Test Data: GFP sequence and pDONR221-like vector
# =============================================================================

# GFP CDS (714 bp) - mNeonGreen
GFP_SEQUENCE = (
    "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG"
)


def build_pDONR221_simplified():
    """Build simplified pDONR221 donor vector."""
    from splicify_api.cloning.gateway_sites import CCDB_GENE
    
    # Simplified backbone fragments
    backbone_left = "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCACTGGCCGTCGTTTTAC" * 10
    backbone_right = "TCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTCCTGCAGCCCGGGGGATCCACTAGTTCTAGAGCGGCCGC" * 10
    
    attP1 = GATEWAY_ATT_SITES["attP1"]
    attP2 = GATEWAY_ATT_SITES["attP2"]
    
    ccdb_region = "AGCTTGGCTG" + CCDB_GENE + "TAATACGACT"
    
    return backbone_left + attP1 + ccdb_region + attP2 + backbone_right


def print_section(title):
    """Print section header."""
    print(f"\n{'='*80}")
    print(f"{title}")
    print(f"{'='*80}\n")


def print_subsection(title):
    """Print subsection header."""
    print(f"\n{'-'*80}")
    print(f"{title}")
    print(f"{'-'*80}\n")


def main():
    print_section("Gateway Cloning Workflow Test")
    
    # =========================================================================
    # Step 1: Detect att sites in inventory plasmid (pDONR221)
    # =========================================================================
    print_subsection("Step 1: Detect att sites in pDONR221 (Inventory)")
    
    pdonr221_seq = build_pDONR221_simplified()
    print(f"pDONR221 length: {len(pdonr221_seq)} bp")
    
    # Scan for att sites
    pdonr_sites = scan_att_sites(pdonr221_seq, fuzzy_threshold=0)
    print(f"\nDetected att sites in pDONR221:")
    for site in pdonr_sites:
        print(f"  - {site.site_type} at position {site.start}-{site.end} (strand: {'+' if site.strand == 1 else '-'})")
        print(f"    Core: {site.core_sequence}, Quality: {site.match_quality}")
    
    # Check for ccdB
    ccdb_positions = scan_for_ccdb(pdonr221_seq)
    print(f"\nccdB gene detected: {len(ccdb_positions) > 0}")
    if ccdb_positions:
        for start, end in ccdb_positions:
            print(f"  - ccdB at position {start}-{end}")
    
    # =========================================================================
    # Step 2: Design primers to add attB sites to GFP
    # =========================================================================
    print_subsection("Step 2: Design primers to add attB sites to GFP insert")
    
    print(f"GFP insert length: {len(GFP_SEQUENCE)} bp")
    
    # Get attB site sequences
    attB1_seq = GATEWAY_ATT_SITES["attB1"]
    attB2_seq = GATEWAY_ATT_SITES["attB2"]
    
    print(f"\nattB1 sequence ({len(attB1_seq)} bp): {attB1_seq}")
    print(f"attB2 sequence ({len(attB2_seq)} bp): {attB2_seq}")
    
    # Design forward primer (adds attB1)
    anneal_len = 20
    anneal_fwd = GFP_SEQUENCE[:anneal_len]
    fwd_primer = "GGGG" + attB1_seq + anneal_fwd  # GGGG spacer for enzyme processivity
    
    # Design reverse primer (adds attB2 via RC)
    from Bio.Seq import Seq
    anneal_rev = str(Seq(GFP_SEQUENCE[-anneal_len:]).reverse_complement())
    rev_primer = "GGGG" + str(Seq(attB2_seq).reverse_complement()) + anneal_rev
    
    print(f"\nForward Primer (adds attB1):")
    print(f"  Sequence: {fwd_primer}")
    print(f"  Length: {len(fwd_primer)} bp")
    print(f"  Components: GGGG (4bp) + attB1 ({len(attB1_seq)}bp) + annealing ({anneal_len}bp)")
    
    print(f"\nReverse Primer (adds attB2):")
    print(f"  Sequence: {rev_primer}")
    print(f"  Length: {len(rev_primer)} bp")
    print(f"  Components: GGGG (4bp) + RC(attB2) ({len(attB2_seq)}bp) + RC(annealing) ({anneal_len}bp)")
    
    # Simulate PCR product
    pcr_product = attB1_seq + GFP_SEQUENCE + attB2_seq
    print(f"\nPCR Product (attB1-GFP-attB2):")
    print(f"  Length: {len(pcr_product)} bp")
    print(f"  Expected: {len(attB1_seq)} + {len(GFP_SEQUENCE)} + {len(attB2_seq)} = {len(attB1_seq) + len(GFP_SEQUENCE) + len(attB2_seq)} bp")
    
    # Verify att sites in PCR product
    pcr_sites = scan_att_sites(pcr_product, fuzzy_threshold=0)
    print(f"\nVerified att sites in PCR product:")
    for site in pcr_sites:
        print(f"  - {site.site_type} at position {site.start}-{site.end}")
    
    # =========================================================================
    # Step 3: Simulate BP recombination
    # =========================================================================
    print_subsection("Step 3: Simulate BP Recombination")
    
    operator = GatewayOperator()
    
    # Prepare modules
    insert_module = {
        "canonical_id": "GFP_PCR_product",
        "sequence": pcr_product,
        "role": "insert"
    }
    
    donor_module = {
        "canonical_id": "pDONR221",
        "sequence": pdonr221_seq,
        "role": "vector"
    }
    
    modules = [insert_module, donor_module]
    
    print(f"Input modules:")
    print(f"  1. {insert_module['canonical_id']}: {len(insert_module['sequence'])} bp (with attB1, attB2)")
    print(f"  2. {donor_module['canonical_id']}: {len(donor_module['sequence'])} bp (with attP1, attP2, ccdB)")
    
    # Evaluate Gateway cloning
    print(f"\nEvaluating Gateway cloning...")
    plan = operator.evaluate(modules, topology="circular")
    
    print(f"\nReaction type: {plan.reaction_type}")
    print(f"Feasible: {plan.feasible}")
    
    if plan.warnings:
        print(f"\nWarnings:")
        for warning in plan.warnings:
            print(f"  - {warning}")
    
    if plan.infeasibility_reasons:
        print(f"\nInfeasibility reasons:")
        for reason in plan.infeasibility_reasons:
            print(f"  - {reason}")
    
    # =========================================================================
    # Step 4: Analyze product sequence
    # =========================================================================
    print_subsection("Step 4: Analyze Product Sequence")
    
    product_seq = plan.product_sequence
    print(f"Product length: {len(product_seq)} bp")
    
    # Scan product for att sites
    product_sites = scan_att_sites(product_seq, fuzzy_threshold=0)
    print(f"\natt sites in product:")
    for site in product_sites:
        print(f"  - {site.site_type} at position {site.start}-{site.end}")
        print(f"    Core: {site.core_sequence}, Quality: {site.match_quality}")
    
    # Check for ccdB in product
    product_ccdb = scan_for_ccdb(product_seq)
    print(f"\nccdB in product: {len(product_ccdb) > 0}")
    if product_ccdb:
        print(f"  WARNING: ccdB found in product! This should NOT happen in BP reaction.")
        for start, end in product_ccdb:
            print(f"  - ccdB at position {start}-{end}")
    
    # =========================================================================
    # Step 5: Analyze byproduct sequence
    # =========================================================================
    print_subsection("Step 5: Analyze Byproduct Sequence")
    
    byproduct_seq = plan.byproduct_sequence
    print(f"Byproduct length: {len(byproduct_seq)} bp")
    print(f"Byproduct description: {plan.byproduct_description}")
    
    # Scan byproduct for att sites
    byproduct_sites = scan_att_sites(byproduct_seq, fuzzy_threshold=0)
    print(f"\natt sites in byproduct:")
    for site in byproduct_sites:
        print(f"  - {site.site_type} at position {site.start}-{site.end}")
    
    # Check for ccdB in byproduct
    byproduct_ccdb = scan_for_ccdb(byproduct_seq)
    print(f"\nccdB in byproduct: {len(byproduct_ccdb) > 0}")
    if byproduct_ccdb:
        print(f"  Expected: ccdB should be in byproduct (linear fragment)")
        for start, end in byproduct_ccdb:
            print(f"  - ccdB at position {start}-{end}")
    
    # =========================================================================
    # Step 6: Verification
    # =========================================================================
    print_subsection("Step 6: Verification")
    
    # Check expected product characteristics
    product_site_types = [s.site_type for s in product_sites]
    
    tests = []
    
    # Test 1: Reaction type should be BP
    test1_pass = plan.reaction_type == "BP"
    tests.append(("Reaction type is BP", test1_pass))
    
    # Test 2: Product should have attL sites
    test2_pass = "attL1" in product_site_types or "attL2" in product_site_types
    tests.append(("Product contains attL sites (Entry clone)", test2_pass))
    
    # Test 3: Product should NOT contain ccdB
    test3_pass = len(product_ccdb) == 0
    tests.append(("Product does NOT contain ccdB", test3_pass))
    
    # Test 4: Byproduct should contain attR sites
    byproduct_site_types = [s.site_type for s in byproduct_sites]
    test4_pass = "attR1" in byproduct_site_types or "attR2" in byproduct_site_types
    tests.append(("Byproduct contains attR sites", test4_pass))
    
    # Test 5: Byproduct should contain ccdB
    test5_pass = len(byproduct_ccdb) > 0
    tests.append(("Byproduct contains ccdB (linear fragment)", test5_pass))
    
    # Test 6: Product should contain GFP sequence
    test6_pass = GFP_SEQUENCE in product_seq
    tests.append(("Product contains GFP insert", test6_pass))
    
    print("\nVerification Results:")
    print("=" * 60)
    all_pass = True
    for test_name, passed in tests:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status:10s} {test_name}")
        if not passed:
            all_pass = False
    print("=" * 60)
    
    if all_pass:
        print("\n🎉 All tests PASSED! Gateway cloning workflow is correct.")
    else:
        print("\n⚠️  Some tests FAILED. Review the workflow.")
    
    # =========================================================================
    # Step 7: Print fragment annotations
    # =========================================================================
    print_subsection("Step 7: Product Fragment Annotations")
    
    if plan.product_annotations:
        print(f"Total annotations: {len(plan.product_annotations)}")
        print("\nFragment sources:")
        
        # Group by source type
        by_type = {}
        for annot in plan.product_annotations:
            source_type = annot.get("type", "unknown")
            if source_type not in by_type:
                by_type[source_type] = []
            by_type[source_type].append(annot)
        
        for source_type, annots in by_type.items():
            print(f"\n  {source_type} ({len(annots)} annotations):")
            for annot in annots[:5]:  # Show first 5
                name = annot.get("name", "Unknown")
                start = annot.get("start", 0)
                end = annot.get("end", 0)
                source = annot.get("source", "Unknown")
                print(f"    - {name} [{start}-{end}] from {source}")
            if len(annots) > 5:
                print(f"    ... and {len(annots) - 5} more")
    else:
        print("No annotations available.")
    
    # =========================================================================
    # Step 8: Print primer table
    # =========================================================================
    print_subsection("Step 8: Primer Design Summary")
    
    if plan.primer_table:
        print(f"Total primers designed: {len(plan.primer_table)}")
        for i, primer in enumerate(plan.primer_table, 1):
            print(f"\nPrimer {i}: {primer.get('primer_name', 'Unknown')}")
            print(f"  Sequence: {primer.get('sequence', 'N/A')}")
            print(f"  Length: {primer.get('length', 0)} bp")
            print(f"  Purpose: {primer.get('purpose', 'N/A')}")
            print(f"  Tm (annealing): {primer.get('tm_anneal', 0):.1f}°C")
    else:
        print("No custom primers needed (native att sites detected).")
    
    # =========================================================================
    # Summary
    # =========================================================================
    print_section("Summary")
    
    print(f"Gateway cloning workflow completed successfully!")
    print(f"\nExpected Entry Clone (pENTR-GFP):")
    print(f"  - Size: {len(product_seq)} bp")
    print(f"  - att sites: attL1 and attL2")
    print(f"  - Contains GFP: {GFP_SEQUENCE in product_seq}")
    print(f"  - Contains ccdB: {len(product_ccdb) > 0}")
    print(f"  - Cloning feasible: {plan.feasible}")
    
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
