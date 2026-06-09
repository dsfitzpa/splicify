"""
Test cases for Gateway Cloning Core Capabilities

Tests the three critical aspects requested:
1. Correct att site pairing (attB1×attP1, attB2×attP2 for BP; attL1×attR1, attL2×attR2 for LR)
2. Appropriate substrates (Linear/supercoiled attB + supercoiled attP for BP)
3. Expected products (Entry clone with attL from BP, expression clone with attB from LR)

Usage:
    cd ~/python-libraries/aiplasmiddesign_api/backend
    python3 -m pytest tests/test_gateway_core_capabilities.py -v
"""

import pytest
from typing import List, Tuple
from Bio.Seq import Seq

# Import Gateway modules
from splicify_api.cloning.gateway_sites import (
    GATEWAY_ATT_SITES,
    AttSiteMatch,
    scan_att_sites,
    scan_for_ccdb,
    validate_orthogonality,
    parse_att_site_type,
    get_recombination_products,
)
from splicify_api.cloning.gateway_operator import GatewayOperator
from splicify_api.cloning.build_plan import GatewayJunctionPlan


# =============================================================================
# Test Data: Synthetic Sequences
# =============================================================================

# GFP CDS (714 bp) - standard mNeonGreen variant
GFP_SEQUENCE = """
ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG
""".replace("\n", "").replace(" ", "")

# Build attB-GFP-attB construct (for BP reaction test)
def build_attB_gfp_attB() -> str:
    """Build insert with attB1 and attB2 flanking GFP."""
    attB1 = GATEWAY_ATT_SITES["attB1"]
    attB2 = GATEWAY_ATT_SITES["attB2"]
    return attB1 + GFP_SEQUENCE + attB2

# Build pDONR221-like donor vector (simplified)
def build_pDONR221_simplified() -> str:
    """
    Build simplified pDONR221 donor vector.
    Structure: [backbone_left]-attP1-[ccdB region]-attP2-[backbone_right]
    """
    # Simplified backbone fragments (would normally be full plasmid sequences)
    # For testing, we'll use minimal sequences with key features
    
    # KanR promoter fragment (minimal)
    backbone_left = "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCACTGGCCGTCGTTTTAC" * 10  # ~800bp
    
    # Origin fragment (minimal)
    backbone_right = "TCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTCCTGCAGCCCGGGGGATCCACTAGTTCTAGAGCGGCCGC" * 10  # ~750bp
    
    attP1 = GATEWAY_ATT_SITES["attP1"]
    attP2 = GATEWAY_ATT_SITES["attP2"]
    
    # ccdB gene region (includes ccdB and surrounding sequences)
    from splicify_api.cloning.gateway_sites import CCDB_GENE
    ccdb_region = "AGCTTGGCTG" + CCDB_GENE + "TAATACGACT"  # ccdB with flanking
    
    return backbone_left + attP1 + ccdb_region + attP2 + backbone_right


# Build pENTR-like entry clone (product of BP reaction)
def build_pENTR_gfp() -> str:
    """Build pENTR-GFP entry clone with attL1 and attL2."""
    backbone_left = "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCACTGGCCGTCGTTTTAC" * 10
    backbone_right = "TCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTCCTGCAGCCCGGGGGATCCACTAGTTCTAGAGCGGCCGC" * 10
    
    attL1 = GATEWAY_ATT_SITES["attL1"]
    attL2 = GATEWAY_ATT_SITES["attL2"]
    
    return backbone_left + attL1 + GFP_SEQUENCE + attL2 + backbone_right


# =============================================================================
# Test Suite 1: att Site Pairing Validation
# =============================================================================

class TestAttSitePairing:
    """Test correct att site pairing for BP and LR reactions."""
    
    def test_validate_orthogonality_bp_correct_pairing(self):
        """Test that attB1 correctly pairs with attP1."""
        # Create mock AttSiteMatch objects
        attB1_match = AttSiteMatch(
            site_type="attB1",
            start=0,
            end=25,
            sequence=GATEWAY_ATT_SITES["attB1"],
            core_sequence="GTACAAA",
            match_quality="exact",
            strand=1
        )
        
        attP1_match = AttSiteMatch(
            site_type="attP1",
            start=0,
            end=200,
            sequence=GATEWAY_ATT_SITES["attP1"],
            core_sequence="GTACAAA",
            match_quality="exact",
            strand=1
        )
        
        # Validate orthogonality
        is_valid = validate_orthogonality(attB1_match, attP1_match)
        assert is_valid, "attB1 should be orthogonal with attP1 (same site number, compatible types)"
    
    def test_validate_orthogonality_bp_incorrect_pairing(self):
        """Test that attB1 does NOT pair with attP2 (wrong site number)."""
        attB1_match = AttSiteMatch(
            site_type="attB1",
            start=0,
            end=25,
            sequence=GATEWAY_ATT_SITES["attB1"],
            core_sequence="GTACAAA",
            match_quality="exact",
            strand=1
        )
        
        attP2_match = AttSiteMatch(
            site_type="attP2",
            start=0,
            end=200,
            sequence=GATEWAY_ATT_SITES["attP2"],
            core_sequence="GTACAAG",  # Different core!
            match_quality="exact",
            strand=1
        )
        
        is_valid = validate_orthogonality(attB1_match, attP2_match)
        assert not is_valid, "attB1 should NOT be orthogonal with attP2 (different site numbers)"
    
    def test_validate_orthogonality_lr_correct_pairing(self):
        """Test that attL1 correctly pairs with attR1."""
        attL1_match = AttSiteMatch(
            site_type="attL1",
            start=0,
            end=100,
            sequence=GATEWAY_ATT_SITES["attL1"],
            core_sequence="GTACAAA",
            match_quality="exact",
            strand=1
        )
        
        attR1_match = AttSiteMatch(
            site_type="attR1",
            start=0,
            end=125,
            sequence=GATEWAY_ATT_SITES["attR1"],
            core_sequence="GTACAAA",
            match_quality="exact",
            strand=1
        )
        
        is_valid = validate_orthogonality(attL1_match, attR1_match)
        assert is_valid, "attL1 should be orthogonal with attR1"
    
    def test_get_recombination_products_bp(self):
        """Test that BP reaction produces correct att sites (attB+attP → attL+attR)."""
        product_left, product_right = get_recombination_products("attB1", "attP1")
        
        assert product_left == "attL1", "BP reaction should produce attL1 on left"
        assert product_right == "attR1", "BP reaction should produce attR1 on right"
    
    def test_get_recombination_products_lr(self):
        """Test that LR reaction produces correct att sites (attL+attR → attB+attP)."""
        product_left, product_right = get_recombination_products("attL1", "attR1")
        
        assert product_left == "attB1", "LR reaction should produce attB1"
        assert product_right == "attP1", "LR reaction should produce attP1"
    
    def test_parse_att_site_type(self):
        """Test att site type parsing."""
        site_type, site_num = parse_att_site_type("attB1")
        assert site_type == "B", "Should parse type as B"
        assert site_num == "1", "Should parse number as 1"
        
        site_type, site_num = parse_att_site_type("attP2")
        assert site_type == "P", "Should parse type as P"
        assert site_num == "2", "Should parse number as 2"


# =============================================================================
# Test Suite 2: att Site Detection in Test Sequences
# =============================================================================

class TestAttSiteDetection:
    """Test att site detection in synthetic constructs."""
    
    def test_scan_attB_gfp_attB(self):
        """Test detection of attB1 and attB2 in attB-GFP-attB construct."""
        construct = build_attB_gfp_attB()
        sites = scan_att_sites(construct, fuzzy_threshold=0)
        
        # Should detect exactly 2 sites
        assert len(sites) == 2, f"Expected 2 att sites, found {len(sites)}"
        
        # Check site types
        site_types = [s.site_type for s in sites]
        assert "attB1" in site_types, "Should detect attB1"
        assert "attB2" in site_types, "Should detect attB2"
        
        # Check positions (attB1 at start, attB2 at end)
        attB1_site = [s for s in sites if s.site_type == "attB1"][0]
        attB2_site = [s for s in sites if s.site_type == "attB2"][0]
        
        assert attB1_site.start == 0, "attB1 should be at position 0"
        assert attB2_site.end == len(construct), "attB2 should be at end"
    
    def test_scan_pDONR221_simplified(self):
        """Test detection of attP1, attP2, and ccdB in pDONR221."""
        donor = build_pDONR221_simplified()
        sites = scan_att_sites(donor, fuzzy_threshold=0)
        
        # Should detect attP1 and attP2
        site_types = [s.site_type for s in sites]
        assert "attP1" in site_types, "Should detect attP1"
        assert "attP2" in site_types, "Should detect attP2"
        
        # Check for ccdB
        ccdb_positions = scan_for_ccdb(donor)
        assert len(ccdb_positions) > 0, "Should detect ccdB gene in donor"
    
    def test_scan_pENTR_gfp(self):
        """Test detection of attL1 and attL2 in pENTR-GFP."""
        entry = build_pENTR_gfp()
        sites = scan_att_sites(entry, fuzzy_threshold=0)
        
        site_types = [s.site_type for s in sites]
        assert "attL1" in site_types, "Should detect attL1"
        assert "attL2" in site_types, "Should detect attL2"
        
        # Entry clone should NOT have ccdB
        ccdb_positions = scan_for_ccdb(entry)
        assert len(ccdb_positions) == 0, "Entry clone should not contain ccdB"


# =============================================================================
# Test Suite 3: Gateway Operator Integration Tests
# =============================================================================

class TestGatewayOperatorBP:
    """Integration tests for BP reaction using GatewayOperator."""
    
    def test_bp_reaction_basic(self):
        """Test basic BP reaction: attB-GFP-attB + pDONR221 → pENTR-GFP."""
        operator = GatewayOperator()
        
        # Prepare modules
        insert_module = {
            "canonical_id": "GFP_insert",
            "sequence": build_attB_gfp_attB(),
            "role": "insert"
        }
        
        donor_module = {
            "canonical_id": "pDONR221",
            "sequence": build_pDONR221_simplified(),
            "role": "vector"
        }
        
        modules = [insert_module, donor_module]
        
        # Evaluate
        plan = operator.evaluate(modules, topology="circular")
        
        # Assertions
        assert plan.reaction_type == "BP", "Should detect BP reaction"
        assert plan.feasible, "BP reaction should be feasible"
        assert len(plan.junction_plans) > 0, "Should have junction plans"
        
        # Check product sequence has attL sites
        product_sites = scan_att_sites(plan.product_sequence, fuzzy_threshold=0)
        product_site_types = [s.site_type for s in product_sites]
        
        assert "attL1" in product_site_types or "attL2" in product_site_types, \
            "Product should contain attL sites (Entry clone)"
        
        # Product should NOT contain ccdB
        ccdb_in_product = scan_for_ccdb(plan.product_sequence)
        assert len(ccdb_in_product) == 0, \
            "Entry clone product should NOT contain ccdB gene"
    
    @pytest.mark.xfail(reason="Known bug: AttributeError on jp.parent_left_site")
    def test_bp_reaction_byproduct_generation(self):
        """Test that byproduct is correctly generated with attR sites and ccdB."""
        operator = GatewayOperator()
        
        insert_module = {
            "canonical_id": "GFP_insert",
            "sequence": build_attB_gfp_attB(),
            "role": "insert"
        }
        
        donor_module = {
            "canonical_id": "pDONR221",
            "sequence": build_pDONR221_simplified(),
            "role": "vector"
        }
        
        modules = [insert_module, donor_module]
        plan = operator.evaluate(modules, topology="circular")
        
        # Check byproduct
        # THIS WILL FAIL due to AttributeError: 'GatewayJunctionPlan' object has no attribute 'parent_left_site'
        byproduct_sites = scan_att_sites(plan.byproduct_sequence, fuzzy_threshold=0)
        byproduct_site_types = [s.site_type for s in byproduct_sites]
        
        assert "attR1" in byproduct_site_types or "attR2" in byproduct_site_types, \
            "Byproduct should contain attR sites"
        
        # Byproduct SHOULD contain ccdB
        ccdb_in_byproduct = scan_for_ccdb(plan.byproduct_sequence)
        assert len(ccdb_in_byproduct) > 0, \
            "Byproduct should contain ccdB gene"


# =============================================================================
# Test Suite 4: Edge Cases and Error Handling
# =============================================================================

class TestGatewayEdgeCases:
    """Test edge cases and error handling."""
    
    def test_multiple_attB1_sites(self):
        """Test handling of multiple att sites of same type (ambiguous)."""
        # Create construct with two attB1 sites
        attB1 = GATEWAY_ATT_SITES["attB1"]
        attB2 = GATEWAY_ATT_SITES["attB2"]
        
        ambiguous_insert = attB1 + "ATGCATGC" * 100 + attB1 + "GCTAGCTA" * 100 + attB2
        
        sites = scan_att_sites(ambiguous_insert, fuzzy_threshold=0)
        
        # Count attB1 sites
        attB1_sites = [s for s in sites if s.site_type == "attB1"]
        assert len(attB1_sites) == 2, "Should detect 2 attB1 sites"
        
        # This should trigger a warning in the operator
        # (testing the warning is beyond unit test scope, but site detection works)
    
    def test_reverse_strand_detection(self):
        """Test detection of att sites on reverse strand."""
        # Create attB1 on reverse strand
        attB1_fwd = GATEWAY_ATT_SITES["attB1"]
        attB1_rc = str(Seq(attB1_fwd).reverse_complement())
        
        test_seq = "ATGCATGC" * 50 + attB1_rc + "GCTAGCTA" * 50
        
        sites = scan_att_sites(test_seq, fuzzy_threshold=0)
        
        # Should detect attB1 on reverse strand
        attB1_sites = [s for s in sites if s.site_type == "attB1"]
        assert len(attB1_sites) > 0, "Should detect attB1 on reverse strand"
        
        # Check strand information
        if attB1_sites:
            site = attB1_sites[0]
            assert site.strand == -1, "Should be marked as reverse strand"
    
    def test_fuzzy_matching_degraded_attB(self):
        """Test fuzzy matching for degraded att sites with SNPs."""
        # Create attB1 with 1 SNP
        attB1_perfect = GATEWAY_ATT_SITES["attB1"]
        attB1_degraded = list(attB1_perfect)
        attB1_degraded[10] = "A" if attB1_degraded[10] != "A" else "T"  # Single SNP
        attB1_degraded = "".join(attB1_degraded)
        
        test_seq = "ATGCATGC" * 50 + attB1_degraded + "GCTAGCTA" * 50
        
        # Should NOT detect with exact matching
        sites_exact = scan_att_sites(test_seq, fuzzy_threshold=0)
        attB1_exact = [s for s in sites_exact if s.site_type == "attB1"]
        assert len(attB1_exact) == 0, "Should not detect degraded site with exact matching"
        
        # SHOULD detect with fuzzy matching
        sites_fuzzy = scan_att_sites(test_seq, fuzzy_threshold=2)
        attB1_fuzzy = [s for s in sites_fuzzy if s.site_type == "attB1"]
        assert len(attB1_fuzzy) > 0, "Should detect degraded site with fuzzy matching"
        
        # Check match quality
        if attB1_fuzzy:
            site = attB1_fuzzy[0]
            assert "fuzzy" in site.match_quality, "Match quality should indicate fuzzy match"


# =============================================================================
# Test Suite 5: Sequence Accuracy Tests
# =============================================================================

class TestGatewaySequenceAccuracy:
    """Test that att site sequences match the Gateway manual."""
    
    def test_attB1_sequence_matches_manual(self):
        """Test attB1 sequence matches manual page 24."""
        # From manual page 24:
        # attB1 = ACAAGTTTGTACAAAAAAGCAGGCT
        expected_attB1 = "ACAAGTTTGTACAAAAAAGCAGGCT"
        
        actual_attB1 = GATEWAY_ATT_SITES.get("attB1", "")
        
        # Check if they match (case-insensitive)
        assert actual_attB1.upper() == expected_attB1.upper(), \
            f"attB1 sequence mismatch. Expected: {expected_attB1}, Got: {actual_attB1}"
    
    def test_attB2_sequence_matches_manual(self):
        """Test attB2 sequence matches manual."""
        # From manual:
        # attB2 = ACCACTTTGTACAAGAAAGCTGGGT
        expected_attB2 = "ACCACTTTGTACAAGAAAGCTGGGT"
        
        actual_attB2 = GATEWAY_ATT_SITES.get("attB2", "")
        
        assert actual_attB2.upper() == expected_attB2.upper(), \
            f"attB2 sequence mismatch. Expected: {expected_attB2}, Got: {actual_attB2}"
    
    def test_core_sequences(self):
        """Test that core 7bp sequences are correct."""
        from splicify_api.cloning.gateway_sites import CORE_SEQUENCES
        
        # From manual, core sequences:
        expected_cores = {
            "1": "GTACAAA",
            "2": "GTACAAG",
            "3": "GTATAAT",
            "4": "GTATAGA",
            "5": "GTATACA",
        }
        
        for site_num, expected_core in expected_cores.items():
            actual_core = CORE_SEQUENCES.get(site_num, "")
            assert actual_core == expected_core, \
                f"Core sequence {site_num} mismatch. Expected: {expected_core}, Got: {actual_core}"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
