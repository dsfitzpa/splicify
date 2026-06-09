#!/usr/bin/env python3
"""
Simple Gateway Cloning Tests (no pytest required)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from splicify_api.cloning.gateway_sites import (
    GATEWAY_ATT_SITES,
    scan_att_sites,
    validate_orthogonality,
    parse_att_site_type,
    get_recombination_products,
    AttSiteMatch,
)


def test_validate_orthogonality_bp():
    """Test that attB1 pairs with attP1."""
    print("Test: attB1 should pair with attP1...", end=" ")
    
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
    
    result = validate_orthogonality(attB1_match, attP1_match)
    assert result == True, "attB1 should be orthogonal with attP1"
    print("✓ PASS")
    return True


def test_validate_orthogonality_bp_mismatch():
    """Test that attB1 does NOT pair with attP2."""
    print("Test: attB1 should NOT pair with attP2...", end=" ")
    
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
        core_sequence="GTACAAG",  # Different core
        match_quality="exact",
        strand=1
    )
    
    result = validate_orthogonality(attB1_match, attP2_match)
    assert result == False, "attB1 should NOT be orthogonal with attP2"
    print("✓ PASS")
    return True


def test_recombination_products_bp():
    """Test BP reaction produces attL and attR."""
    print("Test: BP reaction (attB1+attP1) should produce attL1 and attR1...", end=" ")
    
    product_left, product_right = get_recombination_products("attB1", "attP1")
    
    assert product_left == "attL1", f"Expected attL1, got {product_left}"
    assert product_right == "attR1", f"Expected attR1, got {product_right}"
    print("✓ PASS")
    return True


def test_recombination_products_lr():
    """Test LR reaction produces attB and attP."""
    print("Test: LR reaction (attL1+attR1) should produce attB1 and attP1...", end=" ")
    
    product_left, product_right = get_recombination_products("attL1", "attR1")
    
    assert product_left == "attB1", f"Expected attB1, got {product_left}"
    assert product_right == "attP1", f"Expected attP1, got {product_right}"
    print("✓ PASS")
    return True


def test_parse_att_site_type():
    """Test att site type parsing."""
    print("Test: Parse att site types...", end=" ")
    
    site_type, site_num = parse_att_site_type("attB1")
    assert site_type == "B" and site_num == "1", "Should parse attB1 correctly"
    
    site_type, site_num = parse_att_site_type("attP2")
    assert site_type == "P" and site_num == "2", "Should parse attP2 correctly"
    
    site_type, site_num = parse_att_site_type("attL1")
    assert site_type == "L" and site_num == "1", "Should parse attL1 correctly"
    
    site_type, site_num = parse_att_site_type("attR2")
    assert site_type == "R" and site_num == "2", "Should parse attR2 correctly"
    
    print("✓ PASS")
    return True


def test_att_site_detection():
    """Test att site detection in sequence."""
    print("Test: Detect att sites in sequence...", end=" ")
    
    # Build sequence with attB1 and attB2
    attB1_seq = GATEWAY_ATT_SITES["attB1"]
    attB2_seq = GATEWAY_ATT_SITES["attB2"]
    test_seq = attB1_seq + ("ATGC" * 100) + attB2_seq
    
    sites = scan_att_sites(test_seq, fuzzy_threshold=0)
    site_types = [s.site_type for s in sites]
    
    assert "attB1" in site_types, "Should detect attB1"
    assert "attB2" in site_types, "Should detect attB2"
    
    print("✓ PASS")
    return True


def main():
    print("\n" + "=" * 60)
    print("Gateway Cloning Core Capabilities Tests")
    print("=" * 60 + "\n")
    
    tests = [
        test_validate_orthogonality_bp,
        test_validate_orthogonality_bp_mismatch,
        test_recombination_products_bp,
        test_recombination_products_lr,
        test_parse_att_site_type,
        test_att_site_detection,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"✗ FAIL: {e}")
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60 + "\n")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
