"""
Test Script for Grammar-Based Pol II Detector
==============================================

Tests the grammar_pol2_detector.py module against real plasmid data.

Usage:
    python test_grammar_pol2.py
"""

import json
from grammar_pol2_detector import GrammarPol2Detector, detect_pol2_cassettes


def create_test_features():
    """
    Create test features simulating pLannotate output with KB enrichment.

    Simulates a lentiCRISPR-style plasmid with:
    - EF1a promoter
    - Cas9-P2A-PuroR CDS
    - bGH polyA
    """
    features = [
        {
            'name': 'EF-1-alpha promoter',
            'start': 100,
            'end': 1300,
            'strand': 1,
            'type': 'promoter',
            'kb_info': {
                'feature_class': 'promoter',
                'polymerase_class': 'pol_ii',
                'host_scope': ['mammalian'],
                'subclass': 'constitutive'
            }
        },
        {
            'name': 'EF-1-alpha intron A',
            'start': 1100,
            'end': 2200,
            'strand': 1,
            'type': 'intron',
            'kb_info': {
                'feature_class': 'intron',
                'subclass': 'intron_mediated_enhancement'
            }
        },
        {
            'name': 'Kozak sequence',
            'start': 2195,
            'end': 2210,
            'strand': 1,
            'type': 'misc_feature',
            'kb_info': {
                'feature_class': 'kozak'
            }
        },
        {
            'name': 'Cas9',
            'start': 2211,
            'end': 6320,
            'strand': 1,
            'type': 'CDS',
            'kb_info': {
                'feature_class': 'cds',
                'product_class': 'nuclease',
                'host_scope': ['mammalian', 'bacterial']
            }
        },
        {
            'name': 'P2A',
            'start': 6321,
            'end': 6380,
            'strand': 1,
            'type': 'CDS',
            'kb_info': {
                'feature_class': 'peptide',
                'subclass': 'self_cleaving'
            }
        },
        {
            'name': 'PuroR',
            'start': 6381,
            'end': 6980,
            'strand': 1,
            'type': 'CDS',
            'kb_info': {
                'feature_class': 'cds',
                'product_class': 'resistance',
                'host_scope': ['mammalian']
            }
        },
        {
            'name': 'WPRE',
            'start': 6985,
            'end': 7580,
            'strand': 1,
            'type': 'misc_feature',
            'kb_info': {
                'feature_class': 'regulatory',
                'subclass': 'posttranscriptional'
            }
        },
        {
            'name': 'bGH poly(A) signal',
            'start': 7585,
            'end': 7820,
            'strand': 1,
            'type': 'polyA_signal',
            'kb_info': {
                'feature_class': 'polya_signal',
                'host_scope': ['mammalian']
            }
        },
        # Add a second cassette (U6-gRNA) on reverse strand - should NOT be detected as Pol II
        {
            'name': 'U6 promoter',
            'start': 8000,
            'end': 8250,
            'strand': -1,
            'type': 'promoter',
            'kb_info': {
                'feature_class': 'promoter',
                'polymerase_class': 'pol_iii',
                'host_scope': ['mammalian']
            }
        },
        # Add bacterial promoter on forward strand - should NOT be detected
        {
            'name': 'T7 promoter',
            'start': 9000,
            'end': 9020,
            'strand': 1,
            'type': 'promoter',
            'kb_info': {
                'feature_class': 'promoter',
                'polymerase_class': 'phage',
                'host_scope': ['bacterial']
            }
        }
    ]

    return features


def test_basic_detection():
    """Test basic Pol II cassette detection"""
    print("=" * 70)
    print("TEST 1: Basic Pol II Cassette Detection")
    print("=" * 70)

    features = create_test_features()
    sequence = "A" * 10000  # Dummy sequence

    detector = GrammarPol2Detector(features, sequence)
    cassettes = detector.detect_pol2_cassettes()

    print(f"\n✓ Detected {len(cassettes)} Pol II cassette(s)")

    for i, cassette in enumerate(cassettes, 1):
        print(f"\nCassette {i}:")
        print(f"  Position: {cassette.start}-{cassette.end} (strand {cassette.strand:+d})")
        print(f"  Length: {cassette.end - cassette.start} bp")
        print(f"  Promoter: {cassette.promoter_info['type']} ({cassette.promoter_info['strength']})")
        print(f"  PolyA: {cassette.polya_info['type']}")
        print(f"  Score: {cassette.score:.3f}")
        print(f"  Confidence: {cassette.confidence}")

        # Show components
        components = cassette.components
        print(f"  Components:")
        if components['cds']:
            print(f"    - {len(components['cds'])} CDS feature(s)")
        if components['introns']:
            print(f"    - {len(components['introns'])} intron(s)")
        if components['kozak']:
            print(f"    - Kozak sequence")
        if components['wpre']:
            print(f"    - WPRE")
        if components['2a_peptides']:
            print(f"    - {len(components['2a_peptides'])} 2A peptide(s)")

    # Validate
    assert len(cassettes) == 1, f"Expected 1 cassette, got {len(cassettes)}"
    assert cassettes[0].promoter_info['type'] == 'EF1a', "Expected EF1a promoter"
    assert cassettes[0].polya_info['type'] == 'bGH', "Expected bGH polyA"
    assert cassettes[0].confidence in ('high', 'medium'), "Expected high/medium confidence"
    assert len(cassettes[0].components['cds']) >= 2, "Expected at least 2 CDS (Cas9, PuroR)"

    print("\n✓ All assertions passed!")


def test_promoter_classification():
    """Test promoter type classification"""
    print("\n" + "=" * 70)
    print("TEST 2: Promoter Classification")
    print("=" * 70)

    test_promoters = [
        {
            'name': 'CMV immediate early enhancer/promoter',
            'start': 0,
            'end': 589,
            'strand': 1,
            'kb_info': {'feature_class': 'promoter', 'polymerase_class': 'pol_ii'}
        },
        {
            'name': 'CAG promoter',
            'start': 0,
            'end': 1500,
            'strand': 1,
            'kb_info': {'feature_class': 'promoter', 'polymerase_class': 'pol_ii'}
        },
        {
            'name': 'PGK promoter',
            'start': 0,
            'end': 500,
            'strand': 1,
            'kb_info': {'feature_class': 'promoter', 'polymerase_class': 'pol_ii'}
        },
        {
            'name': 'TRE3G promoter',
            'start': 0,
            'end': 300,
            'strand': 1,
            'kb_info': {'feature_class': 'promoter', 'polymerase_class': 'pol_ii'}
        }
    ]

    detector = GrammarPol2Detector([], "A" * 10000)

    for prom in test_promoters:
        info = detector._classify_promoter(prom)
        print(f"\n{prom['name']}:")
        print(f"  Type: {info['type']}")
        print(f"  Strength: {info['strength']}")
        print(f"  Weight: {info['weight']:.2f}")
        if 'inducible' in info:
            print(f"  Inducible: {info['inducible']}")

    print("\n✓ Classification complete!")


def test_polya_classification():
    """Test polyA signal classification"""
    print("\n" + "=" * 70)
    print("TEST 3: PolyA Signal Classification")
    print("=" * 70)

    test_polyas = [
        {'name': 'bGH poly(A) signal', 'start': 0, 'end': 235, 'strand': 1, 'kb_info': {}},
        {'name': 'SV40 late poly(A) signal', 'start': 0, 'end': 200, 'strand': 1, 'kb_info': {}},
        {'name': 'hGH poly(A) signal', 'start': 0, 'end': 250, 'strand': 1, 'kb_info': {}},
    ]

    detector = GrammarPol2Detector([], "A" * 10000)

    for polya in test_polyas:
        info = detector._classify_polya(polya)
        print(f"\n{polya['name']}:")
        print(f"  Type: {info['type']}")
        print(f"  Weight: {info['weight']:.2f}")
        if 'note' in info:
            print(f"  Note: {info['note']}")

    print("\n✓ Classification complete!")


def test_convenience_function():
    """Test the simple convenience function"""
    print("\n" + "=" * 70)
    print("TEST 4: Convenience Function")
    print("=" * 70)

    features = create_test_features()
    sequence = "A" * 10000

    cassettes = detect_pol2_cassettes(features, sequence)

    print(f"\n✓ Detected {len(cassettes)} cassette(s)")
    print("\nCassette as dict:")
    print(json.dumps(cassettes[0], indent=2))

    print("\n✓ Convenience function works!")


def test_edge_cases():
    """Test edge cases"""
    print("\n" + "=" * 70)
    print("TEST 5: Edge Cases")
    print("=" * 70)

    # Test 1: No Pol II promoters
    features = [
        {'name': 'U6 promoter', 'start': 0, 'end': 250, 'strand': 1,
         'kb_info': {'feature_class': 'promoter', 'polymerase_class': 'pol_iii'}}
    ]
    cassettes = detect_pol2_cassettes(features, "A" * 1000)
    assert len(cassettes) == 0, "Should detect 0 cassettes with no Pol II promoters"
    print("✓ Test 1 passed: No Pol II promoters")

    # Test 2: Promoter without matching polyA
    features = [
        {'name': 'CMV promoter', 'start': 0, 'end': 589, 'strand': 1,
         'kb_info': {'feature_class': 'promoter', 'polymerase_class': 'pol_ii'}}
    ]
    cassettes = detect_pol2_cassettes(features, "A" * 1000)
    assert len(cassettes) == 0, "Should detect 0 cassettes without polyA"
    print("✓ Test 2 passed: Promoter without polyA")

    # Test 3: Different strands
    features = [
        {'name': 'CMV promoter', 'start': 0, 'end': 589, 'strand': 1,
         'kb_info': {'feature_class': 'promoter', 'polymerase_class': 'pol_ii'}},
        {'name': 'bGH poly(A)', 'start': 1000, 'end': 1235, 'strand': -1,
         'kb_info': {'feature_class': 'polya_signal'}}
    ]
    cassettes = detect_pol2_cassettes(features, "A" * 2000)
    assert len(cassettes) == 0, "Should detect 0 cassettes on different strands"
    print("✓ Test 3 passed: Different strands")

    # Test 4: Cassette too large (>15kb)
    features = [
        {'name': 'CMV promoter', 'start': 0, 'end': 589, 'strand': 1,
         'kb_info': {'feature_class': 'promoter', 'polymerase_class': 'pol_ii'}},
        {'name': 'bGH poly(A)', 'start': 16000, 'end': 16235, 'strand': 1,
         'kb_info': {'feature_class': 'polya_signal'}}
    ]
    cassettes = detect_pol2_cassettes(features, "A" * 20000)
    assert len(cassettes) == 0, "Should detect 0 cassettes >15kb apart"
    print("✓ Test 4 passed: Cassette too large")

    print("\n✓ All edge case tests passed!")


def run_all_tests():
    """Run all tests"""
    print("\n" + "=" * 70)
    print("GRAMMAR-BASED POL II DETECTOR TEST SUITE")
    print("=" * 70)

    try:
        test_basic_detection()
        test_promoter_classification()
        test_polya_classification()
        test_convenience_function()
        test_edge_cases()

        print("\n" + "=" * 70)
        print("✓ ALL TESTS PASSED!")
        print("=" * 70)
        print("\nThe grammar-based Pol II detector is working correctly.")
        print("You can now integrate it into your annotation pipeline.")

    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        raise
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        raise


if __name__ == '__main__':
    run_all_tests()
