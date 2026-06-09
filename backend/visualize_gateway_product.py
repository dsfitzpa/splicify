#!/usr/bin/env python3
"""
Gateway Product Visualization Script

Creates a detailed visual representation of the Gateway cloning product
showing all fragments, att sites, and their origins.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from splicify_api.cloning.gateway_sites import (
    GATEWAY_ATT_SITES,
    scan_att_sites,
    scan_for_ccdb,
)
from splicify_api.cloning.gateway_operator import GatewayOperator


def print_plasmid_map(sequence, name="Product"):
    """Print a visual plasmid map with att sites and features."""
    print(f"\n{'='*80}")
    print(f"Plasmid Map: {name}")
    print(f"Total Size: {len(sequence):,} bp")
    print(f"{'='*80}\n")
    
    # Scan for att sites
    att_sites = scan_att_sites(sequence, fuzzy_threshold=0)
    ccdb_sites = scan_for_ccdb(sequence)
    
    # Create position markers
    positions = []
    
    # Add att sites
    for site in att_sites:
        positions.append({
            'start': site.start,
            'end': site.end,
            'name': site.site_type,
            'type': 'att_site',
            'core': site.core_sequence
        })
    
    # Add ccdB
    for start, end in ccdb_sites:
        positions.append({
            'start': start,
            'end': end,
            'name': 'ccdB',
            'type': 'gene'
        })
    
    # Sort by position
    positions.sort(key=lambda x: x['start'])
    
    # Print map
    print("Position Map:")
    print("=" * 80)
    
    if not positions:
        print("  No features detected")
    else:
        prev_end = 0
        for i, feat in enumerate(positions):
            # Print gap
            if feat['start'] > prev_end:
                gap_size = feat['start'] - prev_end
                print(f"  [{prev_end:>6,} - {feat['start']:>6,}]  ({gap_size:>5} bp)  Backbone sequence")
            
            # Print feature
            feat_size = feat['end'] - feat['start']
            if feat['type'] == 'att_site':
                print(f"  [{feat['start']:>6,} - {feat['end']:>6,}]  ({feat_size:>5} bp)  ★ {feat['name']} (core: {feat['core']})")
            elif feat['type'] == 'gene':
                print(f"  [{feat['start']:>6,} - {feat['end']:>6,}]  ({feat_size:>5} bp)  ⚠ {feat['name']} gene")
            
            prev_end = feat['end']
        
        # Print final gap
        if prev_end < len(sequence):
            gap_size = len(sequence) - prev_end
            print(f"  [{prev_end:>6,} - {len(sequence):>6,}]  ({gap_size:>5} bp)  Backbone sequence")
    
    print("=" * 80)
    
    # Summary
    print(f"\nSummary:")
    att_counts = {}
    for site in att_sites:
        site_type = site.site_type[:4]  # attB, attP, attL, attR
        att_counts[site_type] = att_counts.get(site_type, 0) + 1
    
    print(f"  att sites found: {len(att_sites)}")
    for site_type, count in sorted(att_counts.items()):
        print(f"    - {site_type}: {count}")
    
    print(f"  ccdB gene: {'Found' if ccdb_sites else 'Not detected'}")


def print_fragment_diagram(annotations):
    """Print a diagram showing fragment origins."""
    print(f"\n{'='*80}")
    print("Fragment Origins Diagram")
    print(f"{'='*80}\n")
    
    # Group by source type
    by_type = {}
    for annot in annotations:
        source_type = annot.get('type', 'unknown')
        if source_type not in by_type:
            by_type[source_type] = []
        by_type[source_type].append(annot)
    
    # Print each type
    type_symbols = {
        'donor_backbone': '▓',
        'insert': '█',
        'att_site': '▒',
        'recombination_site': '░'
    }
    
    for source_type, annots in sorted(by_type.items()):
        symbol = type_symbols.get(source_type, '·')
        print(f"\n{source_type.upper()} ({symbol}):" )
        for annot in annots:
            name = annot.get('name', 'Unknown')
            start = annot.get('start', 0)
            end = annot.get('end', 0)
            size = end - start
            source = annot.get('source', 'Unknown')
            print(f"  {symbol * 3} [{start:>6,} - {end:>6,}] ({size:>5} bp)  {name}  (from: {source})")


def main():
    print("\n" + "#" * 80)
    print("#" + " " * 78 + "#")
    print("#" + "  Gateway Cloning Product Visualization".center(78) + "#")
    print("#" + " " * 78 + "#")
    print("#" * 80)
    
    # Build test sequences
    GFP_SEQUENCE = (
        "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG"
    )
    
    from splicify_api.cloning.gateway_sites import CCDB_GENE
    
    # Build pDONR221
    backbone_left = "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCACTGGCCGTCGTTTTAC" * 10
    backbone_right = "TCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTCCTGCAGCCCGGGGGATCCACTAGTTCTAGAGCGGCCGC" * 10
    attP1 = GATEWAY_ATT_SITES["attP1"]
    attP2 = GATEWAY_ATT_SITES["attP2"]
    ccdb_region = "AGCTTGGCTG" + CCDB_GENE + "TAATACGACT"
    pdonr221_seq = backbone_left + attP1 + ccdb_region + attP2 + backbone_right
    
    # Build PCR product with attB sites
    attB1 = GATEWAY_ATT_SITES["attB1"]
    attB2 = GATEWAY_ATT_SITES["attB2"]
    pcr_product = attB1 + GFP_SEQUENCE + attB2
    
    # Run Gateway operator
    operator = GatewayOperator()
    
    insert_module = {
        "canonical_id": "GFP_insert",
        "sequence": pcr_product,
        "role": "insert"
    }
    
    donor_module = {
        "canonical_id": "pDONR221",
        "sequence": pdonr221_seq,
        "role": "vector"
    }
    
    modules = [insert_module, donor_module]
    
    print("\n" + "=" * 80)
    print("INPUT PLASMIDS")
    print("=" * 80)
    
    # Visualize inputs
    print_plasmid_map(pcr_product, "GFP Insert (attB1-GFP-attB2)")
    print_plasmid_map(pdonr221_seq, "pDONR221 Donor Vector")
    
    # Run Gateway cloning
    print("\n\n" + "=" * 80)
    print("RUNNING GATEWAY BP REACTION...")
    print("=" * 80)
    
    plan = operator.evaluate(modules, topology="circular")
    
    print(f"\nReaction Type: {plan.reaction_type}")
    print(f"Feasible: {plan.feasible}")
    print(f"Warnings: {len(plan.warnings)}")
    
    # Visualize products
    print("\n\n" + "=" * 80)
    print("OUTPUT PLASMIDS")
    print("=" * 80)
    
    print_plasmid_map(plan.product_sequence, "Entry Clone (pENTR-GFP)")
    
    if plan.byproduct_sequence:
        print_plasmid_map(plan.byproduct_sequence, "Byproduct (Linear Fragment)")
    
    # Show fragment origins
    if plan.product_annotations:
        print_fragment_diagram(plan.product_annotations)
    
    # Final verification
    print("\n\n" + "=" * 80)
    print("VERIFICATION")
    print("=" * 80 + "\n")
    
    product_sites = scan_att_sites(plan.product_sequence)
    product_site_types = [s.site_type for s in product_sites]
    product_ccdb = scan_for_ccdb(plan.product_sequence)
    
    byproduct_sites = scan_att_sites(plan.byproduct_sequence) if plan.byproduct_sequence else []
    byproduct_site_types = [s.site_type for s in byproduct_sites]
    byproduct_ccdb = scan_for_ccdb(plan.byproduct_sequence) if plan.byproduct_sequence else []
    
    checks = [
        ("✓" if plan.reaction_type == "BP" else "✗", "Reaction type is BP"),
        ("✓" if any('attL' in s for s in product_site_types) else "✗", "Product has attL sites (Entry clone)"),
        ("✓" if len(product_ccdb) == 0 else "✗", "Product does NOT contain ccdB"),
        ("✓" if any('attR' in s for s in byproduct_site_types) else "✗", "Byproduct has attR sites"),
        ("✓" if len(byproduct_ccdb) > 0 else "✗", "Byproduct contains ccdB"),
        ("✓" if GFP_SEQUENCE in plan.product_sequence else "✗", "Product contains complete GFP insert"),
    ]
    
    print("Verification Checklist:")
    all_pass = True
    for symbol, check in checks:
        print(f"  {symbol} {check}")
        if symbol == "✗":
            all_pass = False
    
    print("\n" + "=" * 80)
    if all_pass:
        print("\n🎉 SUCCESS! Gateway cloning product is correct.\n")
    else:
        print("\n⚠️  ISSUES DETECTED. Review the product.\n")
    print("=" * 80 + "\n")
    
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
