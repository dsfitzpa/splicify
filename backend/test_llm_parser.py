#!/usr/bin/env python3
"""
Test LLM Module Parser on lentiCRISPR v2 plasmid
"""

import asyncio
import json
from pathlib import Path


async def test_llm_parser():
    """Test the LLM parser with a sample lentiviral plasmid."""
    print("="*80)
    print("Testing LLM Module Parser")
    print("="*80)

    # Import modules
    from splicify_api.llm_module_parser import LLMModuleParser

    # Create sample features (simplified lentiCRISPR v2)
    sample_features = [
        {"name": "5_LTR", "type": "LTR", "start": 0, "end": 634, "strand": 1, "description": "5' long terminal repeat"},
        {"name": "HIV-1_Psi", "type": "misc_feature", "start": 649, "end": 782, "strand": 1, "description": "packaging signal"},
        {"name": "RRE", "type": "misc_feature", "start": 1316, "end": 1551, "strand": 1, "description": "Rev response element"},
        {"name": "cPPT", "type": "misc_feature", "start": 1633, "end": 1751, "strand": 1, "description": "central polypurine tract"},
        {"name": "U6", "type": "promoter", "start": 1900, "end": 2157, "strand": 1, "description": "U6 promoter"},
        {"name": "gRNA_scaffold", "type": "misc_RNA", "start": 2200, "end": 2275, "strand": 1, "description": "guide RNA scaffold"},
        {"name": "CMV", "type": "promoter", "start": 2800, "end": 3389, "strand": 1, "description": "CMV promoter"},
        {"name": "NLS", "type": "misc_feature", "start": 3450, "end": 3470, "strand": 1, "description": "nuclear localization signal"},
        {"name": "Cas9", "type": "CDS", "start": 3471, "end": 7579, "strand": 1, "description": "Cas9 nuclease"},
        {"name": "NLS", "type": "misc_feature", "start": 7580, "end": 7600, "strand": 1, "description": "nuclear localization signal"},
        {"name": "P2A", "type": "misc_feature", "start": 7610, "end": 7675, "strand": 1, "description": "2A self-cleaving peptide"},
        {"name": "Puromycin", "type": "CDS", "start": 7676, "end": 8270, "strand": 1, "description": "puromycin resistance"},
        {"name": "WPRE", "type": "misc_feature", "start": 8400, "end": 8989, "strand": 1, "description": "woodchuck post-transcriptional regulatory element"},
        {"name": "bGH_polyA", "type": "polyA_signal", "start": 9020, "end": 9245, "strand": 1, "description": "bovine growth hormone polyadenylation signal"},
        {"name": "3_LTR", "type": "LTR", "start": 9400, "end": 10034, "strand": 1, "description": "3' long terminal repeat"},
        {"name": "Amp_promoter", "type": "promoter", "start": 10200, "end": 10304, "strand": -1, "description": "ampicillin resistance promoter"},
        {"name": "AmpR", "type": "CDS", "start": 10305, "end": 11165, "strand": -1, "description": "ampicillin resistance"},
        {"name": "ori", "type": "rep_origin", "start": 11336, "end": 11924, "strand": -1, "description": "ColE1 origin of replication"},
    ]

    # Create mock sequence (simplified)
    sequence_length = 12000
    sequence = "A" * sequence_length

    print(f"\nTest plasmid: {sequence_length} bp circular")
    print(f"Features: {len(sample_features)}")
    print("\nFeature breakdown:")
    print(f"  - Viral elements (LTR, Psi, RRE, cPPT, WPRE): 6")
    print(f"  - Promoters (U6, CMV): 2")
    print(f"  - Coding sequences (Cas9, Puromycin): 2")
    print(f"  - Regulatory (P2A, NLS): 3")
    print(f"  - Backbone (AmpR, ori): 3")

    # Initialize parser
    print("\n" + "="*80)
    print("Initializing LLM Module Parser...")
    print("="*80)

    try:
        parser = LLMModuleParser(model="gpt-4o")
        print("✓ Parser initialized with GPT-4o")
    except ValueError as e:
        print(f"✗ Failed to initialize parser: {e}")
        print("\nPlease set OPENAI_API_KEY environment variable:")
        print("  export OPENAI_API_KEY='your-key-here'")
        return

    # Parse modules
    print("\n" + "="*80)
    print("Calling LLM to identify modules...")
    print("="*80)

    result = await parser.parse_modules(
        sequence=sequence,
        plannotate_features=sample_features,
        cds_submodules=None,
        circular=True
    )

    # Display results
    modules = result.get("modules", [])
    print(f"\n✓ LLM identified {len(modules)} modules\n")

    print("Modules detected:")
    print("-" * 80)
    for i, mod in enumerate(modules, 1):
        mod_type = mod['module_type']
        start = mod['start']
        end = mod['end']
        parent = mod.get('parent_module', 'None')
        nested = len(mod.get('nested_modules', []))

        indent = "  " if parent != 'None' else ""
        print(f"{indent}{i}. {mod_type}")
        print(f"{indent}   Position: {start}..{end} ({end-start} bp)")
        print(f"{indent}   Parent: {parent}")
        print(f"{indent}   Nested modules: {nested}")

        # Show key metadata
        metadata = mod.get('metadata', {})
        if metadata:
            if 'promoter_id' in metadata:
                print(f"{indent}   Promoter: {metadata['promoter_id']}")
            if 'polya_id' in metadata:
                print(f"{indent}   PolyA: {metadata['polya_id']}")
            if 'protein_name' in metadata:
                print(f"{indent}   Protein: {metadata['protein_name']}")
        print()

    # Verify hierarchy
    print("\nModule Hierarchy:")
    print("-" * 80)

    # Find top-level modules (no parent)
    top_level = [m for m in modules if not m.get('parent_module')]

    def print_hierarchy(mod, level=0):
        indent = "  " * level
        symbol = "└─" if level > 0 else "●"
        print(f"{indent}{symbol} {mod['module_type']} ({mod['start']}..{mod['end']})")

        # Print nested modules
        for child_id in mod.get('nested_modules', []):
            child = next((m for m in modules if m['module_id'] == child_id), None)
            if child:
                print_hierarchy(child, level + 1)

    for top_mod in top_level:
        print_hierarchy(top_mod)

    # Summary
    print("\n" + "="*80)
    print("Summary:")
    print("="*80)
    summary = result.get('summary', {})
    print(f"Total modules: {summary.get('module_count', 0)}")
    print(f"Module types: {', '.join(summary.get('module_types', []))}")
    print(f"Source: {summary.get('source_pipeline', 'unknown')}")
    print(f"Model: {summary.get('model', 'unknown')}")

    # Check for key biological understanding
    print("\n" + "="*80)
    print("Biological Understanding Check:")
    print("="*80)

    checks = []

    # Check 1: Lentiviral payload detected
    lenti_modules = [m for m in modules if 'lentiviral' in m['module_type'].lower()]
    if lenti_modules:
        checks.append("✓ Lentiviral payload identified")
    else:
        checks.append("✗ Lentiviral payload NOT identified")

    # Check 2: Expression cassettes
    expr_modules = [m for m in modules if 'expression' in m['module_type'].lower()]
    if len(expr_modules) >= 2:
        checks.append(f"✓ Expression cassettes identified ({len(expr_modules)})")
    else:
        checks.append("✗ Expression cassettes NOT identified")

    # Check 3: Pol II vs Pol III distinction
    pol2 = [m for m in modules if 'pol2' in m['module_type'].lower()]
    pol3 = [m for m in modules if 'pol3' in m['module_type'].lower()]
    if pol2 and pol3:
        checks.append("✓ Pol II and Pol III cassettes distinguished")
    else:
        checks.append("✗ Pol II/Pol III NOT distinguished")

    # Check 4: Bacterial backbone
    backbone = [m for m in modules if 'backbone' in m['module_type'].lower()]
    if backbone:
        checks.append("✓ Bacterial backbone identified")
    else:
        checks.append("✗ Bacterial backbone NOT identified")

    for check in checks:
        print(check)

    print("\n" + "="*80)
    print("Test complete!")
    print("="*80)

    return result


if __name__ == "__main__":
    asyncio.run(test_llm_parser())
