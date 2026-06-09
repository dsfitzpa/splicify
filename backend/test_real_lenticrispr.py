#!/usr/bin/env python3
"""
Test LLM Module Parser on real lentiCRISPR v2 plasmid
"""

import asyncio
import json
from pathlib import Path
from Bio import SeqIO


async def test_real_plasmid():
    """Test LLM parser on actual lentiCRISPR v2 .gb file."""
    print("="*80)
    print("Testing LLM Module Parser on Real lentiCRISPR v2")
    print("="*80)

    # Read GenBank file
    gb_path = "/tmp/lenticrispr_v2.gb"

    try:
        record = SeqIO.read(gb_path, "genbank")
    except FileNotFoundError:
        print(f"✗ File not found: {gb_path}")
        print("Please ensure the file exists locally")
        return

    sequence = str(record.seq)
    seq_len = len(sequence)

    print(f"\nPlasmid: {record.name}")
    print(f"Length: {seq_len:,} bp")
    print(f"Features: {len(record.features)}")

    # Convert GenBank features to format expected by parser
    features = []
    for feat in record.features:
        if feat.type == "source":
            continue

        name = feat.qualifiers.get("label", [feat.type])[0]
        start = int(feat.location.start)
        end = int(feat.location.end)
        strand = 1 if feat.location.strand == 1 else -1

        features.append({
            "name": name,
            "type": feat.type,
            "start": start,
            "end": end,
            "strand": strand,
            "description": feat.qualifiers.get("note", [""])[0] if "note" in feat.qualifiers else ""
        })

    print(f"\nFeatures by type:")
    feature_types = {}
    for f in features:
        ftype = f["type"]
        feature_types[ftype] = feature_types.get(ftype, 0) + 1
    for ftype, count in sorted(feature_types.items()):
        print(f"  {ftype}: {count}")

    # Initialize parser
    print("\n" + "="*80)
    print("Calling LLM Module Parser...")
    print("="*80)

    from splicify_api.llm_module_parser import LLMModuleParser

    try:
        parser = LLMModuleParser(model="gpt-4o")
        print("✓ Parser initialized")
    except ValueError as e:
        print(f"✗ Failed: {e}")
        return

    # Parse modules
    result = await parser.parse_modules(
        sequence=sequence,
        plannotate_features=features,
        cds_submodules=None,
        circular=True
    )

    modules = result.get("modules", [])

    print(f"\n✓ LLM identified {len(modules)} modules\n")

    # Display modules organized by hierarchy
    print("="*80)
    print("Module Hierarchy:")
    print("="*80)

    def print_module_tree(mod, level=0):
        indent = "  " * level
        symbol = "└─" if level > 0 else "●"

        mod_name = mod['module_type'].replace("_", " ").title()
        start = mod['start']
        end = mod['end']
        size = end - start

        # Format size
        if size > 1000:
            size_str = f"{size/1000:.1f}kb"
        else:
            size_str = f"{size}bp"

        print(f"{indent}{symbol} {mod_name}")
        print(f"{indent}   {start:,}..{end:,} ({size_str})")

        # Show key metadata
        metadata = mod.get('metadata', {})
        if metadata:
            if metadata.get('promoter_id'):
                print(f"{indent}   Promoter: {metadata['promoter_id']}")
            if metadata.get('polya_id'):
                print(f"{indent}   PolyA: {metadata['polya_id']}")
            if metadata.get('protein_name'):
                print(f"{indent}   Protein: {metadata['protein_name']}")
            if metadata.get('payload_class'):
                print(f"{indent}   Payload: {metadata['payload_class']}")

        # Print nested modules
        for child_id in mod.get('nested_modules', []):
            child = next((m for m in modules if m['module_id'] == child_id), None)
            if child:
                print_module_tree(child, level + 1)

    # Find top-level modules
    top_level = [m for m in modules if not m.get('parent_module')]

    for mod in top_level:
        print_module_tree(mod)
        print()

    # Analysis
    print("="*80)
    print("Biological Analysis:")
    print("="*80)

    # Check for expected modules
    checks = []

    # 1. Lentiviral payload
    lenti = [m for m in modules if 'lentiviral' in m['module_type'].lower()]
    if lenti:
        lenti_mod = lenti[0]
        lenti_size = lenti_mod['end'] - lenti_mod['start']
        checks.append(f"✓ Lentiviral payload: {lenti_mod['start']:,}..{lenti_mod['end']:,} ({lenti_size:,} bp)")
    else:
        checks.append("✗ Lentiviral payload NOT identified")

    # 2. Expression cassettes
    pol2_cassettes = [m for m in modules if 'pol2' in m['module_type'].lower() and 'expression' in m['module_type'].lower()]
    pol3_cassettes = [m for m in modules if 'pol3' in m['module_type'].lower() or 'guide' in m['module_type'].lower()]

    if pol2_cassettes:
        for pc in pol2_cassettes:
            prom = pc.get('metadata', {}).get('promoter_id', 'unknown')
            checks.append(f"✓ Pol II cassette: {prom} promoter")
    else:
        checks.append("✗ Pol II cassette NOT identified")

    if pol3_cassettes:
        for pc in pol3_cassettes:
            prom = pc.get('metadata', {}).get('promoter_id', 'unknown')
            checks.append(f"✓ Pol III/Guide cassette: {prom} promoter")
    else:
        checks.append("✗ Pol III/Guide cassette NOT identified")

    # 3. CDS modules
    cds_mods = [m for m in modules if 'cds' in m['module_type'].lower()]
    if cds_mods:
        checks.append(f"✓ CDS modules: {len(cds_mods)} identified")
    else:
        checks.append("✗ CDS modules NOT identified")

    # 4. Linker modules (2A peptides)
    linker_mods = [m for m in modules if 'linker' in m['module_type'].lower()]
    if linker_mods:
        checks.append(f"✓ Linker modules (2A): {len(linker_mods)} identified")
    else:
        checks.append("✗ Linker modules (2A) NOT identified")

    # 5. NLS modules
    nls_mods = [m for m in modules if 'nls' in m['module_type'].lower()]
    if nls_mods:
        checks.append(f"✓ NLS modules: {len(nls_mods)} identified")
    else:
        checks.append("✗ NLS modules NOT identified")

    # 6. Bacterial backbone
    backbone = [m for m in modules if 'backbone' in m['module_type'].lower()]
    if backbone:
        bb_mod = backbone[0]
        bb_size = bb_mod['end'] - bb_mod['start']
        checks.append(f"✓ Bacterial backbone: {bb_mod['start']:,}..{bb_mod['end']:,} ({bb_size:,} bp)")
    else:
        checks.append("✗ Bacterial backbone NOT identified")

    for check in checks:
        print(check)

    # Module type summary
    print("\n" + "="*80)
    print("Module Type Summary:")
    print("="*80)

    type_counts = {}
    for m in modules:
        mtype = m['module_type']
        type_counts[mtype] = type_counts.get(mtype, 0) + 1

    for mtype, count in sorted(type_counts.items()):
        print(f"  {mtype}: {count}")

    # Save results
    output_file = "/tmp/lenticrispr_v2_llm_modules.json"
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\n✓ Full results saved to: {output_file}")

    print("\n" + "="*80)
    print("Test Complete!")
    print("="*80)

    return result


if __name__ == "__main__":
    asyncio.run(test_real_plasmid())
