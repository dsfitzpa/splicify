#!/usr/bin/env python3
"""
Demonstration of unified pre-design architecture.

This script shows how to use the pre-design components for different
cloning scenarios.
"""

import asyncio
from splicify_api.predesign import (
    DesignRequest,
    PartSpecification,
    TargetSpecification,
    InputSource,
    PartResolver,
    TargetPlasmidBuilder,
    CloningRouter,
    WorkflowMethod,
)


async def demo_gibson_assembly():
    """Demonstrate Gibson assembly workflow."""
    print("\n" + "=" * 70)
    print("DEMO 1: Gibson Assembly - 3 Fragment Cloning")
    print("=" * 70)

    # User wants to build an expression vector
    print("\nUser request: 'Build expression vector with promoter, gene, and terminator'")

    # Step 1: Create design request
    request = DesignRequest(
        session_id="demo_gibson",
        user_message="Build expression vector with CMV, mCherry, and bGH polyA",
        parts=[
            PartSpecification(
                name="CMV promoter",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC" * 20,  # 720 bp
                role="promoter",
            ),
            PartSpecification(
                name="mCherry CDS",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA" * 25,  # 900 bp
                role="CDS",
            ),
            PartSpecification(
                name="bGH polyA",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="TATATATATATATATATATATATATATATATATATA" * 10,  # 360 bp
                role="terminator",
            ),
        ],
    )

    print(f"\n✓ Design request created")
    print(f"  Session: {request.session_id}")
    print(f"  Parts: {len(request.parts)}")

    # Step 2: Resolve parts
    resolver = PartResolver()
    resolved = await resolver.resolve_all(request.parts)

    print(f"\n✓ Parts resolved")
    for i, part in enumerate(resolved, 1):
        print(f"  {i}. {part.name}: {part.length} bp ({part.role})")

    # Step 3: Build target plasmid
    builder = TargetPlasmidBuilder()
    target = builder.build_from_parts(
        parts=resolved,
        assembly_order="listed",
        topology="circular",
    )

    print(f"\n✓ Target plasmid built")
    print(f"  Total size: {target.length} bp")
    print(f"  Topology: {target.topology}")
    print(f"  GC content: {target.metadata['gc_content']:.1%}")

    # Print junctions
    junctions = target.get_part_junctions()
    print(f"  Part junctions: {len(junctions)}")
    for j in junctions[:3]:  # Show first 3
        print(f"    {j['part1_name']} → {j['part2_name']} at position {j['position']}")

    # Step 4: Route to workflows
    router = CloningRouter()
    candidates = await router.route(
        target=target,
        parts=resolved,
        objective="balanced",
    )

    print(f"\n✓ Workflows evaluated")
    print(f"  Compatible: {len([c for c in candidates if c.compatible])}")
    print(f"  Incompatible: {len([c for c in candidates if not c.compatible])}")

    # Show top recommendation
    best = candidates[0]
    print(f"\n🏆 RECOMMENDED WORKFLOW: {best.method.value}")
    print(f"   Estimated cost: ${best.total_cost_usd:.2f}")
    print(f"   Estimated time: {best.total_calendar_days:.0f} days")
    print(f"   Risk score: {best.overall_risk_score:.2f} (0=low, 1=high)")
    print(f"   Confidence: {best.confidence:.0%}")

    # Show alternatives
    print(f"\n   Alternatives:")
    for candidate in candidates[1:]:
        if candidate.compatible:
            print(f"   - {candidate.method.value}: ${candidate.total_cost_usd:.0f}, "
                  f"{candidate.total_calendar_days:.0f} days")
        else:
            print(f"   - {candidate.method.value}: INCOMPATIBLE")


async def demo_golden_gate_assembly():
    """Demonstrate Golden Gate assembly workflow."""
    print("\n" + "=" * 70)
    print("DEMO 2: Golden Gate Assembly - Multi-Fragment Assembly")
    print("=" * 70)

    print("\nUser request: 'Use Golden Gate to assemble 5 fragments'")

    # Create request with 5 fragments
    parts = []
    for i in range(5):
        parts.append(
            PartSpecification(
                name=f"Fragment {i+1}",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGCATGCATGCATGC" * 15,  # 240 bp each
                role="fragment",
            )
        )

    request = DesignRequest(
        session_id="demo_golden_gate",
        user_message="Golden Gate assembly of 5 fragments",
        parts=parts,
    )

    print(f"\n✓ Design request created with {len(request.parts)} parts")

    # Resolve and build
    resolver = PartResolver()
    resolved = await resolver.resolve_all(request.parts)

    builder = TargetPlasmidBuilder()
    target = builder.build_from_parts(resolved, topology="circular")

    print(f"\n✓ Target plasmid: {target.length} bp")

    # Route
    router = CloningRouter()
    candidates = await router.route(target, resolved, objective="time")

    print(f"\n✓ Workflows evaluated (optimized for time)")

    # Show time-optimized ranking
    for i, candidate in enumerate(candidates[:3], 1):
        if candidate.compatible:
            print(f"  {i}. {candidate.method.value}: {candidate.total_calendar_days:.0f} days "
                  f"(${candidate.total_cost_usd:.0f})")


async def demo_uploaded_file():
    """Demonstrate workflow with uploaded plasmid."""
    print("\n" + "=" * 70)
    print("DEMO 3: Uploaded Plasmid Analysis")
    print("=" * 70)

    print("\nUser request: 'Analyze this plasmid and suggest modifications'")

    # Simulate uploaded GenBank file
    genbank_content = """LOCUS       pUC19                   2686 bp    DNA     circular             09-APR-2026
DEFINITION  pUC19 cloning vector
FEATURES             Location/Qualifiers
     promoter        1..150
                     /label="lac promoter"
     CDS             151..1000
                     /label="lacZ alpha"
     misc_feature    1001..1050
                     /label="MCS"
ORIGIN
        1 """ + "atgcatgc" * 335 + """
//
"""

    request = DesignRequest(
        session_id="demo_uploaded",
        user_message="Analyze uploaded plasmid",
        parts=[
            PartSpecification(
                name="pUC19",
                source=InputSource.UPLOADED_FILE,
                file_id="file_puc19",
            )
        ],
    )

    print(f"\n✓ Design request created for uploaded file")

    # Resolve with file cache
    resolver = PartResolver()
    context = {
        "file_cache": {
            "file_puc19": {
                "name": "pUC19.gb",
                "content": genbank_content,
            }
        }
    }

    resolved = await resolver.resolve_all(request.parts, context)

    print(f"\n✓ File parsed")
    print(f"  Sequence length: {resolved[0].length} bp")
    print(f"  Features found: {len(resolved[0].features)}")

    # Build target
    builder = TargetPlasmidBuilder()
    target = builder.build_from_upload(
        sequence=resolved[0].sequence,
        features=resolved[0].features,
        topology="circular",
    )

    print(f"\n✓ Target analyzed")
    print(f"  Restriction sites: {len(target.restriction_sites)}")
    print(f"  Type IIS sites: {len(target.type_iis_sites)}")

    # Show some restriction sites
    if target.restriction_sites:
        print(f"  Available for cloning:")
        for enzyme, positions in list(target.restriction_sites.items())[:5]:
            print(f"    {enzyme}: {len(positions)} site(s)")


async def demo_compatibility_checking():
    """Demonstrate compatibility checking."""
    print("\n" + "=" * 70)
    print("DEMO 4: Workflow Compatibility Checking")
    print("=" * 70)

    # Test different scenarios
    scenarios = [
        {
            "name": "2 fragments (ideal for all)",
            "num_parts": 2,
            "part_size": 500,
        },
        {
            "name": "7 fragments (too many for Gibson)",
            "num_parts": 7,
            "part_size": 300,
        },
        {
            "name": "50 bp fragments (too short for Gibson)",
            "num_parts": 3,
            "part_size": 50,
        },
        {
            "name": "40 fragments (too many for Golden Gate)",
            "num_parts": 40,
            "part_size": 200,
        },
    ]

    router = CloningRouter()

    for scenario in scenarios:
        print(f"\n{scenario['name']}:")

        # Create parts
        parts = []
        for i in range(scenario["num_parts"]):
            parts.append(
                PartSpecification(
                    name=f"Part{i+1}",
                    source=InputSource.DIRECT_SEQUENCE,
                    sequence="ATGC" * (scenario["part_size"] // 4),
                )
            )

        # Resolve and build
        resolver = PartResolver()
        resolved = await resolver.resolve_all(parts)

        builder = TargetPlasmidBuilder()
        target = builder.build_from_parts(resolved, topology="circular")

        # Check compatibility
        candidates = await router.route(target, resolved)

        # Count compatible
        compatible = [c for c in candidates if c.compatible]
        print(f"  Compatible workflows: {len(compatible)}/{len(candidates)}")

        for c in candidates:
            if c.compatible:
                print(f"    ✓ {c.method.value}")
            else:
                print(f"    ✗ {c.method.value}: {c.incompatibility_reasons[0]}")


async def main():
    """Run all demonstrations."""
    print("\n" + "=" * 70)
    print("  UNIFIED PRE-DESIGN ARCHITECTURE - DEMONSTRATION")
    print("=" * 70)

    await demo_gibson_assembly()
    await demo_golden_gate_assembly()
    await demo_uploaded_file()
    await demo_compatibility_checking()

    print("\n" + "=" * 70)
    print("  DEMONSTRATION COMPLETE")
    print("=" * 70)
    print("\nKey Features Demonstrated:")
    print("  ✓ Multi-source part resolution (sequences, files)")
    print("  ✓ Target plasmid building and annotation")
    print("  ✓ Workflow compatibility checking")
    print("  ✓ Multi-objective optimization (cost, time, risk)")
    print("  ✓ Intelligent workflow recommendations")
    print("\nNext Steps:")
    print("  → Integrate with chat.py (Phase 5)")
    print("  → Connect to actual operators for real metrics")
    print("  → Deploy pLannotate knowledge bases")
    print("  → Add homology-based assembly (Phase 6)")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    asyncio.run(main())
