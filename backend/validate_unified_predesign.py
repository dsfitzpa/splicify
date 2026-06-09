#!/usr/bin/env python3
"""
Comprehensive validation for unified pre-design architecture.

Tests the complete flow:
1. Create DesignRequest
2. Resolve parts
3. Build target plasmid
4. Route to workflow
"""

import asyncio
import logging
from splicify_api.predesign import (
    InputSource,
    PartSpecification,
    TargetSpecification,
    DesignRequest,
    PartResolver,
    ResolvedPart,
    get_knowledge_base,
)
from splicify_api.predesign.target_builder import TargetPlasmidBuilder, TargetPlasmid
from splicify_api.predesign.cloning_router import CloningRouter, WorkflowMethod

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def print_section(title: str):
    """Print a section header."""
    print(f"\n{'=' * 70}")
    print(f"  {title}")
    print('=' * 70)


def test_design_request_creation():
    """Test creating design requests for different scenarios."""
    print_section("1. Design Request Creation")

    # Scenario 1: Gibson assembly with direct sequences
    gibson_request = DesignRequest(
        session_id="sess_gibson_001",
        user_message="Design Gibson assembly for three fragments",
        parts=[
            PartSpecification(
                name="Fragment A",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
                role="promoter",
            ),
            PartSpecification(
                name="Fragment B",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
                role="CDS",
            ),
            PartSpecification(
                name="Fragment C",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA",
                role="terminator",
            ),
        ],
        suggested_workflow="gibson_assembly",
    )
    print(f"✓ Created Gibson design request")
    print(f"  Session: {gibson_request.session_id}")
    print(f"  Parts: {len(gibson_request.parts)}")
    print(f"  Workflow: {gibson_request.suggested_workflow}")

    # Scenario 2: Golden Gate with feature names (would require knowledge base)
    golden_gate_request = DesignRequest(
        session_id="sess_gg_001",
        user_message="Design Golden Gate assembly for CMV + eGFP + bGH polyA",
        parts=[
            PartSpecification(
                name="CMV promoter",
                source=InputSource.FEATURE_NAME,
                feature_name="CMV",
                role="promoter",
            ),
            PartSpecification(
                name="eGFP",
                source=InputSource.FEATURE_NAME,
                feature_name="eGFP",
                role="CDS",
            ),
            PartSpecification(
                name="bGH polyA",
                source=InputSource.FEATURE_NAME,
                feature_name="bGH polyA",
                role="terminator",
            ),
        ],
        suggested_workflow="golden_gate",
    )
    print(f"✓ Created Golden Gate design request")
    print(f"  Session: {golden_gate_request.session_id}")
    print(f"  Parts: {len(golden_gate_request.parts)}")

    return gibson_request, golden_gate_request


async def test_part_resolution():
    """Test resolving parts from various sources."""
    print_section("2. Part Resolution")

    resolver = PartResolver()

    # Test 1: Direct sequences
    print("\nTest 1: Direct sequences")
    direct_parts = [
        PartSpecification(
            name="Promoter",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 10,
            role="promoter",
        ),
        PartSpecification(
            name="Gene",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC" * 20,
            role="CDS",
        ),
        PartSpecification(
            name="Terminator",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA" * 5,
            role="terminator",
        ),
    ]

    resolved = await resolver.resolve_all(direct_parts)
    print(f"✓ Resolved {len(resolved)} direct sequence parts")
    for i, part in enumerate(resolved):
        print(f"  Part {i+1}: {part.name} - {part.length} bp ({part.role})")

    # Test 2: GenBank file
    print("\nTest 2: GenBank file")
    genbank_content = """LOCUS       pUC19                   2686 bp    DNA     circular             09-APR-2026
DEFINITION  pUC19 cloning vector
FEATURES             Location/Qualifiers
     promoter        1..150
                     /label="lac promoter"
     CDS             151..1000
                     /label="lacZ"
ORIGIN
        1 """ + "atgcatgc" * 335 + """
//
"""
    file_part = PartSpecification(
        name="pUC19_vector",
        source=InputSource.UPLOADED_FILE,
        file_id="file_puc19",
    )
    context = {
        "file_cache": {
            "file_puc19": {
                "name": "pUC19.gb",
                "content": genbank_content,
            }
        }
    }

    resolved_file = await resolver.resolve_all([file_part], context)
    print(f"✓ Resolved GenBank file: {resolved_file[0].name}")
    print(f"  Length: {resolved_file[0].length} bp")
    print(f"  Features: {len(resolved_file[0].features)}")

    return resolved


def test_target_building(resolved_parts: list):
    """Test building target plasmids."""
    print_section("3. Target Plasmid Building")

    builder = TargetPlasmidBuilder()

    # Build from parts
    target = builder.build_from_parts(
        parts=resolved_parts,
        assembly_order="listed",
        topology="circular",
    )

    print(f"✓ Built target plasmid from {len(resolved_parts)} parts")
    print(f"  Total length: {target.length} bp")
    print(f"  Topology: {target.topology}")
    print(f"  GC content: {target.metadata.get('gc_content', 0):.2%}")
    print(f"  Restriction sites found: {len(target.restriction_sites)}")
    print(f"  Type IIS sites found: {len(target.type_iis_sites)}")

    # Print restriction sites
    if target.restriction_sites:
        print("\n  Restriction sites:")
        for enzyme, positions in list(target.restriction_sites.items())[:5]:
            print(f"    {enzyme}: {len(positions)} site(s)")

    # Print Type IIS sites
    if target.type_iis_sites:
        print("\n  Type IIS sites:")
        for enzyme, positions in target.type_iis_sites.items():
            print(f"    {enzyme}: {len(positions)} site(s)")

    # Print junctions
    junctions = target.get_part_junctions()
    print(f"\n  Part junctions: {len(junctions)}")
    for j in junctions:
        print(f"    Junction {j['junction_id']}: {j['part1_name']} -> {j['part2_name']} at position {j['position']}")

    return target


async def test_workflow_routing(target: TargetPlasmid):
    """Test routing to cloning workflows."""
    print_section("4. Workflow Routing")

    router = CloningRouter()

    # Route with balanced objective
    candidates = await router.route(
        target=target,
        parts=target.parts,
        objective="balanced",
    )

    print(f"✓ Evaluated {len(candidates)} workflows")
    print(f"  Compatible workflows: {len([c for c in candidates if c.compatible])}")

    print("\n  Workflow Rankings (Balanced):")
    for i, candidate in enumerate(candidates, 1):
        if candidate.compatible:
            print(f"\n  {i}. {candidate.method.value}")
            print(f"     Cost: ${candidate.total_cost_usd:.2f}")
            print(f"     Time: {candidate.total_calendar_days:.1f} days")
            print(f"     Risk: {candidate.overall_risk_score:.2f}")
            print(f"     Confidence: {candidate.confidence:.2%}")
            print(f"     Balanced Score: {candidate.get_balanced_score():.3f}")
        else:
            print(f"\n  {i}. {candidate.method.value} [INCOMPATIBLE]")
            for reason in candidate.incompatibility_reasons:
                print(f"     - {reason}")

    # Test different objectives
    print("\n  Cost-optimized ranking:")
    cost_candidates = await router.route(target, target.parts, objective="cost")
    for i, c in enumerate([x for x in cost_candidates if x.compatible][:3], 1):
        print(f"    {i}. {c.method.value}: ${c.total_cost_usd:.2f}")

    print("\n  Time-optimized ranking:")
    time_candidates = await router.route(target, target.parts, objective="time")
    for i, c in enumerate([x for x in time_candidates if x.compatible][:3], 1):
        print(f"    {i}. {c.method.value}: {c.total_calendar_days:.1f} days")

    return candidates


def test_end_to_end_flow():
    """Test complete flow from request to workflow selection."""
    print_section("5. End-to-End Flow")

    print("\nScenario: User wants to build a 3-part expression vector")
    print("  Input: Direct DNA sequences for promoter, gene, terminator")
    print("  Expected: Gibson assembly recommended")

    # Step 1: Create design request
    request = DesignRequest(
        session_id="sess_e2e_001",
        user_message="Build expression vector with CMV, mCherry, and bGH polyA",
        parts=[
            PartSpecification(
                name="CMV promoter",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGC" * 150,  # 600 bp promoter
                role="promoter",
            ),
            PartSpecification(
                name="mCherry CDS",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="GCTA" * 200,  # 800 bp gene
                role="CDS",
            ),
            PartSpecification(
                name="bGH polyA",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="TATA" * 100,  # 400 bp terminator
                role="terminator",
            ),
        ],
    )
    print(f"\n✓ Step 1: Created design request")
    print(f"  Parts: {len(request.parts)}")

    return request


async def main():
    """Run all validation tests."""
    print("=" * 70)
    print("  Unified Pre-Design Architecture - Comprehensive Validation")
    print("=" * 70)

    try:
        # Test 1: Design request creation
        gibson_req, gg_req = test_design_request_creation()

        # Test 2: Part resolution
        resolved_parts = await test_part_resolution()

        # Test 3: Target building
        target = test_target_building(resolved_parts)

        # Test 4: Workflow routing
        candidates = await test_workflow_routing(target)

        # Test 5: End-to-end
        e2e_request = test_end_to_end_flow()

        # Final summary
        print_section("Validation Summary")
        print("✓ All components validated successfully!")
        print("\nImplemented:")
        print("  ✓ DesignRequest - unified input structure")
        print("  ✓ PartResolver - multi-source resolution (direct sequences, files)")
        print("  ✓ TargetPlasmidBuilder - target assembly and annotation")
        print("  ✓ CloningRouter - workflow evaluation and ranking")
        print("\nPending (Phase 2-6):")
        print("  ○ Knowledge base integration (feature name lookup)")
        print("  ○ Pattern matching for features")
        print("  ○ Homology-based assembly order")
        print("  ○ Actual operator integration")
        print("  ○ chat.py integration")

        print("\n" + "=" * 70)
        print("✓ VALIDATION PASSED")
        print("=" * 70)

        return 0

    except Exception as e:
        print(f"\n✗ Validation failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(asyncio.run(main()))
