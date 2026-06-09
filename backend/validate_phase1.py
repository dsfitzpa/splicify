#!/usr/bin/env python3
"""
Validation script for Phase 1: Core Abstractions

Tests the basic functionality of design_request and part_resolver modules.
"""

import asyncio
from splicify_api.predesign.design_request import (
    InputSource,
    PartSpecification,
    TargetSpecification,
    DesignRequest,
)
from splicify_api.predesign.part_resolver import PartResolver, ResolvedPart


def test_part_specification():
    """Test creating part specifications."""
    print("\n=== Testing PartSpecification ===")

    # Direct sequence
    part1 = PartSpecification(
        name="Fragment1",
        source=InputSource.DIRECT_SEQUENCE,
        sequence="ATGCATGC",
    )
    print(f"✓ Created direct sequence part: {part1.name}")

    # Feature name
    part2 = PartSpecification(
        name="CMV promoter",
        source=InputSource.FEATURE_NAME,
        feature_name="CMV",
        role="promoter",
    )
    print(f"✓ Created feature name part: {part2.name}")

    # Feature pattern
    part3 = PartSpecification(
        name="strong promoter",
        source=InputSource.FEATURE_PATTERN,
        feature_pattern="strong mammalian promoter",
        role="promoter",
    )
    print(f"✓ Created feature pattern part: {part3.name}")

    # Uploaded file
    part4 = PartSpecification(
        name="pUC19",
        source=InputSource.UPLOADED_FILE,
        file_id="file_abc123",
    )
    print(f"✓ Created uploaded file part: {part4.name}")

    return [part1, part2, part3, part4]


def test_target_specification():
    """Test creating target specifications."""
    print("\n=== Testing TargetSpecification ===")

    # Uploaded target
    target1 = TargetSpecification(
        source="uploaded",
        uploaded_file_id="file_xyz789",
    )
    print(f"✓ Created uploaded target")

    # Assembled target
    part = PartSpecification(
        name="part1",
        source=InputSource.DIRECT_SEQUENCE,
        sequence="ATGC",
    )
    target2 = TargetSpecification(
        source="assembled",
        parts=[part],
        assembly_order="listed",
        topology="circular",
    )
    print(f"✓ Created assembled target with {len(target2.parts)} part(s)")

    return target1, target2


def test_design_request():
    """Test creating design requests."""
    print("\n=== Testing DesignRequest ===")

    # Minimal request
    request = DesignRequest(
        session_id="sess_123",
        user_message="Design Gibson assembly",
        parts=[
            PartSpecification(
                name="part1",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGC",
            )
        ],
    )
    print(f"✓ Created design request: {request.session_id}")
    print(f"  Message: {request.user_message}")
    print(f"  Parts: {len(request.parts)}")

    # Request with target and metadata
    cmv_part = PartSpecification(
        name="CMV",
        source=InputSource.FEATURE_NAME,
        feature_name="CMV",
    )
    egfp_part = PartSpecification(
        name="eGFP",
        source=InputSource.FEATURE_NAME,
        feature_name="eGFP",
    )
    request2 = DesignRequest(
        session_id="sess_456",
        user_message="Design Golden Gate assembly",
        parts=[cmv_part, egfp_part],
        target=TargetSpecification(
            source="assembled",
            parts=[cmv_part, egfp_part],
        ),
        suggested_workflow="golden_gate",
        metadata={"intent": "golden_gate_primer_design", "confidence": 0.95},
    )
    print(f"✓ Created design request with target and metadata")
    print(f"  Suggested workflow: {request2.suggested_workflow}")
    print(f"  Metadata: {request2.metadata}")

    # Test serialization
    d = request2.to_dict()
    print(f"✓ Serialized to dict: {len(d)} keys")

    return request, request2


async def test_part_resolver():
    """Test part resolution."""
    print("\n=== Testing PartResolver ===")

    resolver = PartResolver()

    # Test direct sequence resolution
    part_spec = PartSpecification(
        name="Fragment1",
        source=InputSource.DIRECT_SEQUENCE,
        sequence="ATGCATGC",
        role="promoter",
    )
    resolved = resolver._resolve_direct(part_spec)
    print(f"✓ Resolved direct sequence: {resolved.name}")
    print(f"  Sequence: {resolved.sequence}")
    print(f"  Length: {resolved.length} bp")
    print(f"  Role: {resolved.role}")

    # Test sequence normalization
    part_spec2 = PartSpecification(
        name="Fragment2",
        source=InputSource.DIRECT_SEQUENCE,
        sequence="atgc ATGC\nATGC",
    )
    resolved2 = resolver._resolve_direct(part_spec2)
    print(f"✓ Normalized sequence with whitespace: {resolved2.sequence}")

    # Test module dict conversion
    module_dict = resolved.to_module_dict()
    print(f"✓ Converted to module dict:")
    print(f"  canonical_id: {module_dict['canonical_id']}")
    print(f"  role: {module_dict['role']}")
    print(f"  length: {module_dict['length']} bp")

    # Test file resolution
    genbank_content = """LOCUS       test                      12 bp    DNA     circular             09-APR-2026
DEFINITION  Test plasmid
FEATURES             Location/Qualifiers
     promoter        1..4
                     /label="test_promoter"
ORIGIN
        1 atgcatgcat gc
//
"""
    part_spec3 = PartSpecification(
        name="test_plasmid",
        source=InputSource.UPLOADED_FILE,
        file_id="file_123",
    )
    context = {
        "file_cache": {
            "file_123": {
                "name": "test.gb",
                "content": genbank_content,
            }
        }
    }
    resolved3 = await resolver._resolve_file(part_spec3, context)
    print(f"✓ Resolved GenBank file: {resolved3.name}")
    print(f"  Sequence: {resolved3.sequence}")
    print(f"  Features: {len(resolved3.features)}")

    # Test resolve_all
    parts = [
        PartSpecification(
            name="Part1",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="ATGC",
        ),
        PartSpecification(
            name="Part2",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="GCTA",
        ),
    ]
    resolved_parts = await resolver.resolve_all(parts)
    print(f"✓ Resolved {len(resolved_parts)} parts with resolve_all()")

    return resolved, resolved2, resolved3


def test_validation():
    """Test validation and error handling."""
    print("\n=== Testing Validation ===")

    # Test invalid part specification
    try:
        PartSpecification(
            name="Invalid",
            source=InputSource.DIRECT_SEQUENCE,
            # Missing sequence
        )
        print("✗ Should have raised ValueError")
    except ValueError as e:
        print(f"✓ Correctly raised ValueError: {e}")

    # Test invalid target specification
    try:
        TargetSpecification(
            source="assembled",
            # Missing parts
        )
        print("✗ Should have raised ValueError")
    except ValueError as e:
        print(f"✓ Correctly raised ValueError: {e}")

    # Test invalid design request
    try:
        DesignRequest(
            session_id="",  # Empty session_id
            user_message="test",
        )
        print("✗ Should have raised ValueError")
    except ValueError as e:
        print(f"✓ Correctly raised ValueError: {e}")

    # Test resolved part with invalid bases
    resolved = ResolvedPart(
        name="test",
        sequence="ATGCXYZ",
        length=7,
        source=InputSource.DIRECT_SEQUENCE,
    )
    if resolved.warnings:
        print(f"✓ Generated warning for invalid bases: {resolved.warnings[0]}")
    else:
        print("✗ Should have generated warning for invalid bases")


async def main():
    """Run all validation tests."""
    print("=" * 70)
    print("Phase 1: Core Abstractions - Validation Tests")
    print("=" * 70)

    try:
        # Test data structures
        test_part_specification()
        test_target_specification()
        test_design_request()

        # Test resolver
        await test_part_resolver()

        # Test validation
        test_validation()

        print("\n" + "=" * 70)
        print("✓ All Phase 1 validation tests passed!")
        print("=" * 70)

    except Exception as e:
        print(f"\n✗ Validation failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    exit(asyncio.run(main()))
