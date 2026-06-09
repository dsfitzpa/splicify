"""
Unit tests for design_request module.

Tests the core data structures: InputSource, PartSpecification,
TargetSpecification, and DesignRequest.
"""

import pytest
from splicify_api.predesign.design_request import (
    InputSource,
    PartSpecification,
    TargetSpecification,
    DesignRequest,
)


class TestInputSource:
    """Test InputSource enum."""

    def test_all_sources_defined(self):
        """Test that all expected sources are defined."""
        expected = {
            "DIRECT_SEQUENCE",
            "UPLOADED_FILE",
            "FEATURE_NAME",
            "FEATURE_PATTERN",
            "HOMOLOGY_DERIVED",
        }
        actual = {source.name for source in InputSource}
        assert actual == expected


class TestPartSpecification:
    """Test PartSpecification dataclass."""

    def test_direct_sequence_valid(self):
        """Test creating a part with direct sequence."""
        part = PartSpecification(
            name="Fragment1",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="ATGCATGC",
        )
        assert part.name == "Fragment1"
        assert part.source == InputSource.DIRECT_SEQUENCE
        assert part.sequence == "ATGCATGC"

    def test_direct_sequence_missing_sequence(self):
        """Test that DIRECT_SEQUENCE requires sequence field."""
        with pytest.raises(ValueError, match="requires sequence field"):
            PartSpecification(
                name="Fragment1",
                source=InputSource.DIRECT_SEQUENCE,
            )

    def test_feature_name_valid(self):
        """Test creating a part with feature name."""
        part = PartSpecification(
            name="CMV promoter",
            source=InputSource.FEATURE_NAME,
            feature_name="CMV",
            role="promoter",
        )
        assert part.name == "CMV promoter"
        assert part.source == InputSource.FEATURE_NAME
        assert part.feature_name == "CMV"
        assert part.role == "promoter"

    def test_feature_name_missing_feature_name(self):
        """Test that FEATURE_NAME requires feature_name field."""
        with pytest.raises(ValueError, match="requires feature_name field"):
            PartSpecification(
                name="CMV",
                source=InputSource.FEATURE_NAME,
            )

    def test_feature_pattern_valid(self):
        """Test creating a part with feature pattern."""
        part = PartSpecification(
            name="strong promoter",
            source=InputSource.FEATURE_PATTERN,
            feature_pattern="strong mammalian promoter",
            role="promoter",
        )
        assert part.name == "strong promoter"
        assert part.source == InputSource.FEATURE_PATTERN
        assert part.feature_pattern == "strong mammalian promoter"

    def test_feature_pattern_missing_pattern(self):
        """Test that FEATURE_PATTERN requires feature_pattern field."""
        with pytest.raises(ValueError, match="requires feature_pattern field"):
            PartSpecification(
                name="promoter",
                source=InputSource.FEATURE_PATTERN,
            )

    def test_uploaded_file_valid(self):
        """Test creating a part from uploaded file."""
        part = PartSpecification(
            name="pUC19",
            source=InputSource.UPLOADED_FILE,
            file_id="file_abc123",
        )
        assert part.name == "pUC19"
        assert part.source == InputSource.UPLOADED_FILE
        assert part.file_id == "file_abc123"

    def test_uploaded_file_missing_file_id(self):
        """Test that UPLOADED_FILE requires file_id field."""
        with pytest.raises(ValueError, match="requires file_id field"):
            PartSpecification(
                name="plasmid",
                source=InputSource.UPLOADED_FILE,
            )

    def test_specified_order(self):
        """Test that specified_order is preserved."""
        part = PartSpecification(
            name="Fragment2",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="ATGC",
            specified_order=2,
        )
        assert part.specified_order == 2


class TestTargetSpecification:
    """Test TargetSpecification dataclass."""

    def test_uploaded_target_valid(self):
        """Test creating target from uploaded file."""
        target = TargetSpecification(
            source="uploaded",
            uploaded_file_id="file_xyz789",
        )
        assert target.source == "uploaded"
        assert target.uploaded_file_id == "file_xyz789"
        assert target.topology == "circular"

    def test_uploaded_target_missing_file_id(self):
        """Test that uploaded target requires file_id."""
        with pytest.raises(ValueError, match="requires uploaded_file_id"):
            TargetSpecification(
                source="uploaded",
            )

    def test_assembled_target_valid(self):
        """Test creating target from assembled parts."""
        parts = [
            PartSpecification(
                name="part1",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGC",
            )
        ]
        target = TargetSpecification(
            source="assembled",
            parts=parts,
            assembly_order="listed",
            topology="circular",
        )
        assert target.source == "assembled"
        assert len(target.parts) == 1
        assert target.assembly_order == "listed"
        assert target.topology == "circular"

    def test_assembled_target_missing_parts(self):
        """Test that assembled target requires parts."""
        with pytest.raises(ValueError, match="requires parts list"):
            TargetSpecification(
                source="assembled",
            )

    def test_invalid_source(self):
        """Test that invalid source raises error."""
        with pytest.raises(ValueError, match="must be one of"):
            TargetSpecification(
                source="invalid",
            )

    def test_invalid_assembly_order(self):
        """Test that invalid assembly_order raises error."""
        parts = [
            PartSpecification(
                name="part1",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGC",
            )
        ]
        with pytest.raises(ValueError, match="Assembly order must be"):
            TargetSpecification(
                source="assembled",
                parts=parts,
                assembly_order="invalid",
            )

    def test_invalid_topology(self):
        """Test that invalid topology raises error."""
        parts = [
            PartSpecification(
                name="part1",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGC",
            )
        ]
        with pytest.raises(ValueError, match="Topology must be"):
            TargetSpecification(
                source="assembled",
                parts=parts,
                topology="invalid",
            )

    def test_linear_topology(self):
        """Test linear topology."""
        parts = [
            PartSpecification(
                name="part1",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGC",
            )
        ]
        target = TargetSpecification(
            source="assembled",
            parts=parts,
            topology="linear",
        )
        assert target.topology == "linear"


class TestDesignRequest:
    """Test DesignRequest dataclass."""

    def test_minimal_valid_request(self):
        """Test creating a minimal valid design request."""
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
        assert request.session_id == "sess_123"
        assert request.user_message == "Design Gibson assembly"
        assert len(request.parts) == 1

    def test_missing_session_id(self):
        """Test that session_id is required."""
        with pytest.raises(ValueError, match="session_id is required"):
            DesignRequest(
                session_id="",
                user_message="test",
                parts=[
                    PartSpecification(
                        name="part1",
                        source=InputSource.DIRECT_SEQUENCE,
                        sequence="ATGC",
                    )
                ],
            )

    def test_missing_user_message(self):
        """Test that user_message is required."""
        with pytest.raises(ValueError, match="user_message is required"):
            DesignRequest(
                session_id="sess_123",
                user_message="",
                parts=[
                    PartSpecification(
                        name="part1",
                        source=InputSource.DIRECT_SEQUENCE,
                        sequence="ATGC",
                    )
                ],
            )

    def test_requires_parts_or_uploaded_target(self):
        """Test that either parts or uploaded target is required."""
        with pytest.raises(ValueError, match="requires either parts or an uploaded target"):
            DesignRequest(
                session_id="sess_123",
                user_message="test",
            )

    def test_with_uploaded_target(self):
        """Test request with uploaded target (no parts needed)."""
        request = DesignRequest(
            session_id="sess_123",
            user_message="Analyze this plasmid",
            target=TargetSpecification(
                source="uploaded",
                uploaded_file_id="file_abc",
            ),
        )
        assert request.target is not None
        assert request.target.source == "uploaded"
        assert len(request.parts) == 0

    def test_with_inventory_files(self):
        """Test request with inventory files."""
        request = DesignRequest(
            session_id="sess_123",
            user_message="Design from inventory",
            parts=[
                PartSpecification(
                    name="part1",
                    source=InputSource.DIRECT_SEQUENCE,
                    sequence="ATGC",
                )
            ],
            inventory_file_ids=["inv_1", "inv_2"],
        )
        assert len(request.inventory_file_ids) == 2

    def test_with_suggested_workflow(self):
        """Test request with suggested workflow."""
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
            suggested_workflow="gibson_assembly",
        )
        assert request.suggested_workflow == "gibson_assembly"

    def test_with_metadata(self):
        """Test request with metadata."""
        request = DesignRequest(
            session_id="sess_123",
            user_message="test",
            parts=[
                PartSpecification(
                    name="part1",
                    source=InputSource.DIRECT_SEQUENCE,
                    sequence="ATGC",
                )
            ],
            metadata={"intent": "gibson_design", "confidence": 0.95},
        )
        assert request.metadata["intent"] == "gibson_design"
        assert request.metadata["confidence"] == 0.95

    def test_to_dict(self):
        """Test serialization to dict."""
        request = DesignRequest(
            session_id="sess_123",
            user_message="test",
            parts=[
                PartSpecification(
                    name="part1",
                    source=InputSource.DIRECT_SEQUENCE,
                    sequence="ATGC",
                )
            ],
        )
        d = request.to_dict()
        assert d["session_id"] == "sess_123"
        assert d["user_message"] == "test"
        assert len(d["parts"]) == 1
        assert d["parts"][0]["name"] == "part1"
        assert d["parts"][0]["source"] == "direct_sequence"

    def test_to_dict_with_target(self):
        """Test serialization with target."""
        parts = [
            PartSpecification(
                name="part1",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGC",
            )
        ]
        request = DesignRequest(
            session_id="sess_123",
            user_message="test",
            parts=parts,
            target=TargetSpecification(
                source="assembled",
                parts=parts,
            ),
        )
        d = request.to_dict()
        assert d["target"] is not None
        assert d["target"]["source"] == "assembled"
        assert d["target"]["topology"] == "circular"
