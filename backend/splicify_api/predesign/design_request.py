"""
Unified input structure for cloning workflow design requests.

This module defines the standard representation of user intent before resolution
and execution. All cloning workflows (Gibson, Golden Gate, SDM, Restriction, Gateway)
use this common structure.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, List, Dict, Any


class InputSource(Enum):
    """Source type for part specifications."""

    DIRECT_SEQUENCE = "direct_sequence"
    UPLOADED_FILE = "uploaded_file"
    FEATURE_NAME = "feature_name"
    FEATURE_PATTERN = "feature_pattern"
    HOMOLOGY_DERIVED = "homology_derived"


@dataclass
class PartSpecification:
    """
    Specification for a single part before resolution.

    Attributes:
        name: User-provided or inferred name for the part
        source: How the part was specified (sequence, file, feature name, etc.)
        sequence: Direct DNA sequence (if source=DIRECT_SEQUENCE)
        feature_name: Exact feature name to lookup (if source=FEATURE_NAME)
        feature_pattern: Pattern/description to search (if source=FEATURE_PATTERN)
        file_id: Uploaded file identifier (if source=UPLOADED_FILE)
        role: Biological role (promoter, CDS, terminator, etc.)
        specified_order: User-specified assembly order (1-indexed, optional)

    Examples:
        # Direct sequence
        PartSpecification(
            name="Fragment1",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="ATGCATGC..."
        )

        # Feature name lookup
        PartSpecification(
            name="CMV promoter",
            source=InputSource.FEATURE_NAME,
            feature_name="CMV",
            role="promoter"
        )

        # Pattern search
        PartSpecification(
            name="strong mammalian promoter",
            source=InputSource.FEATURE_PATTERN,
            feature_pattern="strong mammalian promoter",
            role="promoter"
        )
    """

    name: str
    source: InputSource
    sequence: Optional[str] = None
    feature_name: Optional[str] = None
    feature_pattern: Optional[str] = None
    file_id: Optional[str] = None
    role: Optional[str] = None
    specified_order: Optional[int] = None

    def __post_init__(self):
        """Validate that required fields are present for the source type."""
        if self.source == InputSource.DIRECT_SEQUENCE and not self.sequence:
            raise ValueError(f"Part '{self.name}': DIRECT_SEQUENCE requires sequence field")

        if self.source == InputSource.FEATURE_NAME and not self.feature_name:
            raise ValueError(f"Part '{self.name}': FEATURE_NAME requires feature_name field")

        if self.source == InputSource.FEATURE_PATTERN and not self.feature_pattern:
            raise ValueError(f"Part '{self.name}': FEATURE_PATTERN requires feature_pattern field")

        if self.source == InputSource.UPLOADED_FILE and not self.file_id:
            raise ValueError(f"Part '{self.name}': UPLOADED_FILE requires file_id field")


@dataclass
class TargetSpecification:
    """
    Specification for the target plasmid.

    Attributes:
        source: How the target is defined
            - "uploaded": User uploaded a target plasmid file
            - "assembled": Target is assembled from parts
            - "derived": Target is derived from anchor plasmid + modifications
        uploaded_file_id: File identifier if source="uploaded"
        parts: Parts to assemble (if source="assembled")
        assembly_order: How to order parts
            - "listed": Use the order specified in parts list
            - "homology_based": Detect order from homology regions
        topology: "circular" or "linear"

    Examples:
        # Target from uploaded file
        TargetSpecification(
            source="uploaded",
            uploaded_file_id="pUC19_abc123"
        )

        # Target from assembled parts
        TargetSpecification(
            source="assembled",
            parts=[part1, part2, part3],
            assembly_order="listed",
            topology="circular"
        )
    """

    source: str  # "uploaded" | "assembled" | "derived"
    uploaded_file_id: Optional[str] = None
    parts: List[PartSpecification] = field(default_factory=list)
    assembly_order: str = "listed"  # "listed" | "homology_based"
    topology: str = "circular"  # "circular" | "linear"

    def __post_init__(self):
        """Validate target specification."""
        valid_sources = {"uploaded", "assembled", "derived"}
        if self.source not in valid_sources:
            raise ValueError(f"Target source must be one of {valid_sources}, got '{self.source}'")

        if self.source == "uploaded" and not self.uploaded_file_id:
            raise ValueError("Target source 'uploaded' requires uploaded_file_id")

        if self.source == "assembled" and not self.parts:
            raise ValueError("Target source 'assembled' requires parts list")

        valid_orders = {"listed", "homology_based"}
        if self.assembly_order not in valid_orders:
            raise ValueError(f"Assembly order must be one of {valid_orders}, got '{self.assembly_order}'")

        valid_topologies = {"circular", "linear"}
        if self.topology not in valid_topologies:
            raise ValueError(f"Topology must be one of {valid_topologies}, got '{self.topology}'")


@dataclass
class DesignRequest:
    """
    Unified representation of a cloning design request.

    This structure standardizes input across all cloning workflows and sits between
    intent parsing and operator selection. It captures what the user wants to build
    before resolution and execution.

    Attributes:
        session_id: Unique session identifier
        user_message: Original user message (for context)
        parts: Parts specified by the user (may need resolution)
        target: Target plasmid specification (optional, can be derived from parts)
        inventory_file_ids: Uploaded inventory files for homology-based assembly
        suggested_workflow: User-suggested workflow name (optional override)
        metadata: Additional context (intent classification, confidence scores, etc.)

    Workflow:
        1. Parse user intent → create DesignRequest
        2. Resolve parts → convert PartSpecification to ResolvedPart
        3. Build target → create TargetPlasmid from resolved parts
        4. Route → select optimal workflow (Gibson, Golden Gate, etc.)
        5. Execute → delegate to selected operator

    Example:
        # Gibson assembly with feature names
        DesignRequest(
            session_id="sess_123",
            user_message="Design Gibson assembly for CMV + eGFP + bGH polyA",
            parts=[
                PartSpecification(name="CMV", source=InputSource.FEATURE_NAME, feature_name="CMV"),
                PartSpecification(name="eGFP", source=InputSource.FEATURE_NAME, feature_name="eGFP"),
                PartSpecification(name="bGH polyA", source=InputSource.FEATURE_NAME, feature_name="bGH polyA"),
            ],
            target=TargetSpecification(source="assembled", parts=[...], topology="circular"),
            suggested_workflow="gibson_assembly"
        )
    """

    session_id: str
    user_message: str
    parts: List[PartSpecification] = field(default_factory=list)
    target: Optional[TargetSpecification] = None
    inventory_file_ids: List[str] = field(default_factory=list)
    suggested_workflow: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate design request."""
        if not self.session_id:
            raise ValueError("session_id is required")

        if not self.user_message:
            raise ValueError("user_message is required")

        # At least one of parts or target.uploaded_file_id must be present
        has_parts = len(self.parts) > 0
        has_uploaded_target = (
            self.target is not None
            and self.target.source == "uploaded"
            and self.target.uploaded_file_id
        )

        if not has_parts and not has_uploaded_target:
            raise ValueError("DesignRequest requires either parts or an uploaded target file")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "session_id": self.session_id,
            "user_message": self.user_message,
            "parts": [
                {
                    "name": p.name,
                    "source": p.source.value,
                    "sequence": p.sequence,
                    "feature_name": p.feature_name,
                    "feature_pattern": p.feature_pattern,
                    "file_id": p.file_id,
                    "role": p.role,
                    "specified_order": p.specified_order,
                }
                for p in self.parts
            ],
            "target": {
                "source": self.target.source,
                "uploaded_file_id": self.target.uploaded_file_id,
                "assembly_order": self.target.assembly_order,
                "topology": self.target.topology,
            } if self.target else None,
            "inventory_file_ids": self.inventory_file_ids,
            "suggested_workflow": self.suggested_workflow,
            "metadata": self.metadata,
        }
