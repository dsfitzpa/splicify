"""
Unified pre-design phase for plasmid cloning workflows.

This package provides a standardized layer between intent parsing and operator selection
that handles:
- Input resolution (sequences, files, feature names, patterns)
- Target plasmid building
- Workflow routing and selection

Components:
- design_request: Unified input structure (DesignRequest, PartSpecification, etc.)
- part_resolver: Multi-source part resolution (PartResolver)
- target_builder: Target plasmid construction (TargetPlasmidBuilder)
- cloning_router: Workflow selection and ranking (CloningRouter)
- knowledge_base: pLannotate feature lookup (PlannotateKnowledgeBase)
"""

from .design_request import (
    InputSource,
    PartSpecification,
    TargetSpecification,
    DesignRequest,
)
from .part_resolver import (
    PartResolver,
    ResolvedPart,
)
from .knowledge_base import (
    PlannotateKnowledgeBase,
    get_knowledge_base,
)
from .target_builder import (
    TargetPlasmid,
    TargetPlasmidBuilder,
)
from .cloning_router import (
    CloningRouter,
    WorkflowMethod,
    WorkflowCandidate,
)

__all__ = [
    "InputSource",
    "PartSpecification",
    "TargetSpecification",
    "DesignRequest",
    "PartResolver",
    "ResolvedPart",
    "PlannotateKnowledgeBase",
    "get_knowledge_base",
    "TargetPlasmid",
    "TargetPlasmidBuilder",
    "CloningRouter",
    "WorkflowMethod",
    "WorkflowCandidate",
]
