"""
Cloning workflow router and selector.

This module evaluates all compatible cloning workflows (Gibson, Golden Gate,
Restriction, SDM) and ranks them by objective (cost, time, risk, balanced).
"""

import logging
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from enum import Enum

from .target_builder import TargetPlasmid
from .part_resolver import ResolvedPart

logger = logging.getLogger(__name__)


class WorkflowMethod(Enum):
    """Supported cloning workflow methods."""

    GIBSON = "gibson_assembly"
    GOLDEN_GATE = "golden_gate"
    RESTRICTION = "restriction_cloning"
    SDM = "site_directed_mutagenesis"
    GATEWAY = "gateway_cloning"
    INVENTORY_GIBSON = "inventory_gibson"  # Gibson assembly from existing plasmid inventory


@dataclass
class WorkflowCandidate:
    """
    A candidate cloning workflow with evaluation metrics.

    Attributes:
        method: Workflow method name
        compatible: Whether the workflow is compatible with the target
        build_plan: BuildPlan from operator evaluation (if compatible)
        total_cost_usd: Total estimated cost in USD
        total_calendar_days: Estimated calendar days to completion
        overall_risk_score: Overall risk score (0.0-1.0, higher = riskier)
        confidence: Confidence in success (0.0-1.0)
        warnings: List of warning messages
        incompatibility_reasons: Reasons why workflow is incompatible
        metadata: Additional evaluation metadata
    """

    method: WorkflowMethod
    compatible: bool
    build_plan: Optional[Dict[str, Any]] = None
    total_cost_usd: float = 0.0
    total_calendar_days: float = 0.0
    overall_risk_score: float = 0.0
    confidence: float = 1.0
    warnings: List[str] = field(default_factory=list)
    incompatibility_reasons: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def get_balanced_score(self) -> float:
        """
        Calculate balanced score across all metrics.

        Lower is better. Weights:
        - Cost: 0.3
        - Time: 0.3
        - Risk: 0.4

        Returns:
            Normalized score (0.0-1.0)
        """
        if not self.compatible:
            return 1.0  # Maximum penalty for incompatible

        # Normalize metrics to 0-1 scale
        # These thresholds are based on typical lab costs and timelines
        cost_norm = min(self.total_cost_usd / 500.0, 1.0)  # $500 = max
        time_norm = min(self.total_calendar_days / 14.0, 1.0)  # 14 days = max
        risk_norm = self.overall_risk_score  # Already 0-1

        # Weighted combination
        score = (
            0.3 * cost_norm +
            0.3 * time_norm +
            0.4 * risk_norm
        )

        # Penalty for low confidence
        score *= (1.0 / max(self.confidence, 0.1))

        return score


class CloningRouter:
    """
    Route design requests to optimal cloning workflows.

    This class:
    1. Evaluates workflow compatibility (fast checks)
    2. Calls operators for full evaluation (only for compatible workflows)
    3. Ranks candidates by objective
    4. Returns ordered list of workflow options
    """

    def __init__(self):
        """Initialize cloning router."""
        # Workflow constraints (for fast compatibility checking)
        self.constraints = {
            WorkflowMethod.GIBSON: {
                "min_fragments": 2,
                "max_fragments": 6,
                "min_fragment_length": 100,
                "min_overlap": 15,
                "max_overlap": 50,
            },
            WorkflowMethod.GOLDEN_GATE: {
                "min_fragments": 2,
                "max_fragments": 32,
                "requires_type_iis": True,
            },
            WorkflowMethod.RESTRICTION: {
                "min_fragments": 1,
                "max_fragments": 3,
                "requires_restriction_sites": True,
            },
            WorkflowMethod.SDM: {
                "requires_anchor": True,
                "max_edits": 5,
            },
            WorkflowMethod.GATEWAY: {
                "min_fragments": 2,
                "max_fragments": 8,
            },
        }

    async def route(
        self,
        target: TargetPlasmid,
        parts: Optional[List[ResolvedPart]] = None,
        anchor_plasmid: Optional[TargetPlasmid] = None,
        objective: str = "balanced",
        preferred_method: Optional[WorkflowMethod] = None,
    ) -> List[WorkflowCandidate]:
        """
        Evaluate and rank all compatible workflows.

        Args:
            target: Target plasmid to build
            parts: Resolved parts (if assembly-based)
            anchor_plasmid: Anchor plasmid (for SDM)
            objective: Optimization objective:
                - "balanced": Weighted score across all metrics
            preferred_method: If provided, strongly prefer this workflow method
                - "cost": Minimize cost
                - "time": Minimize time
                - "risk": Minimize risk

        Returns:
            List of workflow candidates, ranked by objective (best first)
        """
        logger.info(f"Routing workflow for target ({target.length} bp), objective={objective}, preferred_method={preferred_method}")

        # Annotate cloning features
        annotations = self._annotate_cloning_features(target, parts or [])

        # Fast compatibility checks
        compatibility = self._check_workflow_compatibility(annotations, anchor_plasmid)

        # Evaluate compatible workflows
        candidates = []
        for method, compat_info in compatibility.items():
            if compat_info["compatible"]:
                candidate = await self._evaluate_workflow(
                    method,
                    target,
                    parts or [],
                    annotations,
                    anchor_plasmid,
                )
                candidates.append(candidate)
            else:
                # Add incompatible candidate for user visibility
                candidates.append(WorkflowCandidate(
                    method=method,
                    compatible=False,
                    incompatibility_reasons=compat_info["reasons"],
                ))

        # Rank candidates by objective
        ranked = self._rank_candidates(candidates, objective, preferred_method)

        logger.info(f"Found {len([c for c in ranked if c.compatible])} compatible workflows")
        if preferred_method and ranked:
            logger.info(f"Preferred method: {preferred_method.value}, selected: {ranked[0].method.value}")

        return ranked

    def _annotate_cloning_features(
        self,
        target: TargetPlasmid,
        parts: List[ResolvedPart],
    ) -> Dict[str, Any]:
        """
        Extract cloning-relevant features from target and parts.

        Args:
            target: Target plasmid
            parts: Resolved parts

        Returns:
            Dict with extracted features for compatibility checking
        """
        annotations = {
            "target_length": target.length,
            "num_parts": len(parts),
            "part_lengths": [len(p.sequence) for p in parts],
            "topology": target.topology,
            "has_restriction_sites": len(target.restriction_sites) > 0,
            "restriction_sites": target.restriction_sites,
            "has_type_iis_sites": len(target.type_iis_sites) > 0,
            "type_iis_sites": target.type_iis_sites,
            "has_homology": len(target.homology_regions) > 0,
            "homology_regions": target.homology_regions,
            "junctions": target.get_part_junctions(),
        }

        return annotations

    def _check_workflow_compatibility(
        self,
        annotations: Dict[str, Any],
        anchor_plasmid: Optional[TargetPlasmid],
    ) -> Dict[WorkflowMethod, Dict[str, Any]]:
        """
        Fast compatibility check for all workflows.

        Args:
            annotations: Cloning feature annotations
            anchor_plasmid: Anchor plasmid (for SDM)

        Returns:
            Dict mapping workflow method to compatibility info
        """
        compatibility = {}

        # Gibson Assembly
        reasons = []
        num_parts = annotations["num_parts"]
        part_lengths = annotations["part_lengths"]

        if num_parts < self.constraints[WorkflowMethod.GIBSON]["min_fragments"]:
            reasons.append(f"Too few parts ({num_parts}, need at least 2)")
        if num_parts > self.constraints[WorkflowMethod.GIBSON]["max_fragments"]:
            reasons.append(f"Too many parts ({num_parts}, max is 6)")

        if part_lengths:
            min_len = min(part_lengths)
            if min_len < self.constraints[WorkflowMethod.GIBSON]["min_fragment_length"]:
                reasons.append(f"Part too short ({min_len} bp, min 100 bp)")

        compatibility[WorkflowMethod.GIBSON] = {
            "compatible": len(reasons) == 0,
            "reasons": reasons,
        }

        # Golden Gate
        reasons = []
        if num_parts < self.constraints[WorkflowMethod.GOLDEN_GATE]["min_fragments"]:
            reasons.append(f"Too few parts ({num_parts}, need at least 2)")
        if num_parts > self.constraints[WorkflowMethod.GOLDEN_GATE]["max_fragments"]:
            reasons.append(f"Too many parts ({num_parts}, max is 32)")

        # Golden Gate requires Type IIS sites for primer design
        # (Note: In practice, primers will ADD Type IIS sites, so this is informational)

        compatibility[WorkflowMethod.GOLDEN_GATE] = {
            "compatible": len(reasons) == 0,
            "reasons": reasons,
        }

        # Restriction Cloning
        reasons = []
        if num_parts < self.constraints[WorkflowMethod.RESTRICTION]["min_fragments"]:
            reasons.append(f"No parts to insert")
        if num_parts > self.constraints[WorkflowMethod.RESTRICTION]["max_fragments"]:
            reasons.append(f"Too many parts ({num_parts}, restriction typically handles 1-3)")

        if not annotations["has_restriction_sites"]:
            reasons.append("No unique restriction sites found at junctions")

        compatibility[WorkflowMethod.RESTRICTION] = {
            "compatible": len(reasons) == 0,
            "reasons": reasons,
        }

        # Site-Directed Mutagenesis
        reasons = []
        if not anchor_plasmid:
            reasons.append("SDM requires an anchor plasmid")
        if num_parts > 0:
            reasons.append("SDM is for modifications, not multi-part assembly")

        compatibility[WorkflowMethod.SDM] = {
            "compatible": len(reasons) == 0,
            "reasons": reasons,
        }

        # Inventory Gibson (evaluated separately when inventory files present)
        # Will be set to compatible if inventory files are provided
        compatibility[WorkflowMethod.INVENTORY_GIBSON] = {
            "compatible": False,
            "reasons": ["No inventory files provided"],
        }

        # Gateway Cloning
        reasons = []
        if num_parts < 2:
            reasons.append(f"Too few parts ({num_parts}, Gateway needs at least 2)")

        # Gateway requires at least one part / target with native att sites
        # (attP donor, attL/attR entry, or attB destination). Adding att
        # sites via PCR primers is only meaningful when paired with an
        # att-bearing partner. Scan the target sequence for att hits.
        try:
            from ..cloning.gateway_sites import scan_att_sites
            target_atts = scan_att_sites(target.sequence, fuzzy_threshold=0)
            if not target_atts:
                # Also scan each part — a Gateway entry/donor in the inputs
                # is enough.
                any_part_has_att = False
                for p in parts:
                    if scan_att_sites(p.sequence, fuzzy_threshold=0):
                        any_part_has_att = True
                        break
                if not any_part_has_att:
                    reasons.append(
                        "No att sites detected in target or parts — Gateway needs an "
                        "attB/attP/attL/attR-bearing donor or destination"
                    )
        except Exception as _exc:
            # Be conservative: if the scanner blows up, do not silently mark
            # gateway as compatible. Surface as a reason.
            reasons.append(f"att-site scan failed: {_exc}")

        compatibility[WorkflowMethod.GATEWAY] = {
            "compatible": len(reasons) == 0,
            "reasons": reasons,
        }

        return compatibility

    async def _evaluate_workflow(
        self,
        method: WorkflowMethod,
        target: TargetPlasmid,
        parts: List[ResolvedPart],
        annotations: Dict[str, Any],
        anchor_plasmid: Optional[TargetPlasmid],
    ) -> WorkflowCandidate:
        """
        Evaluate a compatible workflow using its operator.

        Args:
            method: Workflow method
            target: Target plasmid
            parts: Resolved parts
            annotations: Cloning annotations
            anchor_plasmid: Anchor plasmid (if applicable)

        Returns:
            WorkflowCandidate with evaluation metrics
        """
        logger.debug(f"Evaluating {method.value}...")

        # TODO: Phase 4 - integrate with actual operators
        # For now, return placeholder evaluations

        candidate = WorkflowCandidate(
            method=method,
            compatible=True,
        )

        # Cost estimates use lab_profile.estimate_workflow_cost: primers
        # length-scaled at $0.24/bp, per-fragment Q5 PCR, the workflow
        # reagent bundle, and one ONT sequencing read. Calendar / risk /
        # confidence remain heuristic priors until the operator runs.
        from ..cloning.lab_profile import estimate_workflow_cost
        n = max(1, len(parts))

        if method == WorkflowMethod.GIBSON:
            candidate.total_cost_usd = estimate_workflow_cost("gibson", n)
            candidate.total_calendar_days = 5.0
            candidate.overall_risk_score = 0.2
            candidate.confidence = 0.9

        elif method == WorkflowMethod.GOLDEN_GATE:
            candidate.total_cost_usd = estimate_workflow_cost(
                "golden_gate", n, avg_primer_len=40,
            )
            candidate.total_calendar_days = 3.0
            candidate.overall_risk_score = 0.15
            candidate.confidence = 0.95

        elif method == WorkflowMethod.RESTRICTION:
            candidate.total_cost_usd = estimate_workflow_cost(
                "restriction", n, primers_per_part=0,  # often no engineered primers
            )
            candidate.total_calendar_days = 4.0
            candidate.overall_risk_score = 0.25
            candidate.confidence = 0.85

        elif method == WorkflowMethod.SDM:
            candidate.total_cost_usd = estimate_workflow_cost(
                "sdm", 1, avg_primer_len=45,
            )
            candidate.total_calendar_days = 6.0
            candidate.overall_risk_score = 0.3
            candidate.confidence = 0.8

        elif method == WorkflowMethod.GATEWAY:
            # Actually call Gateway operator
            try:
                from ..cloning.gateway_operator import GatewayOperator
                
                # Convert parts to modules format
                modules = self._parts_to_modules(parts)
                
                operator = GatewayOperator()
                plan = operator.evaluate(modules, topology=annotations["topology"])
                
                candidate.total_cost_usd = plan.metrics.total_cost_usd
                candidate.total_calendar_days = plan.metrics.total_calendar_days
                candidate.overall_risk_score = plan.metrics.overall_risk_score
                candidate.confidence = 1.0 - plan.metrics.overall_risk_score
                candidate.warnings = plan.warnings
                candidate.build_plan = plan.to_dict()
                
                if not plan.feasible:
                    candidate.compatible = False
                    candidate.incompatibility_reasons = plan.infeasibility_reasons
                
            except Exception as e:
                logger.error(f"Gateway evaluation failed: {e}")
                from ..cloning.lab_profile import estimate_workflow_cost
                candidate.total_cost_usd = estimate_workflow_cost(
                    "gateway", max(1, len(parts)), avg_primer_len=45,
                )
                candidate.total_calendar_days = 4.5
                candidate.overall_risk_score = 0.15
                candidate.confidence = 0.9

        return candidate


    def _has_att_sites(self, sequence: str) -> bool:
        """
        Quick check for any Gateway att sites in sequence.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            True if any att sites found
        """
        try:
            from ..cloning.gateway_sites import scan_att_sites
            sites = scan_att_sites(sequence, fuzzy_threshold=2, search_attB_only=False)
            return len(sites) > 0
        except Exception:
            return False
    
    def _parts_to_modules(self, parts: List[ResolvedPart]) -> List[dict]:
        """
        Convert ResolvedPart list to module dict format for operators.
        
        Args:
            parts: List of resolved parts
            
        Returns:
            List of module dicts
        """
        modules = []
        for i, part in enumerate(parts):
            module = {
                "sequence": part.sequence,
                "canonical_id": part.canonical_id or f"Part_{i+1}",
                "role": getattr(part, "role", "insert"),
                "origin": part.origin,
                "confidence": part.confidence,
            }
            modules.append(module)
        
        return modules

    def _rank_candidates(
        self,
        candidates: List[WorkflowCandidate],
        objective: str,
        preferred_method: Optional[WorkflowMethod] = None,
    ) -> List[WorkflowCandidate]:
        """
        Rank candidates by objective.

        Args:
            candidates: List of workflow candidates
            objective: Optimization objective
            preferred_method: If provided, strongly prefer this method

        Returns:
            Sorted list (best first)
        """
        # If a preferred method is specified and it's compatible, put it first
        if preferred_method:
            preferred = [c for c in candidates if c.method == preferred_method and c.compatible]
            others = [c for c in candidates if c.method != preferred_method or not c.compatible]

            if preferred:
                # Sort the preferred by objective (in case there are multiple evaluations)
                if objective == "balanced":
                    preferred.sort(key=lambda c: c.get_balanced_score())
                elif objective == "cost":
                    preferred.sort(key=lambda c: c.total_cost_usd)
                elif objective == "time":
                    preferred.sort(key=lambda c: c.total_calendar_days)
                elif objective == "risk":
                    preferred.sort(key=lambda c: c.overall_risk_score)

                # Sort others by objective
                if objective == "balanced":
                    others.sort(key=lambda c: c.get_balanced_score())
                elif objective == "cost":
                    others.sort(key=lambda c: (not c.compatible, c.total_cost_usd))
                elif objective == "time":
                    others.sort(key=lambda c: (not c.compatible, c.total_calendar_days))
                elif objective == "risk":
                    others.sort(key=lambda c: (not c.compatible, c.overall_risk_score))

                # Preferred method first, then others
                return preferred + others

        # Standard ranking if no preferred method or preferred not compatible
        if objective == "balanced":
            candidates.sort(key=lambda c: c.get_balanced_score())

        elif objective == "cost":
            # Sort by cost, incompatible last
            candidates.sort(key=lambda c: (
                not c.compatible,
                c.total_cost_usd,
            ))

        elif objective == "time":
            # Sort by time, incompatible last
            candidates.sort(key=lambda c: (
                not c.compatible,
                c.total_calendar_days,
            ))

        elif objective == "risk":
            # Sort by risk, incompatible last
            candidates.sort(key=lambda c: (
                not c.compatible,
                c.overall_risk_score,
            ))

        else:
            logger.warning(f"Unknown objective '{objective}', using 'balanced'")
            candidates.sort(key=lambda c: c.get_balanced_score())

        return candidates
