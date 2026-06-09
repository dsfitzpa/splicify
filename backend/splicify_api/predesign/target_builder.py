"""
Target plasmid builder for cloning workflows.

This module constructs explicit target plasmid representations from resolved parts,
including sequence assembly, feature annotation, and cloning site scanning.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
import re

from .part_resolver import ResolvedPart

logger = logging.getLogger(__name__)


@dataclass
class TargetPlasmid:
    """
    Explicit representation of a target plasmid.

    This is the complete target that will be built by the chosen cloning workflow.

    Attributes:
        sequence: Complete target DNA sequence (uppercase)
        length: Sequence length in bp
        topology: "circular" or "linear"
        parts: Ordered list of resolved parts that make up the target
        features: GenBank features (if target from uploaded file)
        restriction_sites: Map of enzyme names to cut positions
        type_iis_sites: Map of Type IIS enzyme names to positions
        homology_regions: Detected homology regions between parts
        metadata: Additional metadata (GC content, complexity, etc.)
    """

    sequence: str
    length: int
    topology: str
    parts: List[ResolvedPart] = field(default_factory=list)
    features: List[Dict[str, Any]] = field(default_factory=list)
    restriction_sites: Dict[str, List[int]] = field(default_factory=dict)
    type_iis_sites: Dict[str, List[int]] = field(default_factory=dict)
    homology_regions: List[Dict[str, Any]] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate and normalize target plasmid."""
        # Normalize sequence
        self.sequence = self.sequence.upper()
        self.length = len(self.sequence)

        # Calculate GC content
        if self.sequence:
            gc_count = self.sequence.count('G') + self.sequence.count('C')
            self.metadata['gc_content'] = gc_count / len(self.sequence) if len(self.sequence) > 0 else 0.0

    def get_part_junctions(self) -> List[Dict[str, Any]]:
        """
        Get junction information between parts.

        Returns:
            List of junction dicts with start, end, part1, part2 info
        """
        junctions = []
        if len(self.parts) < 2:
            return junctions

        current_pos = 0
        for i in range(len(self.parts) - 1):
            part1 = self.parts[i]
            part2 = self.parts[i + 1]

            junction_start = current_pos + len(part1.sequence)
            junction_end = junction_start  # Point junction

            junctions.append({
                "junction_id": i + 1,
                "part1_name": part1.name,
                "part2_name": part2.name,
                "position": junction_start,
                "part1_end": junction_start,
                "part2_start": junction_start,
            })

            current_pos += len(part1.sequence)

        # For circular plasmids, add junction between last and first part
        if self.topology == "circular" and len(self.parts) > 1:
            last_part = self.parts[-1]
            first_part = self.parts[0]

            junctions.append({
                "junction_id": len(self.parts),
                "part1_name": last_part.name,
                "part2_name": first_part.name,
                "position": self.length,
                "part1_end": self.length,
                "part2_start": 0,
                "wraps_origin": True,
            })

        return junctions


class TargetPlasmidBuilder:
    """
    Build target plasmid representations from resolved parts.

    This class handles:
    - Sequence assembly from ordered parts
    - Restriction site scanning
    - Type IIS site detection
    - Homology region identification
    - Feature preservation from uploaded files
    """

    def __init__(self):
        """Initialize target builder."""
        # Type IIS enzymes and their recognition sites
        self.type_iis_enzymes = {
            "BsaI": "GGTCTC",
            "BsmBI": "CGTCTC",
            "BbsI": "GAAGAC",
            "SapI": "GCTCTTC",
        }

    def build_from_parts(
        self,
        parts: List[ResolvedPart],
        assembly_order: str = "listed",
        topology: str = "circular",
    ) -> TargetPlasmid:
        """
        Build target plasmid from ordered resolved parts.

        Args:
            parts: List of resolved parts in assembly order
            assembly_order: "listed" or "homology_based"
            topology: "circular" or "linear"

        Returns:
            TargetPlasmid with assembled sequence and annotations

        Raises:
            ValueError: If parts list is empty or assembly fails
        """
        if not parts:
            raise ValueError("Cannot build target from empty parts list")

        logger.info(f"Building target from {len(parts)} parts (order={assembly_order}, topology={topology})")

        # Order parts based on strategy
        if assembly_order == "homology_based":
            ordered_parts = self._order_by_homology(parts)
        else:  # "listed"
            ordered_parts = parts

        # Assemble sequence by concatenation
        sequence = "".join(part.sequence for part in ordered_parts)

        # Create target
        target = TargetPlasmid(
            sequence=sequence,
            length=len(sequence),
            topology=topology,
            parts=ordered_parts,
        )

        # Annotate cloning features
        self._annotate_restriction_sites(target)
        self._annotate_type_iis_sites(target)
        self._detect_homology_regions(target)

        logger.info(f"Built target plasmid: {target.length} bp, {len(target.parts)} parts")

        return target

    def build_from_upload(
        self,
        sequence: str,
        features: Optional[List[Dict]] = None,
        topology: str = "circular",
        metadata: Optional[Dict[str, Any]] = None,
    ) -> TargetPlasmid:
        """
        Build target plasmid from uploaded file.

        Args:
            sequence: DNA sequence
            features: GenBank features from file
            topology: "circular" or "linear"
            metadata: Additional metadata

        Returns:
            TargetPlasmid with preserved features
        """
        logger.info(f"Building target from uploaded file: {len(sequence)} bp")

        target = TargetPlasmid(
            sequence=sequence,
            length=len(sequence),
            topology=topology,
            features=features or [],
            metadata=metadata or {},
        )

        # Annotate cloning features
        self._annotate_restriction_sites(target)
        self._annotate_type_iis_sites(target)

        return target

    def _order_by_homology(self, parts: List[ResolvedPart]) -> List[ResolvedPart]:
        """
        Order parts by detecting homology regions.

        This will be implemented in Phase 6 using inv_gib k-mer logic.

        Args:
            parts: Unordered parts

        Returns:
            Parts ordered by homology

        Raises:
            NotImplementedError: Phase 6 feature
        """
        raise NotImplementedError(
            "Homology-based part ordering will be implemented in Phase 6"
        )

    def _annotate_restriction_sites(self, target: TargetPlasmid):
        """
        Scan and annotate common restriction enzyme sites.

        This reuses logic from junction.py._scan_restriction_sites()

        Args:
            target: Target plasmid to annotate (modified in place)
        """
        # Common restriction enzymes for cloning
        enzymes = {
            "EcoRI": "GAATTC",
            "BamHI": "GGATCC",
            "HindIII": "AAGCTT",
            "XhoI": "CTCGAG",
            "SalI": "GTCGAC",
            "PstI": "CTGCAG",
            "SmaI": "CCCGGG",
            "KpnI": "GGTACC",
            "SacI": "GAGCTC",
            "XbaI": "TCTAGA",
            "NheI": "GCTAGC",
            "SpeI": "ACTAGT",
            "NotI": "GCGGCCGC",
            "ApaI": "GGGCCC",
        }

        for enzyme, site in enzymes.items():
            positions = self._find_all_positions(target.sequence, site)
            if positions:
                target.restriction_sites[enzyme] = positions

        logger.debug(f"Found {len(target.restriction_sites)} RE sites")

    def _annotate_type_iis_sites(self, target: TargetPlasmid):
        """
        Scan and annotate Type IIS restriction enzyme sites.

        Type IIS enzymes (BsaI, BsmBI, BbsI) are used in Golden Gate assembly.

        Args:
            target: Target plasmid to annotate (modified in place)
        """
        for enzyme, site in self.type_iis_enzymes.items():
            # Search both forward and reverse complement
            positions_fwd = self._find_all_positions(target.sequence, site)
            positions_rev = self._find_all_positions(target.sequence, self._reverse_complement(site))

            positions = []
            for pos in positions_fwd:
                positions.append(pos)
            for pos in positions_rev:
                positions.append(pos)

            if positions:
                target.type_iis_sites[enzyme] = sorted(positions)

        logger.debug(f"Found {len(target.type_iis_sites)} Type IIS sites")

    def _detect_homology_regions(self, target: TargetPlasmid):
        """
        Detect homology regions between adjacent parts.

        Args:
            target: Target plasmid to annotate (modified in place)
        """
        if len(target.parts) < 2:
            return

        # Check for overlaps between adjacent parts
        # This is a simple implementation - Phase 6 will use sophisticated k-mer alignment

        for i in range(len(target.parts) - 1):
            part1 = target.parts[i]
            part2 = target.parts[i + 1]

            # Check for 10-50 bp overlaps at part boundaries
            for overlap_len in range(50, 9, -1):
                if len(part1.sequence) < overlap_len or len(part2.sequence) < overlap_len:
                    continue

                part1_end = part1.sequence[-overlap_len:]
                part2_start = part2.sequence[:overlap_len]

                if part1_end == part2_start:
                    target.homology_regions.append({
                        "part1": part1.name,
                        "part2": part2.name,
                        "length": overlap_len,
                        "sequence": part1_end,
                        "type": "junction_overlap",
                    })
                    break

        logger.debug(f"Found {len(target.homology_regions)} homology regions")

    def _find_all_positions(self, sequence: str, pattern: str) -> List[int]:
        """
        Find all positions of a pattern in sequence.

        Args:
            sequence: DNA sequence to search
            pattern: Pattern to find

        Returns:
            List of 0-indexed positions
        """
        positions = []
        pos = 0
        while True:
            pos = sequence.find(pattern, pos)
            if pos == -1:
                break
            positions.append(pos)
            pos += 1
        return positions

    def _reverse_complement(self, sequence: str) -> str:
        """
        Get reverse complement of DNA sequence.

        Args:
            sequence: DNA sequence

        Returns:
            Reverse complement
        """
        complement = str.maketrans("ATGCRYSWKMBDHVN", "TACGYRSWMKVHDBN")
        return sequence.translate(complement)[::-1]
