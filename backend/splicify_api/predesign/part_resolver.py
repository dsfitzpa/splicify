"""
Multi-source part resolution for cloning workflows.

This module resolves PartSpecification objects to concrete DNA sequences using
multiple strategies:
1. Direct sequence validation
2. Uploaded file parsing (GenBank, FASTA)
3. Feature name lookup (pLannotate knowledge bases)
4. Pattern matching (fuzzy search)
5. Homology-based detection (k-mer alignment)
"""

import logging
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from Bio import SeqIO
from io import StringIO
import re

from .design_request import PartSpecification, InputSource

logger = logging.getLogger(__name__)


@dataclass
class ResolvedPart:
    """
    A part with resolved sequence and metadata.

    Attributes:
        name: Part name
        sequence: Resolved DNA sequence (uppercase)
        length: Sequence length in bp
        source: Original input source type
        canonical_id: Canonical identifier (e.g., "eGFP", "CMV_promoter")
        description: Human-readable description
        role: Biological role (promoter, CDS, terminator, etc.)
        origin: Where the sequence came from (knowledge_base, user_file, etc.)
        source_detail: Additional source information (file name, KB entry, etc.)
        features: GenBank features if parsed from file
        confidence: Resolution confidence score (0.0-1.0)
        warnings: Any warnings during resolution
    """

    name: str
    sequence: str
    length: int
    source: InputSource
    canonical_id: Optional[str] = None
    description: Optional[str] = None
    role: Optional[str] = None
    origin: str = "user_provided"
    source_detail: Optional[str] = None
    features: List[Dict[str, Any]] = field(default_factory=list)
    confidence: float = 1.0
    warnings: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Validate resolved part."""
        # Normalize sequence to uppercase
        self.sequence = self.sequence.upper()

        # Validate DNA alphabet
        valid_bases = set("ATGCRYSWKMBDHVN")
        invalid_bases = set(self.sequence) - valid_bases
        if invalid_bases:
            self.warnings.append(
                f"Sequence contains non-DNA characters: {sorted(invalid_bases)}"
            )

        # Update length
        self.length = len(self.sequence)

        # Set canonical_id to name if not provided
        if not self.canonical_id:
            self.canonical_id = self.name

    def to_module_dict(self) -> Dict[str, Any]:
        """
        Convert to module dict format for operators.

        This format is compatible with existing Gibson/Golden Gate operators.
        """
        return {
            "canonical_id": self.canonical_id,
            "description": self.description or self.name,
            "role": self.role or "unknown",
            "sequence": self.sequence,
            "length": self.length,
            "origin": self.origin,
            "source": self.source_detail or "user_provided",
        }


class PartResolver:
    """
    Resolves part specifications to concrete sequences.

    Resolution strategies:
    1. Direct sequence → validate and wrap
    2. Uploaded file → parse with Bio.SeqIO
    3. Feature name → search knowledge bases
    4. Feature pattern → fuzzy matching
    5. Homology → k-mer alignment (delegated to inv_gib)
    """

    def __init__(self, knowledge_base=None):
        """
        Initialize part resolver.

        Args:
            knowledge_base: Optional PlannotateKnowledgeBase instance
                          (will be lazy-loaded if not provided)
        """
        self.knowledge_base = knowledge_base
        self._kb_loaded = False

    async def resolve_all(
        self,
        parts: List[PartSpecification],
        context: Optional[Dict[str, Any]] = None,
    ) -> List[ResolvedPart]:
        """
        Resolve all part specifications to concrete sequences.

        Args:
            parts: List of part specifications to resolve
            context: Optional context dict containing:
                - file_cache: Dict[file_id, file_content]
                - session_id: Session identifier
                - inventory_files: List of inventory file contents

        Returns:
            List of resolved parts in the same order as input

        Raises:
            ValueError: If a part cannot be resolved
        """
        context = context or {}
        resolved = []

        for i, part in enumerate(parts):
            try:
                logger.info(f"Resolving part {i+1}/{len(parts)}: {part.name} (source: {part.source.value})")
                resolved_part = await self._resolve_single(part, context)
                resolved.append(resolved_part)
            except Exception as e:
                logger.error(f"Failed to resolve part '{part.name}': {e}")
                raise ValueError(f"Failed to resolve part '{part.name}': {e}") from e

        return resolved

    async def _resolve_single(
        self,
        part: PartSpecification,
        context: Dict[str, Any],
    ) -> ResolvedPart:
        """Resolve a single part specification."""
        if part.source == InputSource.DIRECT_SEQUENCE:
            return self._resolve_direct(part)

        elif part.source == InputSource.UPLOADED_FILE:
            return await self._resolve_file(part, context)

        elif part.source == InputSource.FEATURE_NAME:
            return await self._resolve_feature_name(part)

        elif part.source == InputSource.FEATURE_PATTERN:
            return await self._resolve_pattern(part)

        elif part.source == InputSource.HOMOLOGY_DERIVED:
            return await self._resolve_homology(part, context)

        else:
            raise ValueError(f"Unknown input source: {part.source}")

    def _resolve_direct(self, part: PartSpecification) -> ResolvedPart:
        """
        Resolve a direct sequence specification.

        Args:
            part: Part with source=DIRECT_SEQUENCE

        Returns:
            ResolvedPart with validated sequence

        Raises:
            ValueError: If sequence is invalid
        """
        if not part.sequence:
            raise ValueError(f"Part '{part.name}': DIRECT_SEQUENCE requires sequence")

        # Clean sequence (remove whitespace, newlines)
        sequence = re.sub(r'\s+', '', part.sequence)

        # Validate minimum length
        if len(sequence) < 1:
            raise ValueError(f"Part '{part.name}': sequence is empty")

        return ResolvedPart(
            name=part.name,
            sequence=sequence,
            length=len(sequence),
            source=part.source,
            canonical_id=part.name,
            description=f"User-provided sequence ({len(sequence)} bp)",
            role=part.role or "unknown",
            origin="user_input",
            source_detail="direct_sequence",
            confidence=1.0,
        )

    async def _resolve_file(
        self,
        part: PartSpecification,
        context: Dict[str, Any],
    ) -> ResolvedPart:
        """
        Resolve a part from an uploaded file.

        Args:
            part: Part with source=UPLOADED_FILE
            context: Context dict with file_cache

        Returns:
            ResolvedPart with sequence from file

        Raises:
            ValueError: If file cannot be parsed
        """
        if not part.file_id:
            raise ValueError(f"Part '{part.name}': UPLOADED_FILE requires file_id")

        # Get file content from cache
        file_cache = context.get("file_cache", {})
        if part.file_id not in file_cache:
            raise ValueError(f"Part '{part.name}': file_id '{part.file_id}' not found in cache")

        file_content = file_cache[part.file_id]
        file_name = file_content.get("name", "unknown.gb")

        # Detect file format from extension or content
        file_format = self._detect_format(file_name, file_content.get("content", ""))

        try:
            # Parse with Bio.SeqIO
            content_str = file_content.get("content", "")
            handle = StringIO(content_str)
            records = list(SeqIO.parse(handle, file_format))

            if not records:
                raise ValueError(f"No sequences found in file '{file_name}'")

            # Use first record
            record = records[0]
            sequence = str(record.seq)

            # Extract features
            features = []
            for feature in record.features:
                features.append({
                    "type": feature.type,
                    "location": str(feature.location),
                    "qualifiers": dict(feature.qualifiers),
                })

            # Extract description
            description = record.description or f"Sequence from {file_name}"

            return ResolvedPart(
                name=part.name or record.name or file_name,
                sequence=sequence,
                length=len(sequence),
                source=part.source,
                canonical_id=record.name or part.name,
                description=description,
                role=part.role or "plasmid",
                origin="uploaded_file",
                source_detail=file_name,
                features=features,
                confidence=1.0,
            )

        except Exception as e:
            raise ValueError(f"Failed to parse file '{file_name}': {e}") from e

    def _detect_format(self, filename: str, content: str) -> str:
        """
        Detect file format from filename or content.

        Args:
            filename: File name
            content: File content

        Returns:
            Format string for Bio.SeqIO ("genbank", "fasta", etc.)
        """
        # Check extension
        if filename.endswith((".gb", ".gbk", ".genbank")):
            return "genbank"
        elif filename.endswith((".fa", ".fasta", ".fna")):
            return "fasta"
        elif filename.endswith(".dna"):
            return "snapgene"

        # Check content
        if content.strip().startswith("LOCUS"):
            return "genbank"
        elif content.strip().startswith(">"):
            return "fasta"

        # Default to genbank
        return "genbank"

    async def _resolve_feature_name(self, part: PartSpecification) -> ResolvedPart:
        """
        Resolve a part by feature name lookup.

        Args:
            part: Part with source=FEATURE_NAME

        Returns:
            ResolvedPart with sequence from knowledge base

        Raises:
            ValueError: If feature not found
        """
        if not part.feature_name:
            raise ValueError(f"Part '{part.name}': FEATURE_NAME requires feature_name")

        # Lazy-load knowledge base
        if not self._kb_loaded:
            await self._load_knowledge_base()

        # Search knowledge base
        if self.knowledge_base is None:
            raise ValueError("Knowledge base not available for feature lookup")

        results = await self.knowledge_base.search_feature(part.feature_name)

        if not results:
            raise ValueError(
                f"Feature '{part.feature_name}' not found in knowledge base. "
                f"Try using a pattern search or provide the sequence directly."
            )

        # Use best match (first result)
        best_match = results[0]

        return ResolvedPart(
            name=part.name,
            sequence=best_match["sequence"],
            length=len(best_match["sequence"]),
            source=part.source,
            canonical_id=best_match.get("id", part.feature_name),
            description=best_match.get("description", part.feature_name),
            role=part.role or best_match.get("type", "unknown"),
            origin="knowledge_base",
            source_detail=f"pLannotate: {best_match.get('database', 'unknown')}",
            confidence=best_match.get("confidence", 1.0),
        )

    async def _resolve_pattern(self, part: PartSpecification) -> ResolvedPart:
        """
        Resolve a part by pattern/description search.

        Args:
            part: Part with source=FEATURE_PATTERN

        Returns:
            ResolvedPart with sequence from best match

        Raises:
            ValueError: If no matches found
        """
        if not part.feature_pattern:
            raise ValueError(f"Part '{part.name}': FEATURE_PATTERN requires feature_pattern")

        # Lazy-load knowledge base
        if not self._kb_loaded:
            await self._load_knowledge_base()

        if self.knowledge_base is None:
            raise ValueError("Knowledge base not available for pattern search")

        results = await self.knowledge_base.search_pattern(part.feature_pattern)

        if not results:
            raise ValueError(
                f"No features found matching pattern '{part.feature_pattern}'. "
                f"Try refining your search or provide the sequence directly."
            )

        # Use best match (first result)
        best_match = results[0]

        # Add warning if confidence is low
        warnings = []
        if best_match.get("confidence", 1.0) < 0.7:
            warnings.append(
                f"Low confidence match ({best_match['confidence']:.2f}) for pattern '{part.feature_pattern}'"
            )

        return ResolvedPart(
            name=part.name,
            sequence=best_match["sequence"],
            length=len(best_match["sequence"]),
            source=part.source,
            canonical_id=best_match.get("id", part.name),
            description=best_match.get("description", part.feature_pattern),
            role=part.role or best_match.get("type", "unknown"),
            origin="knowledge_base",
            source_detail=f"pLannotate: {best_match.get('database', 'unknown')}",
            confidence=best_match.get("confidence", 1.0),
            warnings=warnings,
        )

    async def _resolve_homology(
        self,
        part: PartSpecification,
        context: Dict[str, Any],
    ) -> ResolvedPart:
        """
        Resolve a part using homology-based detection.

        This delegates to the inv_gib k-mer alignment logic.

        Args:
            part: Part with source=HOMOLOGY_DERIVED
            context: Context dict with inventory_files

        Returns:
            ResolvedPart with sequence detected by homology

        Raises:
            ValueError: If homology detection fails
        """
        # TODO: Phase 6 - integrate with inv_gib k-mer logic
        raise NotImplementedError(
            "Homology-based part resolution will be implemented in Phase 6"
        )

    async def _load_knowledge_base(self):
        """Lazy-load the knowledge base."""
        if self._kb_loaded:
            return

        # Import here to avoid circular dependencies
        try:
            from .knowledge_base import get_knowledge_base

            # Get singleton instance and load if needed
            self.knowledge_base = get_knowledge_base()

            # Load knowledge bases if not already loaded
            try:
                self.knowledge_base.load()
                logger.info("Knowledge base loaded successfully")
            except Exception as e:
                logger.warning(f"Failed to load knowledge base: {e}")
                logger.warning("Feature lookup will be disabled")
                self.knowledge_base = None

            self._kb_loaded = True

        except ImportError as e:
            logger.warning(f"Knowledge base not available: {e}")
            self.knowledge_base = None
            self._kb_loaded = True
