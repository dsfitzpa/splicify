"""
pLannotate database integration for feature lookup.

This module provides access to pLannotate's BLAST databases:
- kb_features: Curated biological features (promoters, terminators, etc.)
- kb_cds: Coding sequences
- kb_rna: RNA sequences
- fpbase: Fluorescent proteins (FPbase)
- swissprot: Protein sequences from SwissProt
- Rfam: RNA families

The BLAST databases contain id, name, and sequence in FASTA format.
Knowledge base JSON files may also be referenced for additional metadata.

Features:
- Exact name lookup
- Fuzzy name matching
- Semantic search by description/pattern
- FASTA parsing from BLAST databases
- Caching for performance
"""

import json
import logging
from typing import List, Dict, Any, Optional
from pathlib import Path
from .. import _data
import re
from difflib import SequenceMatcher
from Bio import SeqIO
from io import StringIO

logger = logging.getLogger(__name__)


class PlannotateKnowledgeBase:
    """
    Interface to pLannotate BLAST databases for feature lookup.

    Loads sequences from BLAST database FASTA files and optionally
    enriches with knowledge base metadata.
    """

    def __init__(self, data_dir: Optional[str] = None, blast_db_dir: Optional[str] = None):
        """
        Initialize knowledge base.

        Args:
            data_dir: Path to pLannotate data directory containing knowledge base JSON files
            blast_db_dir: Path to BLAST_dbs directory containing FASTA files
                         Defaults to data_dir/BLAST_dbs if not specified
        """
        self.data_dir = self._find_data_dir(data_dir)
        self.blast_db_dir = Path(blast_db_dir) if blast_db_dir else self.data_dir / "BLAST_dbs"

        # Canonical feature databases (aligned with feature_db_data/ clean-room tree).
        # Each db_name maps to a FASTA in blast_db_dir (discovered via extensions .fasta/.fa/.faa/.fna).
        self.blast_sequences: Dict[str, List[Dict[str, Any]]] = {
            "feature_reference": [],  # GenoLIB nt seed (replaces kb_features + kb_cds + kb_rna)
            "feature_motifs": [],     # short GenoLIB motifs
            "feature_protein": [],    # GenoLIB CDS translations (replaces kb_cds)
            "fpbase": [],              # FPbase (CC-BY-4.0)
            "swissprot": [],           # SwissProt PE<=3 subset
            # Rfam is a CM, not a FASTA; excluded from predesign's FASTA indexer.
        }

        # Optional: Knowledge base metadata (for enrichment)
        self.kb_metadata: Dict[str, Any] = {}

        self._loaded = False

        # Build indexes for fast lookup
        self._sequence_by_id: Dict[str, Dict] = {}
        self._sequence_by_name: Dict[str, List[Dict]] = {}

    def _find_data_dir(self, data_dir: Optional[str]) -> Path:
        """
        Find pLannotate data directory.

        Args:
            data_dir: Explicit data directory path

        Returns:
            Path to data directory

        Search order:
        1. Explicit data_dir argument
        2. Environment variable PLANNOTATE_DATA_DIR
        3. Standard pLannotate installation path
        4. Fallback placeholder
        """
        if data_dir:
            return Path(data_dir)

        # Check environment variable
        import os
        env_dir = os.environ.get("PLANNOTATE_DATA_DIR")
        if env_dir:
            return Path(env_dir)

        # Check standard feature_db_data (clean-room GenoLIB-seeded tree).
        genolib_path = _data.data_path("feature_db_data")
        if genolib_path.exists():
            return genolib_path

        # Legacy pLannotate install fallback (if still present).
        legacy_path = Path.home() / "python-libraries" / "pLannotate" / "plannotate" / "data"
        if legacy_path.exists():
            return legacy_path

        # Final fallback for development/testing.
        return Path("/tmp/plannotate_data")

    def load(self):
        """
        Load sequences from pLannotate BLAST databases.

        BLAST databases are FASTA files containing:
        - kb_features.fasta: Curated features (promoters, terminators, etc.)
        - kb_cds.fasta: Coding sequences
        - kb_rna.fasta: RNA sequences
        - fpbase.fasta: Fluorescent proteins
        - swissprot.fasta: SwissProt proteins
        - Rfam.fasta: RNA families

        Optionally loads knowledge base JSON for metadata enrichment.

        This should be called once at service startup.

        Raises:
            FileNotFoundError: If BLAST database files not found
        """
        if self._loaded:
            logger.debug("Knowledge base already loaded")
            return

        logger.info(f"Loading pLannotate BLAST databases from {self.blast_db_dir}")

        total_sequences = 0

        # Load each BLAST database
        for db_name in self.blast_sequences.keys():
            # Try common BLAST database file extensions
            fasta_paths = [
                self.blast_db_dir / f"{db_name}.fasta",
                self.blast_db_dir / f"{db_name}.fa",
                self.blast_db_dir / f"{db_name}.faa",  # amino acid
                self.blast_db_dir / f"{db_name}.fna",  # nucleic acid
            ]

            loaded = False
            for fasta_path in fasta_paths:
                if fasta_path.exists():
                    try:
                        sequences = self._load_fasta(fasta_path, db_name)
                        self.blast_sequences[db_name] = sequences
                        total_sequences += len(sequences)
                        logger.info(f"  {db_name}: {len(sequences)} sequences from {fasta_path.name}")
                        loaded = True
                        break
                    except Exception as e:
                        logger.warning(f"  Failed to load {fasta_path}: {e}")

            if not loaded:
                logger.warning(f"  {db_name}: No FASTA file found (tried {[p.name for p in fasta_paths]})")

        logger.info(f"Loaded {total_sequences} total sequences from {len(self.blast_sequences)} databases")

        # Optionally load knowledge base metadata for enrichment
        self._load_kb_metadata()

        # Build indexes
        self._build_indexes()
        self._loaded = True

    def _load_fasta(self, fasta_path: Path, db_name: str) -> List[Dict[str, Any]]:
        """
        Load sequences from a FASTA file.

        FASTA header format is typically:
        >id description
        or
        >id|additional|info description

        Args:
            fasta_path: Path to FASTA file
            db_name: Database name (for source tracking)

        Returns:
            List of sequence dicts with id, name, sequence, database
        """
        sequences = []

        try:
            with open(fasta_path) as f:
                for record in SeqIO.parse(f, "fasta"):
                    # Extract ID and description
                    seq_id = record.id
                    description = record.description

                    # Parse name from description (everything after the ID)
                    # Format: ">id name" or ">id|info name"
                    name = description
                    if " " in description:
                        name = description.split(" ", 1)[1].strip()
                    else:
                        name = seq_id

                    sequences.append({
                        "id": seq_id,
                        "name": name,
                        "sequence": str(record.seq).upper(),
                        "database": db_name,
                        "description": description,
                    })

        except Exception as e:
            logger.error(f"Error parsing FASTA file {fasta_path}: {e}")
            raise

        return sequences

    def _load_kb_metadata(self):
        """
        Optionally load knowledge base JSON files for metadata enrichment.

        This provides additional information beyond id/name/sequence that
        may be useful for pattern matching and descriptions.
        """
        # Try to load feature knowledge base
        kb_paths = [
            self.data_dir / "feature_knowledge_base.json",
            self.data_dir / "data" / "feature_knowledge_base.json",
        ]

        for kb_path in kb_paths:
            if kb_path.exists():
                try:
                    with open(kb_path) as f:
                        self.kb_metadata = json.load(f)
                    logger.info(f"Loaded knowledge base metadata from {kb_path.name}")
                    break
                except Exception as e:
                    logger.warning(f"Failed to load KB metadata from {kb_path}: {e}")

        if not self.kb_metadata:
            logger.debug("No knowledge base metadata loaded (optional)")

    def _build_indexes(self):
        """Build lookup indexes for fast search across all databases."""
        logger.debug("Building sequence indexes...")

        # Index all sequences from all BLAST databases
        for db_name, sequences in self.blast_sequences.items():
            for seq in sequences:
                # Index by ID
                seq_id = seq.get("id", "")
                if seq_id:
                    seq_id_lower = seq_id.lower()
                    if seq_id_lower not in self._sequence_by_id:
                        self._sequence_by_id[seq_id_lower] = seq
                    # If duplicate ID, prefer certain databases
                    else:
                        # Priority: fpbase > kb_features > kb_cds > others
                        priority = ["fpbase", "feature_reference", "feature_motifs", "feature_protein", "swissprot"]
                        current_db = self._sequence_by_id[seq_id_lower].get("database", "")
                        if priority.index(db_name) < priority.index(current_db) if current_db in priority else True:
                            self._sequence_by_id[seq_id_lower] = seq

                # Index by name (multiple sequences can have similar names)
                name = seq.get("name", "")
                if name:
                    name_lower = name.lower()
                    if name_lower not in self._sequence_by_name:
                        self._sequence_by_name[name_lower] = []
                    self._sequence_by_name[name_lower].append(seq)

        logger.debug(f"Indexed {len(self._sequence_by_id)} unique sequence IDs across {len(self.blast_sequences)} databases")

    async def search_feature(self, query: str, databases: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """
        Search for sequences by name or ID (exact or fuzzy match).

        Searches across all BLAST databases or specified databases.

        Args:
            query: Feature name or ID to search for (e.g., "CMV", "eGFP", "mCherry")
            databases: Optional list of database names to search
                      Default: searches all databases
                      Options: kb_features, kb_cds, kb_rna, fpbase, swissprot, Rfam

        Returns:
            List of matching sequences, sorted by relevance (best match first)
            Each result contains: id, name, sequence, database, description, confidence

        Examples:
            >>> kb.search_feature("CMV")
            [{"id": "CMV_promoter", "name": "CMV promoter", "sequence": "ATGC...", ...}]

            >>> kb.search_feature("egfp")  # Case insensitive
            [{"id": "eGFP", "name": "enhanced GFP", "sequence": "ATGC...", ...}]

            >>> kb.search_feature("mCherry", databases=["fpbase"])
            [{"id": "mCherry", "database": "fpbase", ...}]
        """
        if not self._loaded:
            self.load()

        query_lower = query.lower().strip()
        results = []

        # 1. Exact ID match (highest priority)
        if query_lower in self._sequence_by_id:
            seq = self._sequence_by_id[query_lower]
            # Filter by database if specified
            if databases is None or seq.get("database") in databases:
                results.append(self._format_result(seq, confidence=1.0))

        # 2. Exact name match
        if query_lower in self._sequence_by_name:
            for seq in self._sequence_by_name[query_lower]:
                if databases is None or seq.get("database") in databases:
                    # Avoid duplicates
                    if seq not in [r.get("_seq") for r in results]:
                        results.append(self._format_result(seq, confidence=1.0))

        # 3. Fuzzy matching if no exact matches
        if not results:
            results = self._fuzzy_search_sequences(query_lower, databases)

        # Remove internal _seq key used for deduplication
        for r in results:
            r.pop("_seq", None)

        return results

    async def search_pattern(
        self,
        pattern: str,
        databases: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Search for sequences matching a description or pattern.

        Searches in name and description fields across all sequences.

        Args:
            pattern: Description or pattern to match
                    (e.g., "strong mammalian promoter", "fluorescent protein", "GFP")
            databases: Optional list of database names to search

        Returns:
            List of matching sequences, sorted by relevance
            Each result contains: id, name, sequence, database, description, confidence

        Examples:
            >>> kb.search_pattern("strong promoter")
            [{"id": "CMV", ...}, {"id": "EF1a", ...}]

            >>> kb.search_pattern("fluorescent protein", databases=["fpbase"])
            [{"id": "eGFP", ...}, {"id": "mCherry", ...}]

            >>> kb.search_pattern("RFP")  # Red fluorescent proteins
            [{"id": "mCherry", ...}, {"id": "mRuby", ...}]
        """
        if not self._loaded:
            self.load()

        pattern_lower = pattern.lower().strip()
        results = []

        # Get all sequences from specified databases
        sequences_to_search = []
        for db_name, sequences in self.blast_sequences.items():
            if databases is None or db_name in databases:
                sequences_to_search.extend(sequences)

        # Search in name and description fields
        for seq in sequences_to_search:
            score = self._pattern_match_score(pattern_lower, seq)
            if score > 0.3:  # Threshold for relevance
                result = self._format_result(seq, confidence=score)
                results.append(result)

        # Sort by confidence (highest first)
        results.sort(key=lambda x: x["confidence"], reverse=True)

        # Remove internal _seq key
        for r in results:
            r.pop("_seq", None)

        return results

    async def search_protein(self, query: str) -> List[Dict[str, Any]]:
        """
        Search SwissProt protein database specifically.

        This is a convenience method for searching only the swissprot database.

        Args:
            query: Protein name or ID (e.g., "Cas9", "GFP")

        Returns:
            List of matching proteins with sequences

        Examples:
            >>> kb.search_protein("Cas9")
            [{"id": "CAS9_STRP1", "name": "Cas9", "sequence": "ATGC...", ...}]

            >>> kb.search_protein("GFP")
            [{"id": "GFP_AEQVI", "name": "Green fluorescent protein", ...}]
        """
        return await self.search_feature(query, databases=["swissprot"])

    def _fuzzy_search_sequences(
        self,
        query: str,
        databases: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Perform fuzzy string matching on sequence names and IDs.

        Args:
            query: Search query (lowercase)
            databases: Optional list of database names to search

        Returns:
            List of matching sequences with confidence scores
        """
        results = []

        # Get sequences to search
        sequences_to_search = []
        for db_name, sequences in self.blast_sequences.items():
            if databases is None or db_name in databases:
                sequences_to_search.extend(sequences)

        for seq in sequences_to_search:
            name = seq.get("name", "").lower()
            seq_id = seq.get("id", "").lower()

            # Calculate similarity scores
            name_score = SequenceMatcher(None, query, name).ratio()
            id_score = SequenceMatcher(None, query, seq_id).ratio()
            max_score = max(name_score, id_score)

            # Check if query is substring (boost score)
            if query in name or query in seq_id:
                max_score = max(max_score, 0.8)

            if max_score > 0.6:  # Threshold for fuzzy match
                result = self._format_result(seq, confidence=max_score)
                results.append(result)

        # Sort by confidence
        results.sort(key=lambda x: x["confidence"], reverse=True)

        return results

    def _pattern_match_score(self, pattern: str, seq: Dict[str, Any]) -> float:
        """
        Calculate how well a sequence matches a description pattern.

        Args:
            pattern: Search pattern (lowercase)
            seq: Sequence dict from BLAST database

        Returns:
            Match score (0.0-1.0)
        """
        name = seq.get("name", "").lower()
        description = seq.get("description", "").lower()
        seq_id = seq.get("id", "").lower()

        # Split pattern into keywords
        keywords = re.findall(r'\w+', pattern)

        # Count keyword matches
        matches = 0
        for keyword in keywords:
            if keyword in name:
                matches += 2  # Name matches worth more
            elif keyword in description:
                matches += 1.5
            elif keyword in seq_id:
                matches += 1

        # Normalize score
        if not keywords:
            return 0.0

        score = matches / (len(keywords) * 2)  # Max score of 1.0
        return min(score, 1.0)

    def _format_result(
        self,
        seq: Dict[str, Any],
        confidence: float
    ) -> Dict[str, Any]:
        """
        Format a sequence into a standard result dict.

        Args:
            seq: Sequence dict from BLAST database
            confidence: Match confidence score

        Returns:
            Formatted result dict for part resolution
        """
        return {
            "id": seq.get("id", "unknown"),
            "name": seq.get("name", seq.get("id", "")),
            "sequence": seq.get("sequence", ""),
            "type": self._infer_type(seq),
            "description": seq.get("description", seq.get("name", "")),
            "confidence": confidence,
            "database": seq.get("database", "unknown"),
            "_seq": seq,  # Keep reference for deduplication
        }

    def _infer_type(self, seq: Dict[str, Any]) -> str:
        """
        Infer sequence type from database and name.

        Args:
            seq: Sequence dict

        Returns:
            Inferred type (promoter, CDS, terminator, etc.)
        """
        db = seq.get("database", "")
        name = seq.get("name", "").lower()
        seq_id = seq.get("id", "").lower()

        # Database-based inference
        if db == "fpbase":
            return "fluorescent_protein"
        elif db == "kb_cds":
            return "CDS"
        elif db == "kb_rna" or db == "Rfam":
            return "RNA"
        elif db == "swissprot":
            return "protein"

        # Name-based inference
        if "promoter" in name or "promoter" in seq_id:
            return "promoter"
        elif "terminator" in name or "polya" in name:
            return "terminator"
        elif "gfp" in name or "rfp" in name or "yfp" in name or "cfp" in name:
            return "fluorescent_protein"
        elif "gene" in name or "cds" in name:
            return "CDS"
        elif "rna" in name or "rrna" in name or "trna" in name:
            return "RNA"
        elif "origin" in name or "ori" in name:
            return "origin"
        elif "resistance" in name or "marker" in name:
            return "marker"

        return "feature"


    def get_databases(self) -> Dict[str, int]:
        """
        Get statistics about loaded databases.

        Returns:
            Dict mapping database name to number of sequences loaded

        Example:
            >>> kb.get_databases()
            {
                "kb_features": 1523,
                "kb_cds": 2456,
                "kb_rna": 456,
                "fpbase": 234,
                "swissprot": 12456,
                "Rfam": 789
            }
        """
        if not self._loaded:
            self.load()

        return {
            db_name: len(sequences)
            for db_name, sequences in self.blast_sequences.items()
        }

    def get_total_sequences(self) -> int:
        """Get total number of sequences across all databases."""
        if not self._loaded:
            self.load()

        return sum(len(seqs) for seqs in self.blast_sequences.values())


# Singleton instance for global access
_kb_instance: Optional[PlannotateKnowledgeBase] = None


def get_knowledge_base(
    data_dir: Optional[str] = None,
    blast_db_dir: Optional[str] = None
) -> PlannotateKnowledgeBase:
    """
    Get or create the global knowledge base instance.

    Args:
        data_dir: Path to pLannotate data directory (only used on first call)
        blast_db_dir: Path to BLAST_dbs directory (only used on first call)

    Returns:
        PlannotateKnowledgeBase instance

    Example:
        >>> from splicify_api.predesign import get_knowledge_base
        >>> kb = get_knowledge_base()
        >>> kb.load()
        >>> results = await kb.search_feature("CMV")
    """
    global _kb_instance

    if _kb_instance is None:
        _kb_instance = PlannotateKnowledgeBase(data_dir, blast_db_dir)

    return _kb_instance
