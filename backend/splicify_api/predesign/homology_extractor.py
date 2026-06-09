"""
Homology-based part extraction from inventory plasmids.

This module provides k-mer based homology search to find segments of inventory
plasmids that match a target sequence. Generalized from inv_gib.py to support
multiple workflows (Gibson, Golden Gate, Restriction, Gateway).
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


@dataclass
class HomologyMatch:
    """
    A segment of an inventory plasmid that matches the target by homology.

    Attributes:
        inventory_name: Name of the inventory plasmid
        inventory_sequence: Full sequence of inventory plasmid
        inventory_start: Start position in inventory (0-indexed)
        inventory_end: End position in inventory (exclusive)
        inventory_orientation: "+" for forward, "-" for reverse complement
        target_start: Start position in target (normalized to [0, L))
        target_end: End position in target (may be > L for circular wrap)
        match_length: Length of exact match (bp)
        sequence: Matched sequence (in target orientation)
    """

    inventory_name: str
    inventory_sequence: str
    inventory_start: int
    inventory_end: int
    inventory_orientation: str  # "+" or "-"
    target_start: int
    target_end: int
    match_length: int
    sequence: str

    def wraps_target(self, target_length: int) -> bool:
        """Check if this match wraps around the circular target."""
        return self.target_end > target_length

    def __repr__(self) -> str:
        return (
            f"HomologyMatch("
            f"inv={self.inventory_name}, "
            f"inv_pos={self.inventory_start}-{self.inventory_end} ({self.inventory_orientation}), "
            f"tgt_pos={self.target_start}-{self.target_end}, "
            f"len={self.match_length} bp)"
        )


class HomologyExtractor:
    """
    Extract parts from inventory plasmids by homology to target sequence.

    Uses k-mer indexing and exact matching to find inventory segments
    that cover the target. This is a generalization of inv_gib logic.
    """

    def __init__(
        self,
        k: int = 20,
        min_match_bp: int = 200,
        max_kmer_pos_per_seed: int = 200,
        max_total_seed_hits: int = 4000,
        min_new_coverage_bp: int = 50,
    ):
        """
        Initialize homology extractor.

        Args:
            k: K-mer seed size (default 20)
            min_match_bp: Minimum contiguous exact match (default 200)
            max_kmer_pos_per_seed: Cap positions per k-mer to avoid repeats (default 200)
            max_total_seed_hits: Max total seed hits per inventory (default 4000)
            min_new_coverage_bp: Accept hit if it adds this much coverage (default 50)
        """
        self.k = k
        self.min_match_bp = min_match_bp
        self.max_kmer_pos_per_seed = max_kmer_pos_per_seed
        self.max_total_seed_hits = max_total_seed_hits
        self.min_new_coverage_bp = min_new_coverage_bp

    def extract_parts(
        self,
        target_sequence: str,
        inventory_sequences: Dict[str, str],
        topology: str = "circular",
    ) -> List[HomologyMatch]:
        """
        Find inventory segments that cover the target by homology.

        Args:
            target_sequence: Target plasmid sequence
            inventory_sequences: Dict mapping inventory name to sequence
            topology: "circular" or "linear" (affects wrap-around handling)

        Returns:
            List of HomologyMatch objects covering the target
        """
        target = target_sequence.upper()
        L = len(target)

        if L < self.k:
            raise ValueError(f"Target too short for k={self.k}: {L} bp")

        logger.info(f"Extracting parts by homology: target {L} bp, {len(inventory_sequences)} inventory plasmids")

        # Build k-mer index of target (with wrap-around for circular)
        if topology == "circular":
            target2 = target + target  # Double for wrap-around
        else:
            target2 = target

        kmer_index = self._build_kmer_index(target2, self.k, self.max_kmer_pos_per_seed)

        # Find best homology match for each inventory
        all_matches: List[HomologyMatch] = []

        for inv_name, inv_seq in inventory_sequences.items():
            matches = self._find_matches_for_inventory(
                target=target,
                target_length=L,
                inventory_name=inv_name,
                inventory_sequence=inv_seq.upper(),
                kmer_index=kmer_index,
            )
            all_matches.extend(matches)

        # Sort by target position
        all_matches.sort(key=lambda m: m.target_start)

        logger.info(f"Found {len(all_matches)} homology matches")
        return all_matches

    def _build_kmer_index(
        self,
        sequence: str,
        k: int,
        max_pos_per_seed: int,
    ) -> Dict[str, List[int]]:
        """
        Build k-mer index: kmer → list of positions.

        Args:
            sequence: Sequence to index
            k: K-mer length
            max_pos_per_seed: Maximum positions to store per k-mer (avoids repeats)

        Returns:
            Dict mapping k-mer to list of start positions
        """
        index: Dict[str, List[int]] = {}
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if kmer not in index:
                index[kmer] = []
            if len(index[kmer]) < max_pos_per_seed:
                index[kmer].append(i)
        return index

    def _find_matches_for_inventory(
        self,
        target: str,
        target_length: int,
        inventory_name: str,
        inventory_sequence: str,
        kmer_index: Dict[str, List[int]],
    ) -> List[HomologyMatch]:
        """
        Find all matches between an inventory plasmid and the target.

        Args:
            target: Target sequence
            target_length: Length of target (for wrap handling)
            inventory_name: Name of inventory plasmid
            inventory_sequence: Sequence of inventory plasmid
            kmer_index: K-mer index of target

        Returns:
            List of HomologyMatch objects
        """
        matches: List[HomologyMatch] = []

        # Try forward orientation
        fwd_match = self._best_exact_match(
            target=target,
            target_length=target_length,
            inventory_name=inventory_name,
            inventory_sequence=inventory_sequence,
            inventory_orientation="+",
            kmer_index=kmer_index,
        )
        if fwd_match:
            matches.append(fwd_match)

        # Try reverse complement orientation
        inv_rc = self._revcomp(inventory_sequence)
        rc_match = self._best_exact_match(
            target=target,
            target_length=target_length,
            inventory_name=inventory_name,
            inventory_sequence=inv_rc,
            inventory_orientation="-",
            kmer_index=kmer_index,
        )
        if rc_match:
            matches.append(rc_match)

        return matches

    def _best_exact_match(
        self,
        target: str,
        target_length: int,
        inventory_name: str,
        inventory_sequence: str,
        inventory_orientation: str,
        kmer_index: Dict[str, List[int]],
    ) -> Optional[HomologyMatch]:
        """
        Find the best exact match for an inventory sequence in a single orientation.

        Uses k-mer seeds to find candidate regions, then extends to maximal exact match.

        Returns:
            Best HomologyMatch if found, None otherwise
        """
        inv_len = len(inventory_sequence)

        # Collect k-mer seed hits
        seed_hits: List[Tuple[int, int]] = []  # (inv_pos, target_pos)
        for i in range(inv_len - self.k + 1):
            kmer = inventory_sequence[i:i + self.k]
            if kmer in kmer_index:
                for tgt_pos in kmer_index[kmer]:
                    seed_hits.append((i, tgt_pos))

                    if len(seed_hits) >= self.max_total_seed_hits:
                        break
            if len(seed_hits) >= self.max_total_seed_hits:
                break

        if not seed_hits:
            return None

        # Extend each seed to maximal exact match
        best_match: Optional[HomologyMatch] = None
        best_length = 0

        for inv_start, tgt_start in seed_hits:
            # Extend rightward
            inv_end = inv_start
            tgt_end = tgt_start
            while (
                inv_end < inv_len
                and tgt_end < len(target) + target_length  # Allow wrap
                and inventory_sequence[inv_end] == target[tgt_end % target_length]
            ):
                inv_end += 1
                tgt_end += 1

            # Extend leftward
            while (
                inv_start > 0
                and tgt_start > 0
                and inventory_sequence[inv_start - 1] == target[(tgt_start - 1) % target_length]
            ):
                inv_start -= 1
                tgt_start -= 1

            match_len = inv_end - inv_start

            if match_len >= self.min_match_bp and match_len > best_length:
                # Normalize target positions
                tgt_start_norm = tgt_start % target_length
                tgt_end_norm = tgt_start_norm + match_len

                matched_seq = inventory_sequence[inv_start:inv_end]

                best_match = HomologyMatch(
                    inventory_name=inventory_name,
                    inventory_sequence=inventory_sequence,
                    inventory_start=inv_start,
                    inventory_end=inv_end,
                    inventory_orientation=inventory_orientation,
                    target_start=tgt_start_norm,
                    target_end=tgt_end_norm,
                    match_length=match_len,
                    sequence=matched_seq,
                )
                best_length = match_len

        return best_match

    def _revcomp(self, seq: str) -> str:
        """Reverse complement a DNA sequence."""
        comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return seq.translate(comp)[::-1]

    def select_non_overlapping_coverage(
        self,
        matches: List[HomologyMatch],
        target_length: int,
    ) -> List[HomologyMatch]:
        """
        Select a non-overlapping set of matches that covers the target.

        Uses greedy algorithm: sort by match length, select longest non-overlapping.

        Args:
            matches: All potential matches
            target_length: Length of target plasmid

        Returns:
            Non-overlapping subset with best coverage
        """
        # Sort by match length (descending)
        sorted_matches = sorted(matches, key=lambda m: m.match_length, reverse=True)

        selected: List[HomologyMatch] = []
        covered = set()  # Set of covered target positions

        for match in sorted_matches:
            # Calculate positions this match would cover
            positions = set()
            for i in range(match.match_length):
                pos = (match.target_start + i) % target_length
                positions.add(pos)

            # Check overlap with already selected matches
            overlap = positions & covered

            # Accept if it adds enough new coverage
            new_coverage = len(positions - covered)
            if new_coverage >= self.min_new_coverage_bp:
                selected.append(match)
                covered |= positions

        # Sort selected matches by target position
        selected.sort(key=lambda m: m.target_start)

        logger.info(f"Selected {len(selected)} non-overlapping matches covering {len(covered)}/{target_length} bp")

        return selected
