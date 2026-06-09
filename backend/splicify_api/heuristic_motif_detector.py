"""
Heuristic Motif Detector
========================
Detect biological motifs (ORFs, Kozak sequences, RBS, start/stop codons) 
using pattern matching rather than KB-based classification.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
import re


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


@dataclass
class MotifHit:
    """Detected motif with heuristic evidence"""
    motif_type: str  # 'orf', 'kozak', 'rbs', 'start_codon', 'stop_codon', etc.
    name: str
    start: int
    end: int
    strand: int
    sequence: str
    confidence: float
    rules_fired: List[str]
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "motif_type": self.motif_type,
            "name": self.name,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "sequence": self.sequence,
            "confidence": self.confidence,
            "rules_fired": self.rules_fired,
            "metadata": self.metadata
        }


class HeuristicMotifDetector:
    """Detect motifs using heuristic rules rather than KB classification"""

    # Kozak patterns (IUPAC notation converted to regex)
    KOZAK_STRONG = re.compile(r'GCC[AG]CCATGG', re.IGNORECASE)
    KOZAK_ADEQUATE = re.compile(r'[AG]CCATGG', re.IGNORECASE)
    KOZAK_WEAK = re.compile(r'[ACGT]{3}ATGG', re.IGNORECASE)
    
    # RBS/Shine-Dalgarno patterns
    RBS_STRONG = re.compile(r'(AGGAGG|GGAGGA|AGGAG[GA]|[AG]GGAGG)', re.IGNORECASE)
    RBS_ADEQUATE = re.compile(r'(GGAGG|AGGAG|GAGGA|AGGA[GT])', re.IGNORECASE)
    RBS_WEAK = re.compile(r'AGG', re.IGNORECASE)
    
    # Stop codons
    STOP_CODONS = {'TAA', 'TAG', 'TGA'}

    def __init__(self, sequence: str, circular: bool = True):
        self.sequence = sequence.upper()
        self.seq_len = len(sequence)
        self.circular = circular

    def detect_orfs(self, min_length_aa: int = 100) -> List[MotifHit]:
        """
        Detect ORFs >= min_length_aa amino acids on both strands.
        
        Args:
            min_length_aa: Minimum ORF length in amino acids (default 100)
            
        Returns:
            List of MotifHit objects for detected ORFs
        """
        hits = []
        min_length_bp = min_length_aa * 3

        # Forward strand - three reading frames
        for frame in range(3):
            orfs = self._find_orfs_in_frame(self.sequence, frame)
            for orf_start, orf_end, orf_seq in orfs:
                orf_length = orf_end - orf_start
                if orf_length >= min_length_bp:
                    hits.append(MotifHit(
                        motif_type='orf',
                        name=f'ORF ({orf_length // 3} aa)',
                        start=orf_start,
                        end=orf_end,
                        strand=1,
                        sequence=orf_seq[:60] + '...' if len(orf_seq) > 60 else orf_seq,
                        confidence=0.85,
                        rules_fired=['CDS-01'],
                        metadata={
                            'frame': frame,
                            'length_aa': orf_length // 3,
                            'length_bp': orf_length,
                            'start_codon': orf_seq[:3],
                            'stop_codon': orf_seq[-3:] if len(orf_seq) >= 3 else ''
                        }
                    ))

        # Reverse strand - three reading frames
        rc_seq = reverse_complement(self.sequence)
        for frame in range(3):
            orfs = self._find_orfs_in_frame(rc_seq, frame)
            for orf_start, orf_end, orf_seq in orfs:
                orf_length = orf_end - orf_start
                if orf_length >= min_length_bp:
                    # Convert coordinates to forward strand reference
                    genomic_start = self.seq_len - orf_end
                    genomic_end = self.seq_len - orf_start
                    hits.append(MotifHit(
                        motif_type='orf',
                        name=f'ORF ({orf_length // 3} aa)',
                        start=genomic_start,
                        end=genomic_end,
                        strand=-1,
                        sequence=orf_seq[:60] + '...' if len(orf_seq) > 60 else orf_seq,
                        confidence=0.85,
                        rules_fired=['CDS-01'],
                        metadata={
                            'frame': frame,
                            'length_aa': orf_length // 3,
                            'length_bp': orf_length,
                            'start_codon': orf_seq[:3],
                            'stop_codon': orf_seq[-3:] if len(orf_seq) >= 3 else ''
                        }
                    ))

        # Sort by start position
        hits.sort(key=lambda h: (h.start, -h.metadata.get('length_bp', 0)))
        return hits

    def _find_orfs_in_frame(self, seq: str, frame: int) -> List[Tuple[int, int, str]]:
        """Find all ATG...stop ORFs in a specific reading frame."""
        orfs = []
        i = frame
        seq_len = len(seq)
        
        while i <= seq_len - 3:
            codon = seq[i:i+3]
            if codon == 'ATG':
                # Found start codon, search for stop
                orf_start = i
                j = i + 3
                while j <= seq_len - 3:
                    stop_codon = seq[j:j+3]
                    if stop_codon in self.STOP_CODONS:
                        orf_end = j + 3  # Include stop codon
                        orf_seq = seq[orf_start:orf_end]
                        orfs.append((orf_start, orf_end, orf_seq))
                        break
                    j += 3
            i += 3
        
        return orfs

    def detect_kozak_sequences(self) -> List[MotifHit]:
        """
        Detect Kozak consensus sequences on forward strand.
        
        Patterns:
        - Strong: GCC[AG]CCATGG (rules CDS-04)
        - Adequate: [AG]CCATGG (rules CDS-05)
        - Weak: [ACGT]{3}ATGG (rules CDS-06)
        
        Returns:
            List of MotifHit objects with kozak strength in metadata
        """
        hits = []
        
        # Strong Kozak (highest priority)
        for match in self.KOZAK_STRONG.finditer(self.sequence):
            hits.append(MotifHit(
                motif_type='kozak',
                name='Kozak (strong)',
                start=match.start(),
                end=match.end(),
                strand=1,
                sequence=match.group(),
                confidence=0.95,
                rules_fired=['CDS-04'],
                metadata={
                    'strength': 'strong',
                    'atg_position': match.start() + 6  # Position of ATG within pattern
                }
            ))
        
        # Adequate Kozak
        for match in self.KOZAK_ADEQUATE.finditer(self.sequence):
            # Check if this overlaps with a strong Kozak
            overlaps = any(
                h.motif_type == 'kozak' and 
                h.metadata.get('strength') == 'strong' and
                not (match.end() <= h.start or match.start() >= h.end)
                for h in hits
            )
            if not overlaps:
                hits.append(MotifHit(
                    motif_type='kozak',
                    name='Kozak (adequate)',
                    start=match.start(),
                    end=match.end(),
                    strand=1,
                    sequence=match.group(),
                    confidence=0.75,
                    rules_fired=['CDS-05'],
                    metadata={
                        'strength': 'adequate',
                        'atg_position': match.start() + 4
                    }
                ))
        
        # Also detect on reverse strand
        rc_seq = reverse_complement(self.sequence)
        for match in self.KOZAK_STRONG.finditer(rc_seq):
            genomic_end = self.seq_len - match.start()
            genomic_start = self.seq_len - match.end()
            hits.append(MotifHit(
                motif_type='kozak',
                name='Kozak (strong)',
                start=genomic_start,
                end=genomic_end,
                strand=-1,
                sequence=match.group(),
                confidence=0.95,
                rules_fired=['CDS-04'],
                metadata={
                    'strength': 'strong',
                    'atg_position': genomic_start + match.end() - match.start() - 6 - 3
                }
            ))
        
        hits.sort(key=lambda h: h.start)
        return hits

    def detect_shine_dalgarno(self) -> List[MotifHit]:
        """
        Detect Shine-Dalgarno / RBS sequences.
        
        Patterns (upstream of ATG, 4-20bp distance):
        - Strong: AGGAGG, GGAGGA, AGGAG[GA], [AG]GGAGG (CDS-07)
        - Adequate: GGAGG, AGGAG, GAGGA, AGGA[GT] (CDS-08)
        - Weak: AGG (CDS-09)
        
        Returns:
            List of MotifHit objects
        """
        hits = []
        
        # Find all ATG positions first
        atg_positions = []
        for match in re.finditer(r'ATG', self.sequence, re.IGNORECASE):
            atg_positions.append(match.start())
        
        # Strong RBS
        for match in self.RBS_STRONG.finditer(self.sequence):
            rbs_end = match.end()
            # Check if there's an ATG 4-20bp downstream
            for atg_pos in atg_positions:
                distance = atg_pos - rbs_end
                if 4 <= distance <= 20:
                    hits.append(MotifHit(
                        motif_type='rbs',
                        name='RBS (strong)',
                        start=match.start(),
                        end=match.end(),
                        strand=1,
                        sequence=match.group(),
                        confidence=0.95,
                        rules_fired=['CDS-07'],
                        metadata={
                            'strength': 'strong',
                            'distance_to_atg': distance,
                            'atg_position': atg_pos
                        }
                    ))
                    break
        
        # Adequate RBS (avoid duplicates with strong)
        strong_positions = {(h.start, h.end) for h in hits}
        for match in self.RBS_ADEQUATE.finditer(self.sequence):
            if (match.start(), match.end()) in strong_positions:
                continue
            # Overlaps check
            overlaps = any(
                not (match.end() <= h.start or match.start() >= h.end)
                for h in hits if h.metadata.get('strength') == 'strong'
            )
            if overlaps:
                continue
                
            rbs_end = match.end()
            for atg_pos in atg_positions:
                distance = atg_pos - rbs_end
                if 4 <= distance <= 20:
                    hits.append(MotifHit(
                        motif_type='rbs',
                        name='RBS (adequate)',
                        start=match.start(),
                        end=match.end(),
                        strand=1,
                        sequence=match.group(),
                        confidence=0.75,
                        rules_fired=['CDS-08'],
                        metadata={
                            'strength': 'adequate',
                            'distance_to_atg': distance,
                            'atg_position': atg_pos
                        }
                    ))
                    break
        
        hits.sort(key=lambda h: h.start)
        return hits

    def detect_start_codons_in_context(
        self,
        kozak_hits: List[MotifHit],
        rbs_hits: List[MotifHit]
    ) -> List[MotifHit]:
        """
        Detect ATG start codons only in biologically relevant context.
        
        Context requirements:
        - ATG within Kozak sequence -> ATG (Kozak)
        - ATG 4-20bp downstream of RBS -> ATG (RBS)
        
        Args:
            kozak_hits: Previously detected Kozak sequences
            rbs_hits: Previously detected RBS sequences
            
        Returns:
            List of MotifHit for context-validated start codons
        """
        hits = []
        seen_positions = set()
        
        # ATGs from Kozak context
        for kozak in kozak_hits:
            atg_pos = kozak.metadata.get('atg_position', kozak.start + 6)
            if atg_pos >= 0 and atg_pos + 3 <= self.seq_len:
                if atg_pos not in seen_positions:
                    seen_positions.add(atg_pos)
                    hits.append(MotifHit(
                        motif_type='start_codon',
                        name=f'ATG (Kozak-{kozak.metadata.get("strength", "unknown")})',
                        start=atg_pos,
                        end=atg_pos + 3,
                        strand=kozak.strand,
                        sequence='ATG',
                        confidence=kozak.confidence,
                        rules_fired=['CDS-02'],
                        metadata={
                            'context': 'kozak',
                            'kozak_strength': kozak.metadata.get('strength'),
                            'kozak_start': kozak.start
                        }
                    ))
        
        # ATGs from RBS context
        for rbs in rbs_hits:
            atg_pos = rbs.metadata.get('atg_position')
            if atg_pos and atg_pos not in seen_positions:
                seen_positions.add(atg_pos)
                hits.append(MotifHit(
                    motif_type='start_codon',
                    name=f'ATG (RBS-{rbs.metadata.get("strength", "unknown")})',
                    start=atg_pos,
                    end=atg_pos + 3,
                    strand=rbs.strand,
                    sequence='ATG',
                    confidence=rbs.confidence,
                    rules_fired=['CDS-02'],
                    metadata={
                        'context': 'rbs',
                        'rbs_strength': rbs.metadata.get('strength'),
                        'rbs_start': rbs.start,
                        'distance_from_rbs': rbs.metadata.get('distance_to_atg')
                    }
                ))
        
        hits.sort(key=lambda h: h.start)
        return hits

    def detect_stop_codons_in_frame(
        self,
        orf_hits: List[MotifHit]
    ) -> List[MotifHit]:
        """
        Detect stop codons (TAA, TAG, TGA) at ORF boundaries.
        
        Args:
            orf_hits: Previously detected ORFs
            
        Returns:
            List of MotifHit at ORF end positions
        """
        hits = []
        seen_positions = set()
        
        for orf in orf_hits:
            # Stop codon is at the end of the ORF
            stop_start = orf.end - 3
            stop_end = orf.end
            
            if (stop_start, orf.strand) in seen_positions:
                continue
            seen_positions.add((stop_start, orf.strand))
            
            stop_seq = orf.metadata.get('stop_codon', '')
            if not stop_seq and stop_start >= 0 and stop_end <= self.seq_len:
                if orf.strand == 1:
                    stop_seq = self.sequence[stop_start:stop_end]
                else:
                    # For reverse strand, get from reverse complement
                    rc_pos_start = self.seq_len - stop_end
                    rc_pos_end = self.seq_len - stop_start
                    rc_seq = reverse_complement(self.sequence)
                    if 0 <= rc_pos_start < rc_pos_end <= len(rc_seq):
                        stop_seq = rc_seq[rc_pos_start:rc_pos_end]
            
            # Determine stop codon name
            stop_name = {
                'TAA': 'TAA (ochre)',
                'TAG': 'TAG (amber)',
                'TGA': 'TGA (opal)'
            }.get(stop_seq.upper(), f'Stop ({stop_seq})')
            
            hits.append(MotifHit(
                motif_type='stop_codon',
                name=stop_name,
                start=stop_start,
                end=stop_end,
                strand=orf.strand,
                sequence=stop_seq,
                confidence=0.90,
                rules_fired=['CDS-03'],
                metadata={
                    'orf_start': orf.start,
                    'orf_length_aa': orf.metadata.get('length_aa')
                }
            ))
        
        hits.sort(key=lambda h: h.start)
        return hits

    def detect_all_motifs(self, features: List[Dict] = None) -> Dict[str, List[MotifHit]]:
        """
        Run all heuristic motif detection.
        
        Args:
            features: Optional pLannotate features for context
            
        Returns:
            Dict with keys: 'orfs', 'kozak', 'rbs', 'start_codons', 'stop_codons'
        """
        # Detect ORFs first (slowest operation)
        orfs = self.detect_orfs(min_length_aa=100)
        
        # Detect translation initiation context
        kozak = self.detect_kozak_sequences()
        rbs = self.detect_shine_dalgarno()
        
        # Context-validated start codons
        start_codons = self.detect_start_codons_in_context(kozak, rbs)
        
        # Stop codons at ORF ends
        stop_codons = self.detect_stop_codons_in_frame(orfs)
        
        return {
            'orfs': orfs,
            'kozak': kozak,
            'rbs': rbs,
            'start_codons': start_codons,
            'stop_codons': stop_codons
        }


# Convenience function for direct use
def detect_motifs(sequence: str, circular: bool = True) -> Dict[str, List[MotifHit]]:
    """
    Detect all biological motifs in a sequence using heuristic rules.
    
    Args:
        sequence: DNA sequence
        circular: Whether the sequence is circular (default True)
        
    Returns:
        Dictionary of detected motifs by type
    """
    detector = HeuristicMotifDetector(sequence, circular)
    return detector.detect_all_motifs()
