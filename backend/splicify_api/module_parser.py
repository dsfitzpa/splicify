"""
Module Parser - Motif and Boundary Detection

Provides:
- Motif detection (start/stop codons, 2A peptides, Kozak sequences)
- Module boundary detection
- Sequence analysis utilities

This module is part of the annotation pipeline restructuring to create
cleaner separation between motif/boundary detection, validation/linting,
and inter-module analysis.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
import re
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


# =============================================================================
# CONSTANTS
# =============================================================================

# Standard codon table
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Only ATG for start codon detection (not alternative start codons)
START_CODON = 'ATG'
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# 2A peptide signatures (conserved C-terminal motif)
PEPTIDE_2A_SIGNATURES = {
    'T2A': 'EGRGSLLTCGDVEENPGP',
    'P2A': 'GSGATNFSLLKQAGDVEENPGP',
    'E2A': 'QCTNYALLKLAGDVESNPGP',
    'F2A': 'VKQTLNFDLLKLAGDVESNPGP',
}

# Consensus 2A C-terminal motif
PEPTIDE_2A_MOTIF = re.compile(r'[A-Z]{5,}GDVE[ES]NPGP', re.IGNORECASE)

# Kozak patterns
KOZAK_STRONG_PATTERN = re.compile(r'GCC[AG]CCATGG', re.IGNORECASE)
KOZAK_ADEQUATE_PATTERN = re.compile(r'[AG]CCATGG', re.IGNORECASE)
KOZAK_WEAK_PATTERN = re.compile(r'[ACGT]{3}ATGG', re.IGNORECASE)


# Shine-Dalgarno / RBS patterns with tiered strength
# Strong: Full consensus or near-consensus (5-6bp)
RBS_STRONG_PATTERN = re.compile(r"AGGAGG|GGAGGA|AGGAG[GA]|[AG]GGAGG", re.IGNORECASE)
# Adequate: Partial consensus (4-5bp)
RBS_ADEQUATE_PATTERN = re.compile(r"GGAGG|AGGAG|GAGGA|AGGA[GT]", re.IGNORECASE)
# Weak: Minimal core (3bp, less specific)
RBS_WEAK_PATTERN = re.compile(r"AGG|GAG", re.IGNORECASE)
# Legacy combined pattern (for backwards compatibility)
RBS_PATTERN = re.compile(r"AGGAGG|GGAGGA|AGGAG[GA]|[AG]GGAGG|GGAGG", re.IGNORECASE)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class MotifHit:
    """Represents a detected motif in the sequence"""
    motif_type: str  # 'start_codon', 'stop_codon', '2a_peptide', 'kozak', 'internal_stop'
    name: str  # e.g., 'ATG', 'TAA', 'P2A', 'Kozak (strong)'
    start: int  # 0-indexed start position
    end: int  # 0-indexed end position (exclusive)
    strand: int  # 1 = forward, -1 = reverse
    sequence: str  # The actual sequence matched
    description: str = ""  # Human-readable description
    confidence: float = 1.0  # Confidence score (0-1)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'motif_type': self.motif_type,
            'name': self.name,
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'sequence': self.sequence,
            'description': self.description,
            'confidence': self.confidence,
            'metadata': self.metadata,
        }


@dataclass
class ModuleBoundary:
    """Represents a detected module boundary"""
    module_type: str  # e.g., 'pol2_expression', 'pol3_expression', 'transcript'
    boundary_type: str  # 'start' or 'end'
    position: int  # 0-indexed position
    strand: int
    feature_id: Optional[str] = None  # Associated feature ID
    feature_name: Optional[str] = None
    confidence: float = 1.0


@dataclass
class TranscriptBoundary:
    """Represents transcript start/end boundaries"""
    start_position: int  # Transcript start (promoter end + Kozak/ATG)
    end_position: int  # Transcript end (stop codon + polyA signal)
    strand: int
    promoter_id: Optional[str] = None
    terminator_id: Optional[str] = None


@dataclass
class ProteinBoundary:
    """Represents protein-level CDS boundaries (removing introns/2A)"""
    start_position: int
    end_position: int
    strand: int
    cds_feature_ids: List[str] = field(default_factory=list)
    has_2a_split: bool = False
    intron_positions: List[Tuple[int, int]] = field(default_factory=list)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def translate_sequence(seq: str, frame: int = 0) -> str:
    """Translate DNA sequence to amino acid sequence"""
    seq = seq.upper().replace('U', 'T')
    protein = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, 'X')
            protein.append(aa)
    return ''.join(protein)


def find_start_codon(seq: str) -> Optional[int]:
    """
    Find position of first ATG start codon (0-indexed).
    Only detects ATG (not alternative start codons CTG/GTG/TTG).
    """
    seq = seq.upper()
    pos = seq.find(START_CODON)
    return pos if pos >= 0 else None


def find_stop_codons(seq: str, frame: int = 0) -> List[int]:
    """Find positions of all stop codons in frame (0-indexed)"""
    seq = seq.upper()
    stops = []
    for i in range(frame, len(seq) - 2, 3):
        if seq[i:i+3] in STOP_CODONS:
            stops.append(i)
    return stops


def detect_2a_peptides(protein_seq: str) -> List[Tuple[str, int, int]]:
    """
    Detect 2A peptide sequences in protein.
    
    Returns list of (name, start, end) tuples where positions are
    amino acid indices in the protein sequence.
    """
    detected = []
    protein_upper = protein_seq.upper()

    for name, signature in PEPTIDE_2A_SIGNATURES.items():
        sig_upper = signature.upper()
        start = 0
        while True:
            pos = protein_upper.find(sig_upper, start)
            if pos < 0:
                break
            detected.append((name, pos, pos + len(signature)))
            start = pos + 1

    # Also check generic motif if no specific 2A found
    if not detected:
        for match in PEPTIDE_2A_MOTIF.finditer(protein_upper):
            detected.append(("unknown_2A", match.start(), match.end()))

    return detected


def get_sequence_slice(sequence: str, start: int, end: int, strand: int, circular: bool = True) -> str:
    """
    Extract sequence for a feature, handling circularity and strand.
    Coordinates are 0-based, end is exclusive.
    """
    if not sequence:
        return ""

    seq_len = len(sequence)

    if circular and start > end:
        # Feature wraps around origin
        subseq = sequence[start:] + sequence[:end]
    else:
        start = max(0, min(start, seq_len))
        end = max(0, min(end, seq_len))
        subseq = sequence[start:end]

    # Reverse complement if on minus strand
    if strand == -1:
        subseq = reverse_complement(subseq)

    return subseq


# =============================================================================
# MOTIF DETECTOR CLASS
# =============================================================================

class MotifDetector:
    """Detect short motifs for precise module boundaries"""

    # Maximum distance to search for start codon after promoter
    MAX_PROMOTER_TO_START_BP = 500

    def __init__(self, sequence: str, circular: bool = True):
        self.sequence = sequence.upper()
        self.seq_len = len(sequence)
        self.circular = circular

    def detect_start_codons(
        self,
        features: List[Dict[str, Any]],
        search_after_promoters: bool = True
    ) -> List[MotifHit]:
        """
        Find ATG start codons within reasonable distance after promoter end.
        
        Args:
            features: List of feature dicts with 'start', 'end', 'strand', 'role' keys
            search_after_promoters: If True, only search downstream of promoters
            
        Returns:
            List of MotifHit objects for detected start codons
        """
        hits = []
        
        if search_after_promoters:
            # Find promoters
            promoters = [
                f for f in features 
                if 'promoter' in f.get('role', '').lower() or 
                   'promoter' in f.get('feature_type', '').lower()
            ]
            
            for promoter in promoters:
                strand = promoter.get('strand', 1)
                
                if strand == 1:
                    # Forward strand: search after promoter end
                    search_start = promoter.get('end', 0)
                    search_end = min(search_start + self.MAX_PROMOTER_TO_START_BP, self.seq_len)
                    region = self.sequence[search_start:search_end]
                else:
                    # Reverse strand: search before promoter start
                    search_end = promoter.get('start', 0)
                    search_start = max(0, search_end - self.MAX_PROMOTER_TO_START_BP)
                    region = reverse_complement(self.sequence[search_start:search_end])
                
                # Find all ATG in region
                pos = 0
                while True:
                    atg_pos = region.find(START_CODON, pos)
                    if atg_pos < 0:
                        break
                    
                    # Calculate genomic position
                    if strand == 1:
                        genomic_pos = search_start + atg_pos
                    else:
                        genomic_pos = search_end - atg_pos - 3
                    
                    hits.append(MotifHit(
                        motif_type='start_codon',
                        name='ATG',
                        start=genomic_pos,
                        end=genomic_pos + 3,
                        strand=strand,
                        sequence='ATG',
                        description=f'Start codon {self.MAX_PROMOTER_TO_START_BP}bp downstream of {promoter.get("feature_name", "promoter")}',
                        metadata={
                            'promoter_id': promoter.get('instance_id'),
                            'distance_from_promoter': atg_pos,
                        }
                    ))
                    pos = atg_pos + 1
        else:
            # Search entire sequence
            pos = 0
            while True:
                atg_pos = self.sequence.find(START_CODON, pos)
                if atg_pos < 0:
                    break
                hits.append(MotifHit(
                    motif_type='start_codon',
                    name='ATG',
                    start=atg_pos,
                    end=atg_pos + 3,
                    strand=1,
                    sequence='ATG',
                    description='Start codon (ATG)',
                ))
                pos = atg_pos + 1
        
        return hits

    def detect_start_loss(
        self,
        promoters: List[Dict[str, Any]],
        cds_features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """
        Detect missing start codon within reasonable distance after promoter.
        
        Flags if no ATG found within expected range after promoter end.
        """
        hits = []
        
        for promoter in promoters:
            strand = promoter.get('strand', 1)
            promoter_end = promoter.get('end', 0) if strand == 1 else promoter.get('start', 0)
            
            # Check if there's a CDS that starts within the expected range
            expected_cds = None
            for cds in cds_features:
                if cds.get('strand', 1) != strand:
                    continue
                cds_start = cds.get('start', 0) if strand == 1 else cds.get('end', 0)
                distance = abs(cds_start - promoter_end)
                
                if distance < self.MAX_PROMOTER_TO_START_BP:
                    expected_cds = cds
                    break
            
            if not expected_cds:
                continue
            
            # Search for ATG in the region
            if strand == 1:
                search_start = promoter_end
                search_end = min(promoter_end + self.MAX_PROMOTER_TO_START_BP, self.seq_len)
                region = self.sequence[search_start:search_end]
            else:
                search_end = promoter.get('start', 0)
                search_start = max(0, search_end - self.MAX_PROMOTER_TO_START_BP)
                region = reverse_complement(self.sequence[search_start:search_end])
            
            if START_CODON not in region:
                hits.append(MotifHit(
                    motif_type='start_loss',
                    name='Missing ATG',
                    start=promoter_end if strand == 1 else search_start,
                    end=search_end if strand == 1 else promoter.get('start', 0),
                    strand=strand,
                    sequence='',
                    description=f'No start codon (ATG) found within {self.MAX_PROMOTER_TO_START_BP}bp after promoter',
                    confidence=0.8,
                    metadata={
                        'promoter_id': promoter.get('instance_id'),
                        'expected_cds_id': expected_cds.get('instance_id'),
                    }
                ))
        
        return hits

    def detect_stop_codons(
        self,
        cds_features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """Find TAA/TAG/TGA stop codons at CDS ends"""
        hits = []
        
        for cds in cds_features:
            strand = cds.get('strand', 1)
            start = cds.get('start', 0)
            end = cds.get('end', 0)
            
            # Get CDS sequence
            cds_seq = get_sequence_slice(self.sequence, start, end, strand, self.circular)
            
            if len(cds_seq) < 3:
                continue
            
            # Check for stop codon at end
            last_codon = cds_seq[-3:].upper()
            if last_codon in STOP_CODONS:
                # Calculate genomic position of stop codon
                if strand == 1:
                    stop_start = end - 3
                    stop_end = end
                else:
                    stop_start = start
                    stop_end = start + 3
                
                hits.append(MotifHit(
                    motif_type='stop_codon',
                    name=last_codon,
                    start=stop_start,
                    end=stop_end,
                    strand=strand,
                    sequence=last_codon,
                    description=f'Stop codon at end of {cds.get("feature_name", "CDS")}',
                    metadata={'cds_id': cds.get('instance_id')}
                ))
        
        return hits

    def detect_2a_peptides_in_cds(
        self,
        cds_features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """Identify P2A, T2A, E2A, F2A peptide sequences in CDS features"""
        hits = []
        
        for cds in cds_features:
            strand = cds.get('strand', 1)
            start = cds.get('start', 0)
            end = cds.get('end', 0)
            
            # Get CDS sequence and translate
            cds_seq = get_sequence_slice(self.sequence, start, end, strand, self.circular)
            protein_seq = translate_sequence(cds_seq)
            
            # Detect 2A peptides
            peptides = detect_2a_peptides(protein_seq)
            
            for name, aa_start, aa_end in peptides:
                # Convert amino acid positions to nucleotide positions
                nt_start = aa_start * 3
                nt_end = aa_end * 3
                
                # Convert to genomic coordinates
                if strand == 1:
                    genomic_start = start + nt_start
                    genomic_end = start + nt_end
                else:
                    genomic_start = end - nt_end
                    genomic_end = end - nt_start
                
                hits.append(MotifHit(
                    motif_type='2a_peptide',
                    name=name,
                    start=genomic_start,
                    end=genomic_end,
                    strand=strand,
                    sequence=protein_seq[aa_start:aa_end],
                    description=f'{name} self-cleaving peptide',
                    metadata={
                        'cds_id': cds.get('instance_id'),
                        'aa_start': aa_start,
                        'aa_end': aa_end,
                    }
                ))
        
        return hits

    def detect_internal_stops(
        self,
        cds_features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """Find premature stop codons within CDS"""
        hits = []
        
        for cds in cds_features:
            strand = cds.get('strand', 1)
            start = cds.get('start', 0)
            end = cds.get('end', 0)
            
            # Get CDS sequence
            cds_seq = get_sequence_slice(self.sequence, start, end, strand, self.circular)
            
            # Find start codon to determine frame
            start_pos = find_start_codon(cds_seq)
            frame = start_pos % 3 if start_pos is not None else 0
            
            # Find all stop codons in frame
            stops = find_stop_codons(cds_seq, frame)
            
            # Exclude terminal stop (last 3bp)
            internal_stops = [s for s in stops if s < len(cds_seq) - 3]
            
            for stop_pos in internal_stops:
                codon = cds_seq[stop_pos:stop_pos+3]
                
                # Convert to genomic coordinates
                if strand == 1:
                    genomic_start = start + stop_pos
                    genomic_end = genomic_start + 3
                else:
                    genomic_end = end - stop_pos
                    genomic_start = genomic_end - 3
                
                hits.append(MotifHit(
                    motif_type='internal_stop',
                    name=codon,
                    start=genomic_start,
                    end=genomic_end,
                    strand=strand,
                    sequence=codon,
                    description=f'Premature stop codon in {cds.get("feature_name", "CDS")}',
                    confidence=1.0,
                    metadata={
                        'cds_id': cds.get('instance_id'),
                        'position_in_cds': stop_pos,
                    }
                ))
        
        return hits

    def detect_introns_from_kb(
        self,
        features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """
        Detect introns using KB feature_class annotations only.
        Does not use splice site motif detection.
        """
        hits = []
        
        for feature in features:
            # Check if KB marks this as an intron
            feature_class = feature.get('kb_feature_class', '')
            if feature_class == 'intron':
                hits.append(MotifHit(
                    motif_type='intron',
                    name=feature.get('feature_name', 'intron'),
                    start=feature.get('start', 0),
                    end=feature.get('end', 0),
                    strand=feature.get('strand', 1),
                    sequence='',  # Don't include full intron sequence
                    description='Intron (from KB annotation)',
                    metadata={
                        'feature_id': feature.get('instance_id'),
                        'source': 'kb_feature_class',
                    }
                ))
        
        return hits

    def detect_kozak_sequences(
        self,
        cds_features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """Find Kozak consensus sequences at CDS starts"""
        hits = []
        
        for cds in cds_features:
            strand = cds.get('strand', 1)
            start = cds.get('start', 0)
            end = cds.get('end', 0)
            
            # Get context around CDS start (need ~10bp upstream)
            if strand == 1:
                context_start = max(0, start - 10)
                context_end = min(start + 7, self.seq_len)
                context_seq = self.sequence[context_start:context_end]
                cds_start_in_context = start - context_start
            else:
                context_start = max(0, end - 7)
                context_end = min(end + 10, self.seq_len)
                context_seq = reverse_complement(self.sequence[context_start:context_end])
                cds_start_in_context = context_end - end
            
            # Look for ATG
            atg_pos = context_seq.find('ATG')
            if atg_pos < 0:
                continue
            
            # Check Kozak patterns (need sequence around ATG)
            kozak_context = context_seq[max(0, atg_pos-6):atg_pos+4]
            
            strength = 'unknown'
            kozak_seq = ''
            
            if len(kozak_context) >= 7:
                if KOZAK_STRONG_PATTERN.search(kozak_context):
                    strength = 'strong'
                    match = KOZAK_STRONG_PATTERN.search(kozak_context)
                    kozak_seq = match.group() if match else kozak_context
                elif KOZAK_ADEQUATE_PATTERN.search(kozak_context):
                    strength = 'adequate'
                    match = KOZAK_ADEQUATE_PATTERN.search(kozak_context)
                    kozak_seq = match.group() if match else kozak_context
                elif KOZAK_WEAK_PATTERN.search(kozak_context):
                    strength = 'weak'
                    match = KOZAK_WEAK_PATTERN.search(kozak_context)
                    kozak_seq = match.group() if match else kozak_context
            
            if strength != 'unknown':
                # Calculate genomic position
                if strand == 1:
                    kozak_start = context_start + atg_pos - 3
                    kozak_end = context_start + atg_pos + 4
                else:
                    kozak_end = context_end - atg_pos + 3
                    kozak_start = context_end - atg_pos - 4
                
                hits.append(MotifHit(
                    motif_type='kozak',
                    name=f'Kozak ({strength})',
                    start=max(0, kozak_start),
                    end=min(self.seq_len, kozak_end),
                    strand=strand,
                    sequence=kozak_seq,
                    description=f'{strength.title()} Kozak consensus at {cds.get("feature_name", "CDS")} start',
                    confidence=1.0 if strength == 'strong' else 0.8 if strength == 'adequate' else 0.6,
                    metadata={
                        'cds_id': cds.get('instance_id'),
                        'strength': strength,
                    }
                ))
        
        return hits

    def _is_bacterial_promoter_name(self, name: str) -> bool:
        """Check if promoter name indicates bacterial/phage promoter."""
        name_lower = (name or '').lower()
        bacterial_indicators = ['t7', 't5', 't3', 'sp6', 'lac', 'tac', 'trc',
                               'arabad', 'pbad', 'rha', 'tet', 'ampr', 'bla']
        return any(ind in name_lower for ind in bacterial_indicators)

    def detect_rbs_sequences(
        self,
        features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """
        Detect RBS/Shine-Dalgarno sequences:
        1. Features with kb_feature_class = 'RBS' or 'ribosome_binding_site'
        2. Pattern search in 40bp window around bacterial/bacteriophage promoter ends
           (20bp before to 20bp after promoter end)
        3. Pattern search upstream of bacterial marker CDS (existing logic)

        Uses tiered strength detection:
        - Strong: AGGAGG or near-consensus (5-6bp) - high confidence
        - Adequate: GGAGG, AGGAG (4-5bp) - moderate confidence
        - Weak: AGG or GAG (3bp) - low confidence, may indicate partial RBS
        """
        hits = []

        # 1. Detect RBS from KB feature_class
        for f in features:
            kb_class = f.get('kb_feature_class', '').lower()
            if kb_class in ('rbs', 'ribosome_binding_site', 'shine_dalgarno'):
                hits.append(MotifHit(
                    motif_type='rbs',
                    name='RBS (annotated)',
                    start=f.get('start', 0),
                    end=f.get('end', 0),
                    strand=f.get('strand', 1),
                    sequence=self.sequence[f.get('start', 0):f.get('end', 0)],
                    description='Ribosome binding site from annotation',
                    confidence=1.0,
                    metadata={'source': 'kb_feature_class'}
                ))

        # 2. Search 40bp window around bacterial/bacteriophage promoter ends
        bacterial_promoters = [
            f for f in features
            if (f.get('role', '') == 'promoter'  # Generic promoter (not pol2/pol3)
                and f.get('role', '') not in ('pol2_promoter', 'pol3_promoter'))
            or f.get('kb_polymerase_class', '') == 'bacterial_or_phage'
            or self._is_bacterial_promoter_name(f.get('feature_name', ''))
        ]

        for prom in bacterial_promoters:
            strand = prom.get('strand', 1)
            prom_end = prom.get('end', 0) if strand == 1 else prom.get('start', 0)

            if strand == 1:
                search_start = max(0, prom_end - 20)
                search_end = min(prom_end + 50, self.seq_len)
                region = self.sequence[search_start:search_end].upper()
            else:
                search_start = max(0, prom_end - 50)
                search_end = min(prom_end + 20, self.seq_len)
                region = reverse_complement(self.sequence[search_start:search_end])

            # Search for RBS patterns (strongest first)
            for pattern, name, strength, confidence in [
                (RBS_STRONG_PATTERN, 'RBS (strong)', 'strong', 0.95),
                (RBS_ADEQUATE_PATTERN, 'RBS (adequate)', 'adequate', 0.75),
            ]:
                match = pattern.search(region)
                if match:
                    if strand == 1:
                        genomic_start = search_start + match.start()
                        genomic_end = search_start + match.end()
                    else:
                        genomic_end = search_end - match.start()
                        genomic_start = search_end - match.end()

                    distance_from_prom = genomic_start - prom_end if strand == 1 else prom_end - genomic_end

                    # Avoid duplicates
                    if not any(h.start == genomic_start and h.end == genomic_end for h in hits):
                        hits.append(MotifHit(
                            motif_type='rbs',
                            name=name,
                            start=genomic_start,
                            end=genomic_end,
                            strand=strand,
                            sequence=match.group(),
                            description=f'{strength.title()} RBS near {prom.get("feature_name", "promoter")} ({distance_from_prom:+d}bp)',
                            confidence=confidence,
                            metadata={
                                'strength': strength,
                                'promoter': prom.get('feature_name'),
                                'distance_from_promoter_end': distance_from_prom
                            }
                        ))
                    break  # Only use strongest match

        # 3. Existing: Search upstream of bacterial markers
        hits.extend(self._detect_rbs_upstream_of_markers(features))

        return hits

    def _detect_rbs_upstream_of_markers(
        self,
        features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """
        Find RBS/Shine-Dalgarno sequences upstream of bacterial CDS starts.

        RBS must be 4-20bp upstream of start codon.
        """
        hits = []
        cds_features = [
            f for f in features
            if f.get('role', '') in ('expression_payload', 'editing_payload', 'reporter_payload', 'selection_payload', 'bacterial_marker')
            or f.get('feature_type', '').lower() == 'cds'
        ]

        for cds in cds_features:
            strand = cds.get("strand", 1)
            start = cds.get("start", 0)
            end = cds.get("end", 0)
            role = cds.get("role", "")

            # Only look for RBS near bacterial markers
            feature_name = cds.get("feature_name", "").lower()
            bacterial_marker_names = ["ampr", "kanr", "camr", "cmr", "specr", "tetr", "genr",
                                      "bla", "neo", "amp", "kan", "cat", "spec", "tet", "strep"]
            is_bacterial_marker = (
                role in ("bacterial_marker", "selection_payload") or
                any(marker in feature_name for marker in bacterial_marker_names)
            )
            if not is_bacterial_marker:
                continue

            # Get 30bp upstream of CDS start
            if strand == 1:
                upstream_start = max(0, start - 30)
                upstream_end = start
                region = self.sequence[upstream_start:upstream_end]
            else:
                upstream_start = end
                upstream_end = min(end + 30, self.seq_len)
                region = reverse_complement(self.sequence[upstream_start:upstream_end])

            # Try patterns in order of specificity (strongest first)
            found_positions = set()  # Track positions to avoid duplicates
            # Use a tuple of (start, end, strand) as unique CDS identifier
            cds_id = cds.get("instance_id") or (start, end, strand)

            # Strong RBS (5-6bp consensus)
            for match in RBS_STRONG_PATTERN.finditer(region):
                rbs_pos = match.start()
                rbs_seq = match.group()
                distance_to_start = len(region) - rbs_pos - len(rbs_seq)

                if not (4 <= distance_to_start <= 20):
                    continue

                if strand == 1:
                    genomic_start = upstream_start + rbs_pos
                    genomic_end = genomic_start + len(rbs_seq)
                else:
                    genomic_end = upstream_end - rbs_pos
                    genomic_start = genomic_end - len(rbs_seq)

                found_positions.add((genomic_start, genomic_end))
                hits.append(MotifHit(
                    motif_type="rbs",
                    name="RBS (strong)",
                    start=genomic_start,
                    end=genomic_end,
                    strand=strand,
                    sequence=rbs_seq,
                    description=f"Strong Shine-Dalgarno sequence {distance_to_start}bp upstream of CDS",
                    confidence=0.95,
                    metadata={
                        "cds_id": cds_id,
                        "distance_to_start": distance_to_start,
                        "strength": "strong",
                    }
                ))

            # Adequate RBS (4-5bp partial consensus)
            for match in RBS_ADEQUATE_PATTERN.finditer(region):
                rbs_pos = match.start()
                rbs_seq = match.group()
                distance_to_start = len(region) - rbs_pos - len(rbs_seq)

                if not (4 <= distance_to_start <= 20):
                    continue

                if strand == 1:
                    genomic_start = upstream_start + rbs_pos
                    genomic_end = genomic_start + len(rbs_seq)
                else:
                    genomic_end = upstream_end - rbs_pos
                    genomic_start = genomic_end - len(rbs_seq)

                # Skip if overlaps with stronger match
                if any(gs <= genomic_start < ge or gs < genomic_end <= ge
                       for gs, ge in found_positions):
                    continue

                found_positions.add((genomic_start, genomic_end))
                hits.append(MotifHit(
                    motif_type="rbs",
                    name="RBS (adequate)",
                    start=genomic_start,
                    end=genomic_end,
                    strand=strand,
                    sequence=rbs_seq,
                    description=f"Adequate Shine-Dalgarno sequence {distance_to_start}bp upstream of CDS",
                    confidence=0.75,
                    metadata={
                        "cds_id": cds_id,
                        "distance_to_start": distance_to_start,
                        "strength": "adequate",
                    }
                ))

            # Weak RBS (3bp minimal core) - only if no stronger match found for this CDS
            # Use position-based comparison to avoid None==None matching
            cds_has_strong = any(
                h.metadata.get("cds_id") == cds_id and cds_id is not None
                for h in hits
            ) if cds_id else False
            if not cds_has_strong:
                for match in RBS_WEAK_PATTERN.finditer(region):
                    rbs_pos = match.start()
                    rbs_seq = match.group()
                    distance_to_start = len(region) - rbs_pos - len(rbs_seq)

                    if not (4 <= distance_to_start <= 20):
                        continue

                    if strand == 1:
                        genomic_start = upstream_start + rbs_pos
                        genomic_end = genomic_start + len(rbs_seq)
                    else:
                        genomic_end = upstream_end - rbs_pos
                        genomic_start = genomic_end - len(rbs_seq)

                    # Skip if overlaps with any existing match
                    if any(gs <= genomic_start < ge or gs < genomic_end <= ge
                           for gs, ge in found_positions):
                        continue

                    found_positions.add((genomic_start, genomic_end))
                    hits.append(MotifHit(
                        motif_type="rbs",
                        name="RBS (weak)",
                        start=genomic_start,
                        end=genomic_end,
                        strand=strand,
                        sequence=rbs_seq,
                        description=f"Weak Shine-Dalgarno sequence {distance_to_start}bp upstream of CDS",
                        confidence=0.5,
                        metadata={
                            "cds_id": cds_id,
                            "distance_to_start": distance_to_start,
                            "strength": "weak",
                        }
                    ))
                    break  # Only report one weak RBS per CDS

        return hits

    def detect_start_codons_in_context(
        self,
        features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """
        Detect start codons (ATG) in biologically relevant contexts:
        1. Within a Kozak sequence (for eukaryotic genes)
        2. 4-20bp downstream of an RBS (for prokaryotic genes)
        """
        hits = []
        cds_features = [
            f for f in features
            if f.get('role', '') in ('expression_payload', 'editing_payload', 'reporter_payload', 'selection_payload', 'bacterial_marker')
            or f.get('feature_type', '').lower() == 'cds'
        ]

        # First, detect Kozak sequences and extract ATG positions from them
        kozak_hits = self.detect_kozak_sequences(cds_features)
        kozak_atg_positions = set()

        for kozak in kozak_hits:
            # Find ATG position within Kozak sequence
            kozak_seq = self.sequence[kozak.start:kozak.end].upper()
            atg_pos = kozak_seq.find('ATG')
            if atg_pos >= 0:
                genomic_atg_start = kozak.start + atg_pos
                kozak_atg_positions.add(genomic_atg_start)
                hits.append(MotifHit(
                    motif_type='start_codon',
                    name='ATG (Kozak)',
                    start=genomic_atg_start,
                    end=genomic_atg_start + 3,
                    strand=kozak.strand,
                    sequence='ATG',
                    description=f'Start codon within {kozak.name}',
                    confidence=kozak.confidence,
                    metadata={'context': 'kozak', 'kozak_strength': kozak.metadata.get('strength')}
                ))

        # Then, detect ATG 4-20bp after RBS
        rbs_hits = self.detect_rbs_sequences(features)

        for rbs in rbs_hits:
            strand = rbs.strand
            if strand == 1:
                # Forward strand: search 4-20bp after RBS end
                search_start = rbs.end + 4
                search_end = min(rbs.end + 21, self.seq_len)
                if search_start >= search_end:
                    continue
                region = self.sequence[search_start:search_end].upper()

                atg_pos = region.find('ATG')
                if atg_pos >= 0 and atg_pos <= 16:  # Within 4-20bp window
                    genomic_atg_start = search_start + atg_pos
                    if genomic_atg_start not in kozak_atg_positions:
                        hits.append(MotifHit(
                            motif_type='start_codon',
                            name='ATG (RBS)',
                            start=genomic_atg_start,
                            end=genomic_atg_start + 3,
                            strand=strand,
                            sequence='ATG',
                            description=f'Start codon {atg_pos + 4}bp after RBS',
                            confidence=rbs.confidence * 0.9,
                            metadata={'context': 'rbs', 'distance_from_rbs': atg_pos + 4}
                        ))
            else:
                # Reverse strand: search 4-20bp before RBS start
                search_end = rbs.start - 4
                search_start = max(rbs.start - 21, 0)
                if search_start >= search_end:
                    continue
                region = self.sequence[search_start:search_end].upper()

                # Look for CAT (reverse complement of ATG)
                cat_pos = region.rfind('CAT')
                if cat_pos >= 0:
                    distance = len(region) - cat_pos - 3
                    if 4 <= distance <= 20:
                        genomic_atg_start = search_start + cat_pos
                        if genomic_atg_start not in kozak_atg_positions:
                            hits.append(MotifHit(
                                motif_type='start_codon',
                                name='ATG (RBS)',
                                start=genomic_atg_start,
                                end=genomic_atg_start + 3,
                                strand=strand,
                                sequence='CAT',  # Reverse complement
                                description=f'Start codon {distance}bp after RBS (complement)',
                                confidence=rbs.confidence * 0.9,
                                metadata={'context': 'rbs', 'distance_from_rbs': distance}
                            ))

        return hits

    def detect_all_motifs(
        self,
        features: List[Dict[str, Any]]
    ) -> List[MotifHit]:
        """
        Run all motif detection and return combined results.

        Detects:
        - Kozak sequences (strong/adequate/weak)
        - RBS/Shine-Dalgarno sequences
        - Start codons in context (within Kozak or after RBS)

        NOTE: 2A peptides, stop codons, internal stops, and introns are
        handled by CDS submodule parsing, not motif detection.
        """
        all_hits = []

        # Separate feature types
        promoters = [
            f for f in features
            if 'promoter' in f.get('role', '').lower()
        ]
        cds_features = [
            f for f in features
            if f.get('role', '') in ('expression_payload', 'editing_payload', 'reporter_payload', 'selection_payload', 'bacterial_marker')
            or f.get('feature_type', '').lower() == 'cds' or any(x in f.get('feature_name', '').lower() for x in ['cas9', 'cas12', 'puro', 'bla', 'neo', 'hygro', 'zeo'])
        ]

        # Run detections
        all_hits.extend(self.detect_kozak_sequences(cds_features))
        all_hits.extend(self.detect_rbs_sequences(features))  # Pass all features for promoter lookup
        all_hits.extend(self.detect_start_codons_in_context(features))  # Context-aware start codons

        if promoters and cds_features:
            all_hits.extend(self.detect_start_loss(promoters, cds_features))

        return all_hits


# =============================================================================
# MODULE BOUNDARY DETECTOR CLASS
# =============================================================================

class ModuleBoundaryDetector:
    """Use hierarchical_annotator functions for module parsing"""

    def __init__(self, sequence: str, features: List[Dict[str, Any]], circular: bool):
        self.sequence = sequence
        self.features = features
        self.circular = circular
        self.seq_len = len(sequence)

    def detect_expression_modules(self) -> List[ModuleBoundary]:
        """
        Detect pol2/pol3 expression module boundaries.
        
        Uses promoter and terminator/polyA features to define module boundaries.
        """
        boundaries = []
        
        # Find promoters by type
        pol2_promoters = [
            f for f in self.features
            if f.get('role') == 'pol2_promoter' or
               f.get('kb_polymerase_class') == 'pol_ii'
        ]
        
        pol3_promoters = [
            f for f in self.features
            if f.get('role') == 'pol3_promoter' or
               f.get('kb_polymerase_class') == 'pol_iii'
        ]
        
        # Find terminators
        terminators = [
            f for f in self.features
            if f.get('role') in ('terminator', 'polya') or
               f.get('kb_feature_class') in ('terminator', 'polyA_signal')
        ]
        
        # Create boundaries for Pol II promoters
        for promoter in pol2_promoters:
            boundaries.append(ModuleBoundary(
                module_type='pol2_expression',
                boundary_type='start',
                position=promoter.get('start', 0),
                strand=promoter.get('strand', 1),
                feature_id=promoter.get('instance_id'),
                feature_name=promoter.get('feature_name'),
            ))
        
        # Create boundaries for Pol III promoters
        for promoter in pol3_promoters:
            boundaries.append(ModuleBoundary(
                module_type='pol3_expression',
                boundary_type='start',
                position=promoter.get('start', 0),
                strand=promoter.get('strand', 1),
                feature_id=promoter.get('instance_id'),
                feature_name=promoter.get('feature_name'),
            ))
        
        # Create end boundaries from terminators
        for term in terminators:
            boundaries.append(ModuleBoundary(
                module_type='pol2_expression',  # Default to pol2
                boundary_type='end',
                position=term.get('end', 0),
                strand=term.get('strand', 1),
                feature_id=term.get('instance_id'),
                feature_name=term.get('feature_name'),
            ))
        
        return boundaries

    def detect_transcript_boundaries(self) -> List[TranscriptBoundary]:
        """
        Detect transcript start/end from motifs.
        
        Start: promoter end + Kozak/ATG
        End: stop codon + polyA signal
        """
        boundaries = []
        
        # Find promoters and terminators
        promoters = [
            f for f in self.features
            if 'promoter' in f.get('role', '')
        ]
        
        terminators = [
            f for f in self.features
            if f.get('role') in ('terminator', 'polya')
        ]
        
        # Sort by position
        promoters.sort(key=lambda f: f.get('start', 0))
        terminators.sort(key=lambda f: f.get('end', 0))
        
        # Match promoters to terminators on same strand
        for promoter in promoters:
            strand = promoter.get('strand', 1)
            promoter_end = promoter.get('end', 0)
            
            # Find nearest downstream terminator on same strand
            downstream_terms = [
                t for t in terminators
                if t.get('strand', 1) == strand and
                   ((strand == 1 and t.get('start', 0) > promoter_end) or
                    (strand == -1 and t.get('end', 0) < promoter.get('start', 0)))
            ]
            
            if downstream_terms:
                if strand == 1:
                    nearest_term = min(downstream_terms, key=lambda t: t.get('start', 0))
                    boundaries.append(TranscriptBoundary(
                        start_position=promoter_end,
                        end_position=nearest_term.get('end', 0),
                        strand=strand,
                        promoter_id=promoter.get('instance_id'),
                        terminator_id=nearest_term.get('instance_id'),
                    ))
                else:
                    nearest_term = max(downstream_terms, key=lambda t: t.get('end', 0))
                    boundaries.append(TranscriptBoundary(
                        start_position=nearest_term.get('start', 0),
                        end_position=promoter.get('start', 0),
                        strand=strand,
                        promoter_id=promoter.get('instance_id'),
                        terminator_id=nearest_term.get('instance_id'),
                    ))
        
        return boundaries

    def detect_protein_boundaries(self) -> List[ProteinBoundary]:
        """
        Detect protein-level CDS boundaries (removing introns/2A).
        
        Each protein after 2A/IRES becomes a separate protein boundary.
        """
        boundaries = []
        
        # Find CDS features
        cds_features = [
            f for f in self.features
            if f.get('role') in ('expression_payload', 'editing_payload', 'reporter_payload', 'selection_payload')
            or f.get('feature_type', '').lower() == 'cds' or any(x in f.get('feature_name', '').lower() for x in ['cas9', 'cas12', 'puro', 'bla', 'neo', 'hygro', 'zeo'])
        ]
        
        # Find 2A peptides and IRES elements
        linkers = [
            f for f in self.features
            if any(x in f.get('feature_name', '').lower() for x in ['2a', 'ires', 't2a', 'p2a', 'e2a', 'f2a'])
        ]
        
        # Group CDS by strand
        for strand in [1, -1]:
            strand_cds = [c for c in cds_features if c.get('strand', 1) == strand]
            strand_cds.sort(key=lambda c: c.get('start', 0))
            
            if not strand_cds:
                continue
            
            # Check for 2A/IRES between CDS
            current_start = strand_cds[0].get('start', 0)
            current_ids = [strand_cds[0].get('instance_id')]
            has_2a = False
            
            for i in range(len(strand_cds) - 1):
                cds1 = strand_cds[i]
                cds2 = strand_cds[i + 1]
                
                # Check if there's a linker between them
                gap_start = cds1.get('end', 0)
                gap_end = cds2.get('start', 0)
                
                linker_between = any(
                    l.get('start', 0) >= gap_start and l.get('end', 0) <= gap_end
                    for l in linkers
                )
                
                if linker_between:
                    # Create protein boundary up to this point
                    boundaries.append(ProteinBoundary(
                        start_position=current_start,
                        end_position=cds1.get('end', 0),
                        strand=strand,
                        cds_feature_ids=current_ids.copy(),
                        has_2a_split=True,
                    ))
                    
                    # Start new protein
                    current_start = cds2.get('start', 0)
                    current_ids = [cds2.get('instance_id')]
                    has_2a = True
                else:
                    current_ids.append(cds2.get('instance_id'))
            
            # Add final protein
            if strand_cds:
                boundaries.append(ProteinBoundary(
                    start_position=current_start,
                    end_position=strand_cds[-1].get('end', 0),
                    strand=strand,
                    cds_feature_ids=current_ids,
                    has_2a_split=has_2a,
                ))
        
        return boundaries


# =============================================================================
# GENBANK ANNOTATION FUNCTIONS
# =============================================================================

def add_motif_annotations_to_genbank(
    record: SeqRecord,
    motifs: List[MotifHit]
) -> SeqRecord:
    """
    Add detected motifs as feature annotations to plasmid.
    
    Adds features with type='motif' and appropriate qualifiers.
    """
    for motif in motifs:
        # Create feature location
        location = FeatureLocation(
            motif.start,
            motif.end,
            strand=motif.strand
        )
        
        # Create feature
        feature = SeqFeature(
            location=location,
            type='motif',
            qualifiers={
                'label': motif.name,
                'motif_type': motif.motif_type,
                'sequence': motif.sequence,
                'note': motif.description,
            }
        )
        
        # Add metadata as additional qualifiers
        if motif.metadata:
            for key, value in motif.metadata.items():
                if isinstance(value, (str, int, float)):
                    feature.qualifiers[key] = str(value)
        
        record.features.append(feature)
    
    return record


# =============================================================================
# CDS SEQUENCE ANALYSIS (from plasmid_analyzer.py)
# =============================================================================

def analyze_cds_sequence(
    sequence: str,
    start: int,
    end: int,
    strand: int,
    circular: bool = True
) -> Dict[str, Any]:
    """
    Analyze CDS sequence for start/stop codons, Kozak strength, internal stops, 2A peptides.
    
    Returns dict with analysis results.
    """
    result = {
        'has_start_codon': False,
        'has_stop_codon': False,
        'kozak_strength': 'unknown',
        'internal_stops': 0,
        'contains_2a': False,
        'detected_2a_peptides': [],
        'codon_position_5p': 0,
    }
    
    # Extract feature sequence
    subseq = get_sequence_slice(sequence, start, end, strand, circular)
    if not subseq or len(subseq) < 3:
        return result
    
    subseq = subseq.upper()
    
    # Check for start codon
    start_pos = find_start_codon(subseq)
    result['has_start_codon'] = start_pos is not None and start_pos < 10
    
    if start_pos is not None:
        result['codon_position_5p'] = start_pos % 3
    
    # Check for stop codon at end
    if len(subseq) >= 3:
        last_codon = subseq[-3:]
        result['has_stop_codon'] = last_codon in STOP_CODONS
    
    # Find internal stop codons (excluding the terminal one)
    frame = result['codon_position_5p'] if result['has_start_codon'] else 0
    stops = find_stop_codons(subseq, frame)
    
    # Exclude terminal stop
    if result['has_stop_codon'] and stops and stops[-1] >= len(subseq) - 3:
        stops = stops[:-1]
    
    result['internal_stops'] = len(stops)
    
    # Translate and look for 2A peptides
    protein = translate_sequence(subseq, frame)
    detected_2a = detect_2a_peptides(protein)
    result['contains_2a'] = len(detected_2a) > 0
    result['detected_2a_peptides'] = [name for name, _, _ in detected_2a]
    
    # Assess Kozak strength if we have start codon
    if start_pos is not None and start_pos >= 3:
        # Get context around ATG
        context_start = max(0, start_pos - 6)
        context_end = min(len(subseq), start_pos + 7)
        context = subseq[context_start:context_end]
        
        if KOZAK_STRONG_PATTERN.search(context):
            result['kozak_strength'] = 'strong'
        elif KOZAK_ADEQUATE_PATTERN.search(context):
            result['kozak_strength'] = 'adequate'
        elif KOZAK_WEAK_PATTERN.search(context):
            result['kozak_strength'] = 'weak'
    
    return result
