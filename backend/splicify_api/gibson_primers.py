from __future__ import annotations

from typing import List, Union, Optional, Literal, Dict, Any, Tuple
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field, ConfigDict
import math
import primer3

from .utils import (
    normalize_dna,
    reverse_complement,
    safe_tm,
    ensure_session_id,
    is_valid_dna,
)

router = APIRouter(tags=["gibson_primers"])


# ----------------------------
# Models
# ----------------------------

class GibsonFragment(BaseModel):
    name: str
    sequence: str


class GibsonRequest(BaseModel):
    model_config = ConfigDict(populate_by_name=True, extra="ignore")

    fragments: List[Union[GibsonFragment, str]]
    assembly: Literal["linear", "circular"] = "circular"
    
    anneal_tm_target: float = 60.0
    overlap_tm_min: float = 52.0
    overlap_tm_max: float = 70.0
    
    anneal_min_len: int = 18
    anneal_max_len: int = 30
    overlap_min_len: int = 20
    overlap_max_len: int = 40
    
    # For existing overlap detection
    existing_overlap_min: int = 15
    existing_overlap_max: int = 50
    
    overlap_len: Optional[int] = Field(default=None, ge=10, le=60)
    homology_window: int = Field(40, ge=10, le=200)
    homology_min: int = Field(15, ge=8, le=60)
    
    session_id: Optional[str] = Field(default=None, alias="sessionId")
    include_ai_explanation: Optional[bool] = Field(default=None)

# ----------------------------
# Thermodynamic & Sequence Analysis
# ----------------------------

class ThermodynamicCalculator:
    """Thermodynamic calculations using primer3"""
    
    @staticmethod
    def calculate_tm(sequence: str, salt_conc: float = 50.0, dna_conc: float = 250.0) -> float:
        """Calculate Tm using primer3"""
        if not sequence or len(sequence) < 2:
            return 0.0
        try:
            tm = primer3.calcTm(sequence, mv_conc=salt_conc, dv_conc=0.0, dntp_conc=0.0, dna_conc=dna_conc)
            return float(tm)
        except Exception:
            return 0.0
    
    @classmethod
    def calculate_hairpin_dg(cls, sequence: str) -> float:
        if not sequence:
            return 0.0
        try:
            result = primer3.calcHairpin(sequence)
            return float(result.dg / 1000.0)
        except Exception:
            return 0.0
    
    @classmethod
    def calculate_homodimer_dg(cls, sequence: str) -> float:
        if not sequence:
            return 0.0
        try:
            result = primer3.calcHomodimer(sequence)
            return float(result.dg / 1000.0)
        except Exception:
            return 0.0
    
    @classmethod
    def calculate_heterodimer_dg(cls, seq1: str, seq2: str) -> float:
        if not seq1 or not seq2:
            return 0.0
        try:
            result = primer3.calcHeterodimer(seq1, seq2)
            return float(result.dg / 1000.0)
        except Exception:
            return 0.0
    
    @classmethod
    def calculate_end_stability(cls, sequence: str) -> float:
        if len(sequence) < 5:
            return 0.0
        end_seq = sequence[-5:]
        result = primer3.calcHomodimer(end_seq)
        return float(result.dg / 1000.0)


class SequenceAnalyzer:
    """Analyze sequences for problematic features"""
    
    @staticmethod
    def find_homopolymers(sequence: str) -> List[Tuple[int, int, str]]:
        homopolymers = []
        seq = sequence.upper()
        i = 0
        while i < len(seq):
            current_base = seq[i]
            run_length = 1
            while i + run_length < len(seq) and seq[i + run_length] == current_base:
                run_length += 1
            if run_length >= 4:
                homopolymers.append((i, run_length, current_base))
            i += run_length
        return homopolymers
    
    @staticmethod
    def find_dinucleotide_repeats(sequence: str) -> List[Tuple[int, int, str]]:
        repeats = []
        seq = sequence.upper()
        i = 0
        while i < len(seq) - 1:
            dinuc = seq[i:i+2]
            repeat_count = 1
            j = i + 2
            while j + 1 < len(seq) and seq[j:j+2] == dinuc:
                repeat_count += 1
                j += 2
            if repeat_count >= 4:
                repeats.append((i, repeat_count, dinuc))
                i = j
            else:
                i += 1
        return repeats
    
    @staticmethod
    def sequence_similarity(seq1: str, seq2: str) -> float:
        if not seq1 or not seq2:
            return 0.0
        min_len = min(len(seq1), len(seq2))
        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
        return matches / max(len(seq1), len(seq2))
    
    @staticmethod
    def check_3prime_composition(sequence: str) -> Dict[str, Any]:
        if not sequence:
            return {'has_gc_clamp': False, 'ends_with_t': False, 'gc_percent': 0.0}
        end = sequence[-5:] if len(sequence) >= 5 else sequence
        gc_count = end.count('G') + end.count('C')
        return {
            'has_gc_clamp': sequence[-1] in 'GC',
            'ends_with_t': sequence[-1] == 'T',
            'gc_percent': gc_count / len(end) if end else 0.0
        }

class UniquenessCalculator:
    """Calculate uniqueness score based on off-target Tm differences"""

    def __init__(self):
        self.thermo_calc = ThermodynamicCalculator()
        self.seq_analyzer = SequenceAnalyzer()

    def calculate_uniqueness_score(
        self,
        overlap_seq: str,
        construct_seq: str,
        similarity_threshold: float = 0.66,
        tm_window_c: float = 10.0,
        penalty_per_site: float = 10.0,
    ) -> Tuple[float, List[Dict[str, Any]]]:
        """Estimate overlap uniqueness by scanning the full construct for similar sites,
        including origin-spanning windows on circular constructs.
        """
        if not overlap_seq or not construct_seq:
            return 100.0, []

        overlap_seq = overlap_seq.upper()
        construct_seq = construct_seq.upper()

        target_tm = self.thermo_calc.calculate_tm(overlap_seq)
        overlap_len = len(overlap_seq)

        n = len(construct_seq)
        if n < overlap_len:
            return 100.0, []

        # Circular scan string: allows windows to wrap the origin
        circular_seq = construct_seq + construct_seq[:overlap_len - 1]

        candidates: List[Dict[str, Any]] = []
        exact_match_positions: List[int] = []

        # Scan exactly n start positions, so every circular window is checked once
        for i in range(n):
            window = circular_seq[i:i + overlap_len]
            similarity = self.seq_analyzer.sequence_similarity(overlap_seq, window)
            if similarity < similarity_threshold:
                continue

            spans_origin = (i + overlap_len > n)
            is_exact = (window == overlap_seq)
            if is_exact:
                exact_match_positions.append(i)

            off_tm = self.thermo_calc.calculate_tm(window)
            tm_delta = float(target_tm - off_tm)  # positive => off-target is lower Tm

            candidates.append({
                "position": i,                 # 0..n-1, start index on circular construct
                "spans_origin": spans_origin,  # True if window wraps end->start
                "sequence": window,
                "similarity": float(similarity),
                "tm": float(off_tm),
                "tm_delta_from_target": tm_delta,
                "is_exact_match": is_exact,
            })

        # No >=66% identity sites => uniqueness 100
        if not candidates:
            return 100.0, []

        # Ignore one exact match occurrence (assumed on-target overlap).
        ignored_exact_pos = min(exact_match_positions) if exact_match_positions else None

        off_targets: List[Dict[str, Any]] = []
        for c in candidates:
            if c["is_exact_match"] and ignored_exact_pos is not None and c["position"] == ignored_exact_pos:
                continue
            off_targets.append(c)

        # Only the on-target exact occurrence existed
        if not off_targets:
            return 100.0, []

        score = 100.0
        for ot in off_targets:
            off_tm = ot["tm"]

            # High-risk if within tm_window_c below target (or higher than target).
            if off_tm >= (target_tm - tm_window_c):
                ot["risk"] = "high"
                ot["rule_trigger"] = f"off_tm >= target_tm - {tm_window_c}"
                score = 0.0
                break
            else:
                ot["risk"] = "moderate"
                ot["rule_trigger"] = f"off_tm < target_tm - {tm_window_c}"
                score -= penalty_per_site

        score = max(0.0, min(100.0, score))

        for ot in off_targets:
            ot["target_tm"] = float(target_tm)

        return score, off_targets


class OverlapScorer:
    """Comprehensive overlap quality scoring"""
    
    def __init__(self):
        self.thermo_calc = ThermodynamicCalculator()
        self.seq_analyzer = SequenceAnalyzer()
        self.uniqueness_calc = UniquenessCalculator()
    
    def score_overlap(
        self,
        overlap_seq: str,
        all_overlaps: List[str] = None,
        construct_seq: str = None
    ) -> Dict[str, Any]:
        if not overlap_seq:
            return self._empty_score()
        
        overlap_seq = overlap_seq.upper()
        
        tm = self.thermo_calc.calculate_tm(overlap_seq)
        gc_content = (overlap_seq.count('G') + overlap_seq.count('C')) / len(overlap_seq)
        hairpin_dg = self.thermo_calc.calculate_hairpin_dg(overlap_seq)
        self_dimer_dg = self.thermo_calc.calculate_homodimer_dg(overlap_seq)
        end_dg = self.thermo_calc.calculate_end_stability(overlap_seq)
        
        homopolymers = self.seq_analyzer.find_homopolymers(overlap_seq)
        repeats = self.seq_analyzer.find_dinucleotide_repeats(overlap_seq)
        end_comp = self.seq_analyzer.check_3prime_composition(overlap_seq)
        
        length_score = self._score_length(len(overlap_seq))
        tm_score = self._score_tm(tm)
        gc_score = self._score_gc(gc_content)
        hairpin_score = self._score_hairpin(hairpin_dg)
        dimer_score = self._score_dimer(self_dimer_dg)
        end_score = self._score_end_stability(end_dg)
        homopoly_score = self._score_homopolymers(homopolymers, overlap_seq)
        repeat_score = self._score_repeats(repeats)
        
        cross_dimer_score = 100.0
        worst_cross_dg = 0.0
        if all_overlaps:
            for other in all_overlaps:
                if other != overlap_seq:
                    cross_dg = self.thermo_calc.calculate_heterodimer_dg(overlap_seq, other)
                    if cross_dg < worst_cross_dg:
                        worst_cross_dg = cross_dg
            cross_dimer_score = self._score_dimer(worst_cross_dg)
        
        tm_uniformity_score = 100.0
        if all_overlaps and len(all_overlaps) > 1:
            tms = [self.thermo_calc.calculate_tm(seq) for seq in all_overlaps]
            tm_range = max(tms) - min(tms)
            tm_uniformity_score = 100.0 * (1.0 - min(tm_range / 8.0, 1.0))
        
        uniqueness_score = 100.0
        off_targets = []
        if construct_seq:
            uniqueness_score, off_targets = self.uniqueness_calc.calculate_uniqueness_score(
                overlap_seq, construct_seq
            )
        
        total_score = (
            length_score * 0.12 +
            tm_score * 0.15 +
            tm_uniformity_score * 0.12 +
            gc_score * 0.08 +
            hairpin_score * 0.10 +
            dimer_score * 0.08 +
            cross_dimer_score * 0.08 +
            end_score * 0.05 +
            homopoly_score * 0.04 +
            repeat_score * 0.03 +
            uniqueness_score * 0.15
        )
        
        warnings = []
        if tm < 48:
            warnings.append(f"Low Tm ({tm:.1f}°C)")
        if hairpin_dg < -3:
            warnings.append(f"Hairpin ΔG={hairpin_dg:.1f}")
        if self_dimer_dg < -6:
            warnings.append(f"Self-dimer ΔG={self_dimer_dg:.1f}")
        if worst_cross_dg < -6:
            warnings.append(f"Cross-dimer ΔG={worst_cross_dg:.1f}")
        if not end_comp['has_gc_clamp']:
            warnings.append("No GC clamp")
        if uniqueness_score < 50:
            warnings.append(f"Low uniqueness ({uniqueness_score:.0f}/100)")
        if homopolymers:
            max_homo = max(h[1] for h in homopolymers)
            if max_homo > 4:
                warnings.append(f"{max_homo}bp homopolymer")
        
        return {
            'total_score': round(total_score, 1),
            'overlap_length': len(overlap_seq),
            'tm': round(tm, 1),
            'gc_content': round(gc_content, 3),
            'hairpin_dg': round(hairpin_dg, 1),
            'self_dimer_dg': round(self_dimer_dg, 1),
            'cross_dimer_dg': round(worst_cross_dg, 1),
            'end_stability_dg': round(end_dg, 1),
            'homopolymer_count': len(homopolymers),
            'max_homopolymer_length': max((h[1] for h in homopolymers), default=0),
            'dinuc_repeat_count': len(repeats),
            'has_gc_clamp': end_comp['has_gc_clamp'],
            'ends_with_t': end_comp['ends_with_t'],
            'uniqueness_score': round(uniqueness_score, 1),
            'off_target_count': len(off_targets),
            'length_score': round(length_score, 1),
            'tm_score': round(tm_score, 1),
            'tm_uniformity_score': round(tm_uniformity_score, 1),
            'gc_score': round(gc_score, 1),
            'hairpin_score': round(hairpin_score, 1),
            'dimer_score': round(dimer_score, 1),
            'cross_dimer_score': round(cross_dimer_score, 1),
            'end_stability_score': round(end_score, 1),
            'homopolymer_score': round(homopoly_score, 1),
            'repeat_score': round(repeat_score, 1),
            'warnings': warnings,
        }
    
    def _empty_score(self) -> Dict[str, Any]:
        return {'total_score': 0.0, 'overlap_length': 0, 'warnings': ['No overlap']}
    
    def _score_length(self, length: int, target: int = 30) -> float:
        if length < 20 or length > 40:
            return 0.0
        deviation = abs(length - target)
        return 100.0 * math.exp(-(deviation ** 2) / (2 * (5 ** 2)))
    
    def _score_tm(self, tm: float, target: float = 62.0) -> float:
        if tm < 48:
            return max(0.0, 100.0 - (48 - tm) * 10)
        if tm > 75:
            return max(0.0, 100.0 - (tm - 75) * 5)
        deviation = abs(tm - target)
        return 100.0 * (1.0 - deviation / 14.0)
    
    def _score_gc(self, gc: float) -> float:
        if gc < 0.40 or gc > 0.60:
            deviation = max(0.40 - gc, gc - 0.60, 0)
            return max(0.0, 100.0 - deviation * 200)
        deviation = abs(gc - 0.50)
        return 100.0 * (1.0 - deviation / 0.10)
    
    def _score_hairpin(self, dg: float) -> float:
        if dg >= 0:
            return 100.0
        if dg >= -3:
            return max(50.0, 100.0 - abs(dg) * 10)
        return max(0.0, 50.0 - (abs(dg) - 3) * 20)
    
    def _score_dimer(self, dg: float) -> float:
        if dg >= 0:
            return 100.0
        if dg >= -6:
            return max(60.0, 100.0 - abs(dg) * 8)
        return max(0.0, 60.0 - (abs(dg) - 6) * 15)
    
    def _score_end_stability(self, dg: float) -> float:
        if dg >= -9:
            return 100.0
        return max(0.0, 100.0 - (abs(dg) - 9) * 8)
    
    def _score_homopolymers(self, homopolymers: List, sequence: str) -> float:
        if not homopolymers:
            return 100.0
        score = 100.0
        for start, length, base in homopolymers:
            is_3prime = (start + length >= len(sequence) - 3)
            max_len = 3 if is_3prime else 4
            if length > max_len:
                penalty = (length - max_len) * 15
                if is_3prime:
                    penalty *= 2
                score -= penalty
        return max(0.0, score)
    
    def _score_repeats(self, repeats: List) -> float:
        if not repeats:
            return 100.0
        score = 100.0
        for start, repeat_count, dinuc in repeats:
            if repeat_count > 4:
                score -= (repeat_count - 4) * 12
        return max(0.0, score)


# ----------------------------
# Helper Functions
# ----------------------------
def _max_suffix_prefix_overlap(left: str, right: str, max_k: int) -> int:
    left = normalize_dna(left)
    right = normalize_dna(right)
    max_k = min(max_k, len(left), len(right))
    for k in range(max_k, 0, -1):
        if left[-k:] == right[:k]:
            return k
    return 0


def assemble_construct_with_existing_overlaps(
    frags: List[GibsonFragment],
    assembly: str,
    existing_overlaps: Dict[tuple, Dict[str, Any]],
) -> Tuple[str, List[Dict[str, Any]]]:
    """
    Assemble construct sequence without duplicating pre-existing overlaps.
    - For each junction A->B with an existing overlap k, we trim k bases off the *prefix of B*.
    - For circular, we do NOT append the first fragment again; we just avoid duplicating internal overlaps.
    Returns (assembled_seq, annotations).
    """
    if not frags:
        return "", []

    n = len(frags)

    def next_index(i: int) -> Optional[int]:
        if i < n - 1:
            return i + 1
        return 0 if assembly == "circular" else None

    # how much to trim from the start of each fragment due to overlap from prev->this
    trim_prefix: Dict[str, int] = {f.name: 0 for f in frags}

    for i, f in enumerate(frags):
        ni = next_index(i)
        if ni is None:
            continue
        g = frags[ni]
        key = (f.name, g.name)
        ov = existing_overlaps.get(key)
        if ov:
            k = int(ov.get("overlap_length") or 0)
            if k > 0:
                # Trim overlap off the beginning of the RIGHT fragment
                trim_prefix[g.name] = max(trim_prefix[g.name], k)

    # Assemble + annotate
    assembled_parts: List[str] = []
    annotations: List[Dict[str, Any]] = []
    pos = 0

    for f in frags:
        k = trim_prefix.get(f.name, 0)
        seq = f.sequence[k:] if k > 0 else f.sequence

        assembled_parts.append(seq)

        annotations.append({
            "name": f.name,
            "start": pos,
            "end": pos + len(seq),
            "direction": 1,
            "trimmed_prefix_bp": k,
        })
        pos += len(seq)

    return "".join(assembled_parts), annotations


def normalize_fragments(frags: List[Union[GibsonFragment, str]]) -> List[GibsonFragment]:
    out: List[GibsonFragment] = []
    for i, f in enumerate(frags):
        if isinstance(f, GibsonFragment):
            name = f.name or f"Fragment_{i+1}"
            seq = normalize_dna(f.sequence)
        elif isinstance(f, str):
            name = f"Fragment_{i+1}"
            seq = normalize_dna(f)
        else:
            raise HTTPException(status_code=400, detail="Invalid fragments format")

        if not is_valid_dna(seq):
            raise HTTPException(status_code=400, detail=f"Invalid DNA in fragment {name}")

        out.append(GibsonFragment(name=name, sequence=seq))
    return out


def find_existing_overlap(
    left_seq: str,
    right_seq: str,
    min_overlap: int = 15,
    max_overlap: int = 50
) -> Optional[Dict[str, Any]]:
    """
    Find existing suffix/prefix overlap (15-50bp range for Gibson)
    """
    left = normalize_dna(left_seq)
    right = normalize_dna(right_seq)
    if not left or not right:
        return None

    # Check range 15-50bp
    best_k = 0
    for k in range(min_overlap, min(max_overlap + 1, len(left) + 1, len(right) + 1)):
        if left[-k:] == right[:k]:
            best_k = k

    if best_k >= min_overlap:
        seq = right[:best_k]
        thermo_calc = ThermodynamicCalculator()
        return {
            "overlap_length": best_k,
            "overlap_sequence": seq,
            "overlap_tm": thermo_calc.calculate_tm(seq),
        }
    return None


def design_spanning_overlap(
    left_frag_seq: str,
    right_frag_seq: str,
    min_len: int,
    max_len: int,
    tm_min: float,
    tm_max: float,
    all_overlaps: List[str] = None,
    construct_seq: str = None
) -> Tuple[Optional[str], Optional[int], Optional[int], Optional[Dict[str, Any]]]:
    """
    Design overlap that SPANS the junction between two fragments.
    
    Returns: (overlap_seq, left_bp, right_bp, score_details)
    where:
    - overlap_seq: the full overlap sequence
    - left_bp: number of bp from left fragment
    - right_bp: number of bp from right fragment
    - score_details: scoring metrics
    
    Example: If overlap is 30bp with 15bp from each:
    - overlap_seq = left[-15:] + right[:15]
    - left_bp = 15
    - right_bp = 15
    """
    left = normalize_dna(left_frag_seq)
    right = normalize_dna(right_frag_seq)
    
    if not left or not right:
        return None, None, None, None
    
    scorer = OverlapScorer()
    thermo_calc = ThermodynamicCalculator()
    
    best_overlap = None
    best_left_bp = None
    best_right_bp = None
    best_score_val = -1
    best_score_details = None
    
    # Try different split points across the junction
    # For a 30bp overlap, try: (30,0), (29,1), (28,2), ..., (15,15), ..., (1,29), (0,30)
    for total_len in range(min_len, max_len + 1):
        for left_bp in range(0, total_len + 1):
            right_bp = total_len - left_bp
            
            # Must have at least some from each side for most cases
            # But allow edge cases like (total_len, 0) or (0, total_len)
            if left_bp > len(left) or right_bp > len(right):
                continue
            
            # Build candidate overlap
            left_portion = left[-left_bp:] if left_bp > 0 else ""
            right_portion = right[:right_bp] if right_bp > 0 else ""
            candidate = left_portion + right_portion
            
            if len(candidate) != total_len:
                continue
            
            # Calculate Tm
            tm = thermo_calc.calculate_tm(candidate)
            
            # Pre-filter on Tm
            if tm < tm_min - 5 or tm > tm_max + 5:
                continue
            
            # Score
            score_dict = scorer.score_overlap(
                candidate,
                all_overlaps=all_overlaps,
                construct_seq=construct_seq
            )
            score = score_dict['total_score']
            
            if score > best_score_val:
                best_score_val = score
                best_overlap = candidate
                best_left_bp = left_bp
                best_right_bp = right_bp
                best_score_details = score_dict
    
    return best_overlap, best_left_bp, best_right_bp, best_score_details


def design_annealing_sequence(
    template: str,
    min_len: int,
    max_len: int,
    target_tm: float
) -> Tuple[Optional[str], Optional[float]]:
    """Design annealing sequence targeting specific Tm"""
    template = normalize_dna(template)
    if not template:
        return None, None
    
    thermo_calc = ThermodynamicCalculator()
    
    max_len = min(max_len, len(template))
    min_len = min(min_len, max_len)
    
    best_seq = None
    best_tm = None
    best_score = None
    
    for L in range(min_len, max_len + 1):
        cand = template[:L]
        tm = thermo_calc.calculate_tm(cand)
        
        score = abs(tm - target_tm)
        
        if best_score is None or score < best_score:
            best_score = score
            best_seq = cand
            best_tm = tm
    
    return best_seq, best_tm


# ----------------------------
# Main Endpoint
# ----------------------------

@router.post("/design-gibson-primers")
def design_gibson_primers(req: GibsonRequest):
    """
    Gibson Assembly primer design with proper junction-spanning overlaps
    
    Algorithm:
    1. Check for existing overlaps (15-50bp) at each junction
    2. For junctions without existing overlaps, design spanning overlaps
    3. Design primers:
       - If both ends have existing overlap: no primers needed
       - If one end existing: just annealing on that side, full primer on other
       - If no existing: full primers both sides
    4. Extensions come from the OTHER fragment at the junction
    """
    frags = normalize_fragments(req.fragments)
    if len(frags) < 2:
        raise HTTPException(status_code=400, detail="Need at least 2 fragments")

    sid = ensure_session_id(req.session_id)

    overlap_max_len = int(req.overlap_max_len)
    if req.overlap_len is not None:
        overlap_max_len = min(overlap_max_len, int(req.overlap_len))

    n = len(frags)
    
    def prev_index(i: int) -> Optional[int]:
        if i > 0:
            return i - 1
        return (n - 1) if req.assembly == "circular" else None

    def next_index(i: int) -> Optional[int]:
        if i < n - 1:
            return i + 1
        return 0 if req.assembly == "circular" else None

    # ==== STEP 1: Check for existing overlaps (15-50bp) ====
    existing_overlaps = {}  # (from_name, to_name) -> overlap_info
    junctions_needing_design = []
    
    for i, frag in enumerate(frags):
        ni = next_index(i)
        if ni is not None:
            next_frag = frags[ni]
            
            overlap = find_existing_overlap(
                frag.sequence,
                next_frag.sequence,
                min_overlap=int(req.existing_overlap_min),
                max_overlap=int(req.existing_overlap_max)
            )
            
            if overlap and overlap["overlap_length"] >= req.existing_overlap_min:
                existing_overlaps[(frag.name, next_frag.name)] = overlap
            else:
                junctions_needing_design.append((i, ni, frag.name, next_frag.name))
    
    # Build initial construct for uniqueness checking
    initial_construct, _ = assemble_construct_with_existing_overlaps(
        frags=frags,
        assembly=req.assembly,
        existing_overlaps=existing_overlaps,
    )

    
    # ==== STEP 2: Design spanning overlaps ====
    designed_overlaps = {}  # (from_name, to_name) -> (overlap_seq, left_bp, right_bp, score)
    all_designed_seqs = []
    
    # First pass: design
    for from_idx, to_idx, from_name, to_name in junctions_needing_design:
        left_frag = frags[from_idx]
        right_frag = frags[to_idx]
        
        overlap_seq, left_bp, right_bp, score_details = design_spanning_overlap(
            left_frag.sequence,
            right_frag.sequence,
            min_len=int(req.overlap_min_len),
            max_len=overlap_max_len,
            tm_min=float(req.overlap_tm_min),
            tm_max=float(req.overlap_tm_max),
            all_overlaps=None,
            construct_seq=None  # skip uniqueness scan in first pass; second pass handles it
        )
        
        if overlap_seq:
            designed_overlaps[(from_name, to_name)] = (overlap_seq, left_bp, right_bp, score_details)
            all_designed_seqs.append(overlap_seq)
    
    # Second pass: re-score with cross-dimer and Tm uniformity
    scorer = OverlapScorer()
    for key, (overlap_seq, left_bp, right_bp, _) in designed_overlaps.items():
        final_score = scorer.score_overlap(
            overlap_seq,
            all_overlaps=all_designed_seqs,
            construct_seq=initial_construct
        )
        designed_overlaps[key] = (overlap_seq, left_bp, right_bp, final_score)
    
    # ==== STEP 3: Build junctions list ====
    junctions = []
    
    for i, frag in enumerate(frags):
        ni = next_index(i)
        if ni is not None:
            next_frag = frags[ni]
            junction_key = (frag.name, next_frag.name)
            
            if junction_key in existing_overlaps:
                overlap_info = existing_overlaps[junction_key]
                junctions.append({
                    "from": frag.name,
                    "to": next_frag.name,
                    "overlap_sequence": overlap_info["overlap_sequence"],
                    "overlap_length": overlap_info["overlap_length"],
                    "overlap_tm": overlap_info["overlap_tm"],
                    "left_bp": overlap_info["overlap_length"],  # All from left
                    "right_bp": 0,  # None from right (existing homology)
                    "source": "existing"
                })
            elif junction_key in designed_overlaps:
                overlap_seq, left_bp, right_bp, score_dict = designed_overlaps[junction_key]
                junctions.append({
                    "from": frag.name,
                    "to": next_frag.name,
                    "overlap_sequence": overlap_seq,
                    "overlap_length": len(overlap_seq),
                    "overlap_tm": score_dict.get('tm'),
                    "overlap_score": score_dict.get('total_score'),
                    "left_bp": left_bp,
                    "right_bp": right_bp,
                    "source": "designed",
                    **{f"overlap_{k}": v for k, v in score_dict.items() if k not in ['warnings']}
                })
    
    # ==== STEP 4: Design primers ====
    primers_by_fragment = []
    
    for i, frag in enumerate(frags):
        pi = prev_index(i)
        ni = next_index(i)
        
        # Check junction status
        left_junction_key = (frags[pi].name, frag.name) if pi is not None else None
        right_junction_key = (frag.name, frags[ni].name) if ni is not None else None
        
        left_has_existing = left_junction_key in existing_overlaps
        right_has_existing = right_junction_key in existing_overlaps
        
        left_has_designed = left_junction_key in designed_overlaps
        right_has_designed = right_junction_key in designed_overlaps
        
        # Fragment at linear end
        left_is_end = pi is None
        right_is_end = ni is None
        
        # Determine if primers needed
        needs_left_primer = not (left_has_existing or left_is_end)
        needs_right_primer = not (right_has_existing or right_is_end)
        
        if not needs_left_primer and not needs_right_primer:
            primers_by_fragment.append({
                "fragment": frag.name,
                "needs_primers": False,
                "reason": "preexisting_homology_both_sides"
            })
            continue
        
        # Design annealing sequences
        f_anneal_seq, f_anneal_tm = design_annealing_sequence(
            template=frag.sequence[:int(req.anneal_max_len)],
            min_len=int(req.anneal_min_len),
            max_len=int(req.anneal_max_len),
            target_tm=float(req.anneal_tm_target)
        )

        rev_template = reverse_complement(frag.sequence)[:int(req.anneal_max_len)]
        r_anneal_seq, r_anneal_tm = design_annealing_sequence(
            template=rev_template,
            min_len=int(req.anneal_min_len),
            max_len=int(req.anneal_max_len),
            target_tm=float(req.anneal_tm_target)
        )

        if not f_anneal_seq or not r_anneal_seq:
            raise HTTPException(500, detail=f"Failed to design annealing for {frag.name}")

        # Design forward primer
        f_extension = ""
        f_extension_score = None
        
        if needs_left_primer and left_has_designed:
            # Get the LEFT portion of the junction (from previous fragment)
            overlap_seq, left_bp, right_bp, score_dict = designed_overlaps[left_junction_key]
            f_extension = overlap_seq[:left_bp]  # Portion from left (previous) fragment
            f_extension_score = score_dict
        
        forward_primer = f_extension + f_anneal_seq
        
        # Design reverse primer
        r_extension = ""
        r_extension_score = None
        
        if needs_right_primer and right_has_designed:
            # Get the RIGHT portion of the junction (from next fragment), then RC
            overlap_seq, left_bp, right_bp, score_dict = designed_overlaps[right_junction_key]
            right_portion = overlap_seq[left_bp:]  # Portion from right (next) fragment
            r_extension = reverse_complement(right_portion)  # CRITICAL: RC for reverse primer
            r_extension_score = score_dict
        
        reverse_primer = r_extension + r_anneal_seq
        
        # Calculate extension Tms
        thermo_calc = ThermodynamicCalculator()
        f_extension_tm = thermo_calc.calculate_tm(f_extension) if f_extension else None
        r_extension_tm = thermo_calc.calculate_tm(reverse_complement(r_extension) if r_extension else "") if r_extension else None
        
        primer_record = {
            "fragment": frag.name,
            "needs_primers": True,
            "forward_primer": forward_primer,
            "reverse_primer": reverse_primer,
            "forward_anneal_seq": f_anneal_seq,
            "reverse_anneal_seq": r_anneal_seq,
            "forward_anneal_tm": f_anneal_tm,
            "reverse_anneal_tm": r_anneal_tm,
            "forward_extension_seq": f_extension,
            "reverse_extension_seq": r_extension,  # This is RC of the homology
            "forward_extension_tm": f_extension_tm,
            "reverse_extension_tm": r_extension_tm,
            # Aliases for n8n viz node compatibility
            "forward_tail_seq": f_extension,
            "reverse_tail_seq": r_extension,
        }
        
        # Add scoring details
        if f_extension_score:
            for k, v in f_extension_score.items():
                if k != 'warnings':
                    primer_record[f"forward_extension_{k}"] = v
            primer_record["forward_extension_warnings"] = "; ".join(f_extension_score.get('warnings', []))
        
        if r_extension_score:
            for k, v in r_extension_score.items():
                if k != 'warnings':
                    primer_record[f"reverse_extension_{k}"] = v
            primer_record["reverse_extension_warnings"] = "; ".join(r_extension_score.get('warnings', []))
        
        primers_by_fragment.append(primer_record)
    
    # Build construct (simplified - just concatenate for now)
    # Build construct (collapse existing overlaps so we don't duplicate junction homology)
    assembled_seq, annotations = assemble_construct_with_existing_overlaps(
        frags=frags,
        assembly=req.assembly,
        existing_overlaps=existing_overlaps,
    )

    # --- Add overlap + primer binding annotations to viz/construct ---

    # fragment positions on assembled construct
    frag_pos = {a["name"]: a for a in annotations if "name" in a and "start" in a and "end" in a}

    extra_ann = []

    # 1) Overlap annotations
    for j in junctions:
        left_name = j["from"]
        right_name = j["to"]
        left = frag_pos.get(left_name)
        right = frag_pos.get(right_name)
        if not left or not right:
            continue

        src = j.get("source", "unknown")
        ov_len = int(j.get("overlap_length") or 0)
        left_bp = int(j.get("left_bp") or 0)
        right_bp = int(j.get("right_bp") or 0)

        if src == "existing" and ov_len > 0:
            # existing overlap exists ONCE in assembled construct:
            # it ends the LEFT fragment because the prefix of RIGHT was trimmed away
            s = max(left["start"], left["end"] - ov_len)
            e = left["end"]
            extra_ann.append({
                "name": f"{left_name}→{right_name} overlap ({ov_len}bp, existing)",
                "start": s,
                "end": e,
                "direction": 1,
                "type": "overlap",
                "source": "existing",
            })
        else:
            # designed overlap: single annotation spanning the junction
            ov_seq = j.get("overlap_sequence", "")
            span_len = len(ov_seq) if ov_seq else (left_bp + right_bp)
            s = max(left["start"], left["end"] - left_bp)
            e = min(right["end"], right["start"] + right_bp)
            extra_ann.append({
                "name": f"{left_name}→{right_name} overlap ({span_len}bp, designed)",
                "start": s,
                "end": e,
                "direction": 1,
                "type": "overlap",
                "source": "designed",
            })

    # 2) Primer binding (anneal) annotations - REMOVED
    # The n8n visualization node rebuilds primer annotations from primers_by_fragment data,
    # so we don't need to add anneal-only annotations here that would need filtering later.
    # This prevents duplicate/conflicting primer annotations in the final viz.

    # combine
    annotations = annotations + extra_ann
    construct = {
        "assembly": req.assembly,
        "sequence": assembled_seq,
        "annotations": annotations,
        "type": "gibson",
    }

    return {
        "sessionId": sid,
        "include_ai_explanation": req.include_ai_explanation,
        "assembly": req.assembly,
        "overlap_len": req.overlap_len,
        "homology_window": req.homology_window,
        "homology_min": req.homology_min,
        "junctions": junctions,
        "primers_by_fragment": primers_by_fragment,
        "primers": primers_by_fragment,
        "construct": construct,
        "viz": {
            "type": "gibson",
            "sequence": assembled_seq,
            "annotations": annotations,
        },
    }
