"""
SDMOperator — Q5-style site-directed mutagenesis primer designer.

Primer Design Strategies:
1. Deletions (any size): Back-to-back primers flanking the deletion
2. Insertions/Substitutions ≤6 bp: Single mutagenic primer with mismatch Tm calculation
3. Insertions/Substitutions >6 bp: Overlapping primers, both contain the mutation

Tm Calculations (PrimerX formula):
- Substitutions: Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch
- Insertions/Deletions: Tm = 81.5 + 0.41(%GC) - 675/N (no mismatch penalty)

Target annealing Tm: 61°C (Primer3 calculates ~7°C below NEB Q5 annealing temp of ~68°C)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Tuple

from ..utils import reverse_complement, normalize_dna
from .build_plan import BuildStep, OperatorMetrics, SDMBuildPlan, SDMPrimerDesign
from .lab_profile import DEFAULT_LAB_PROFILE, LabProfile

# Q5 targets 61°C Tm (Primer3 calculates ~7°C below NEB's recommended annealing temp)
_ANNEAL_TM_TARGET = 61.0
_ANNEAL_MIN = 18
_ANNEAL_MAX = 28

# Threshold for single-primer vs overlapping-primer strategy
_SINGLE_PRIMER_MAX_EDIT = 6

# SDM constraints
_MAX_INSERTION_BP = 41
_MAX_SUBSTITUTION_BP = 41
# No constraint on deletion size


def _primerx_tm(
    sequence: str,
    mismatch_count: int = 0,
    is_indel: bool = False,
) -> float:
    """
    Calculate Tm using PrimerX formula.
    
    For substitutions: includes mismatch penalty
    For insertions/deletions: no mismatch penalty, N excludes indel bases
    
    Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch
    """
    seq = normalize_dna(sequence)
    if not seq or len(seq) < 2:
        return 0.0
    
    n = len(seq)
    gc_count = seq.count("G") + seq.count("C")
    gc_percent = (gc_count / n) * 100
    
    tm = 81.5 + 0.41 * gc_percent - 675 / n
    
    # Apply mismatch penalty only for substitutions, not indels
    if not is_indel and mismatch_count > 0:
        mismatch_percent = (mismatch_count / n) * 100
        tm -= mismatch_percent
    
    return tm


def _primer3_tm(sequence: str) -> Optional[float]:
    """Calculate Tm using primer3 for comparison/validation."""
    try:
        import primer3
        seq = normalize_dna(sequence)
        if not seq:
            return None
        return float(primer3.calcTm(seq))
    except Exception:
        return None


def _gc_content(sequence: str) -> float:
    """Calculate GC content as fraction."""
    seq = normalize_dna(sequence)
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


@dataclass
class SDMPrimerResult:
    """Result of SDM primer design."""
    strategy: str  # "single_primer" | "overlapping" | "back_to_back"
    forward_primer: str
    reverse_primer: str
    forward_tm: float
    reverse_tm: float
    forward_anneal_seq: str
    reverse_anneal_seq: str
    overlap_tm: Optional[float] = None  # For overlapping strategy
    overlap_seq: Optional[str] = None
    warnings: List[str] = field(default_factory=list)


class SDMOperator:
    """Design Q5-style primers for site-directed mutagenesis."""

    def __init__(self, lab_profile: Optional[LabProfile] = None) -> None:
        self.lab = lab_profile or DEFAULT_LAB_PROFILE

    def evaluate(
        self,
        template_seq: str,
        old_seq: str,
        new_seq: str,
        template_name: str = "template_plasmid",
        features: Optional[List[Dict[str, Any]]] = None,
        insertion_position: Optional[int] = None,
    ) -> SDMBuildPlan:
        """
        Design SDM primers for the specified edit.
        
        Args:
            template_seq: Full plasmid sequence
            old_seq: Sequence to be replaced/deleted (empty for pure insertion)
            new_seq: Replacement sequence (empty for deletion)
            template_name: Name for the template plasmid
            features: Optional list of plasmid features for context
            insertion_position: Position for pure insertions (when old_seq is empty)
        """
        plan = SDMBuildPlan(template_name=template_name)

        t = normalize_dna(template_seq)
        old = normalize_dna(old_seq)
        new = normalize_dna(new_seq)

        # Validate inputs
        if not t:
            plan.feasible = False
            plan.infeasibility_reasons.append("template_seq is required")
            return plan
        
        # Determine mutation type and validate
        if not old and new:
            mutation_type = "insertion"
            if insertion_position is None:
                plan.feasible = False
                plan.infeasibility_reasons.append(
                    "Pure insertion requires insertion_position parameter"
                )
                return plan
            start = insertion_position
            edit_end = start
        elif old and not new:
            mutation_type = "deletion"
            start = t.find(old)
            if start < 0:
                plan.feasible = False
                plan.infeasibility_reasons.append("old_seq not found in template sequence")
                return plan
            edit_end = start + len(old)
        else:
            mutation_type = "substitution"
            start = t.find(old)
            if start < 0:
                plan.feasible = False
                plan.infeasibility_reasons.append("old_seq not found in template sequence")
                return plan
            edit_end = start + len(old)
            
            # Check for multiple occurrences
            if t.find(old, start + 1) >= 0:
                plan.warnings.append(
                    "old_seq appears multiple times in template; using first occurrence"
                )

        # Check SDM constraints
        if mutation_type == "insertion" and len(new) > _MAX_INSERTION_BP:
            plan.feasible = False
            plan.infeasibility_reasons.append(
                f"Insertion of {len(new)} bp exceeds {_MAX_INSERTION_BP} bp limit; "
                "use Gibson assembly instead"
            )
            return plan
        
        if mutation_type == "substitution" and len(new) > _MAX_SUBSTITUTION_BP:
            plan.feasible = False
            plan.infeasibility_reasons.append(
                f"Substitution with {len(new)} bp exceeds {_MAX_SUBSTITUTION_BP} bp limit; "
                "use Gibson assembly instead"
            )
            return plan

        # Design primers based on mutation type and size
        edit_size = max(len(old), len(new))
        
        if mutation_type == "deletion":
            # Deletions always use back-to-back primers
            primer_result = self._design_back_to_back_primers(t, start, edit_end, new)
        elif edit_size <= _SINGLE_PRIMER_MAX_EDIT:
            # Small insertions/substitutions: single mutagenic primer
            primer_result = self._design_single_primer(t, start, edit_end, old, new, mutation_type)
        else:
            # Larger insertions/substitutions: overlapping primers
            primer_result = self._design_overlapping_primers(t, start, edit_end, old, new)

        plan.warnings.extend(primer_result.warnings)
        
        # Build mutated sequence for visualization
        mutated_seq = t[:start] + new + t[edit_end:]

        plan.primer_design = SDMPrimerDesign(
            template_name=template_name,
            edit_start=start,
            edit_end=edit_end,
            old_sequence=old,
            new_sequence=new,
            forward_primer=primer_result.forward_primer,
            reverse_primer=primer_result.reverse_primer,
            forward_tm=round(primer_result.forward_tm, 1),
            reverse_tm=round(primer_result.reverse_tm, 1),
            warnings=primer_result.warnings,
        )
        
        # Store additional data for visualization and files
        plan.mutated_sequence = mutated_seq
        plan.mutation_type = mutation_type
        plan.primer_strategy = primer_result.strategy
        plan.fwd_anneal_seq = primer_result.forward_anneal_seq
        plan.rev_anneal_seq = primer_result.reverse_anneal_seq
        plan.overlap_tm = primer_result.overlap_tm
        plan.overlap_seq = primer_result.overlap_seq

        plan.bom = [
            "Q5 High-Fidelity 2X Master Mix",
            "KLD Enzyme Mix (kinase, ligase, DpnI)",
            "Competent cells (e.g., NEB 5-alpha)",
            "LB agar + selection antibiotic",
            "Sanger sequencing primers",
        ]

        plan.steps = [
            BuildStep(1, "pcr", 
                f"Q5 SDM PCR ({primer_result.strategy}): {mutation_type} {len(old)} bp -> {len(new)} bp",
                ["Q5 master mix", "SDM primers"], 1.5, 0.0),
            BuildStep(2, "kld", 
                "KLD treatment (5 min RT): phosphorylate, ligate, digest template",
                ["KLD enzyme mix"], 0.25, 0.0),
            BuildStep(3, "transformation", 
                "Transform 5 uL KLD product into competent cells",
                ["Competent cells", "LB + antibiotic plates"], 0.5, 1.0),
            BuildStep(4, "colony_screening", 
                "Pick 2-4 colonies, miniprep, sequence verify",
                ["Sequencing primers"], 1.5, 1.0),
        ]

        plan.metrics = self._compute_metrics(
            len(primer_result.forward_primer), 
            len(primer_result.reverse_primer)
        )
        
        tm_info = f"Fwd Tm={primer_result.forward_tm:.1f}C, Rev Tm={primer_result.reverse_tm:.1f}C"
        if primer_result.overlap_tm:
            tm_info += f", Overlap Tm={primer_result.overlap_tm:.1f}C"
            
        plan.summary = (
            f"Q5 SDM ({primer_result.strategy}) for {mutation_type}: "
            f"{len(old)} bp -> {len(new)} bp at position {start}; "
            f"{tm_info}; "
            f"estimated ${plan.metrics.total_cost_usd:.2f}, {plan.metrics.total_calendar_days:.1f} days"
        )

        return plan

    def _design_back_to_back_primers(
        self,
        template: str,
        start: int,
        edit_end: int,
        new_seq: str,
    ) -> SDMPrimerResult:
        """
        Design back-to-back primers for deletions.
        
        Forward primer: [new_seq (if any)][downstream annealing]
        Reverse primer: reverse_complement([upstream annealing])
        """
        downstream_seq = template[edit_end:]
        upstream_seq = template[:start]
        
        fwd_anneal = self._pick_anneal(downstream_seq, reverse_side=False)
        rev_anneal_template = self._pick_anneal(upstream_seq, reverse_side=True)
        
        fwd_primer = f"{new_seq}{fwd_anneal}"
        rev_primer = reverse_complement(rev_anneal_template)
        
        # Tm for annealing portions (no mismatch for deletions)
        fwd_tm = _primerx_tm(fwd_anneal, mismatch_count=0, is_indel=True)
        rev_tm = _primerx_tm(rev_anneal_template, mismatch_count=0, is_indel=True)
        
        warnings = []
        if abs(fwd_tm - rev_tm) > 4.0:
            warnings.append(
                f"Annealing Tm mismatch: forward={fwd_tm:.1f}C, reverse={rev_tm:.1f}C"
            )
        
        return SDMPrimerResult(
            strategy="back_to_back",
            forward_primer=fwd_primer,
            reverse_primer=rev_primer,
            forward_tm=fwd_tm,
            reverse_tm=rev_tm,
            forward_anneal_seq=fwd_anneal,
            reverse_anneal_seq=rev_anneal_template,
            warnings=warnings,
        )

    def _design_single_primer(
        self,
        template: str,
        start: int,
        edit_end: int,
        old_seq: str,
        new_seq: str,
        mutation_type: str,
    ) -> SDMPrimerResult:
        """
        Design single mutagenic primer for small (<=6 bp) insertions/substitutions.
        
        The mutagenic primer contains the mutation flanked by template-matching regions.
        The reverse primer is a simple template-matching primer on the opposite strand.
        
        Forward primer: [upstream flank][mutation][downstream flank]
        Reverse primer: reverse_complement of template region ~150-200 bp away
        """
        # Get flanking regions for the mutagenic primer
        upstream_flank_len = 12  # bp of template before mutation
        downstream_flank_len = 12  # bp of template after mutation
        
        upstream_start = max(0, start - upstream_flank_len)
        upstream_flank = template[upstream_start:start]
        
        downstream_end = min(len(template), edit_end + downstream_flank_len)
        downstream_flank = template[edit_end:downstream_end]
        
        # Build mutagenic primer
        fwd_primer = f"{upstream_flank}{new_seq}{downstream_flank}"
        
        # Calculate mismatch count for Tm
        if mutation_type == "substitution":
            mismatch_count = sum(1 for a, b in zip(old_seq, new_seq) if a != b)
            is_indel = False
        else:
            mismatch_count = 0
            is_indel = True
        
        # For the mutagenic primer, calculate Tm excluding the mutation
        # The annealing Tm considers only the matching portions
        anneal_seq = upstream_flank + downstream_flank
        fwd_tm = _primerx_tm(anneal_seq, mismatch_count=mismatch_count, is_indel=is_indel)
        
        # Design reverse primer ~150-200 bp downstream
        rev_start = min(edit_end + 150, len(template) - _ANNEAL_MIN)
        if rev_start < edit_end + 50:
            rev_start = edit_end + 50
        
        # Handle circular template
        if rev_start >= len(template):
            rev_start = 50  # Wrap around to start
        
        rev_anneal_template = self._pick_anneal(template[rev_start:], reverse_side=False)
        rev_primer = reverse_complement(rev_anneal_template)
        rev_tm = _primerx_tm(rev_anneal_template, mismatch_count=0, is_indel=True)
        
        warnings = []
        if abs(fwd_tm - rev_tm) > 5.0:
            warnings.append(
                f"Tm difference: mutagenic={fwd_tm:.1f}C, reverse={rev_tm:.1f}C"
            )
        
        return SDMPrimerResult(
            strategy="single_primer",
            forward_primer=fwd_primer,
            reverse_primer=rev_primer,
            forward_tm=fwd_tm,
            reverse_tm=rev_tm,
            forward_anneal_seq=anneal_seq,
            reverse_anneal_seq=rev_anneal_template,
            warnings=warnings,
        )

    def _design_overlapping_primers(
        self,
        template: str,
        start: int,
        edit_end: int,
        old_seq: str,
        new_seq: str,
    ) -> SDMPrimerResult:
        """
        Design overlapping primers for larger (>6 bp) insertions/substitutions.
        
        Both primers contain the mutation. They overlap at the mutation site,
        with the overlap Tm calculated for the region where they anneal to each other.
        
        For insertions: split the insertion between both primers
        Forward primer: [upstream_anneal][first_half_of_insertion]
        Reverse primer: RC([second_half_of_insertion][downstream_anneal])
        """
        # Get annealing regions
        upstream_seq = template[:start]
        downstream_seq = template[edit_end:]
        
        upstream_anneal = self._pick_anneal(upstream_seq, reverse_side=True)
        downstream_anneal = self._pick_anneal(downstream_seq, reverse_side=False)
        
        # Split the new sequence between primers
        if len(new_seq) > 0:
            split_point = len(new_seq) // 2
            fwd_mutation_part = new_seq[:split_point]
            rev_mutation_part = new_seq[split_point:]
        else:
            fwd_mutation_part = ""
            rev_mutation_part = ""
        
        # Build primers
        # Forward: upstream_anneal + first half of mutation
        fwd_primer = f"{upstream_anneal}{fwd_mutation_part}"
        
        # Reverse: RC(second half of mutation + downstream_anneal)
        rev_template_seq = f"{rev_mutation_part}{downstream_anneal}"
        rev_primer = reverse_complement(rev_template_seq)
        
        # Calculate Tms
        fwd_tm = _primerx_tm(upstream_anneal, mismatch_count=0, is_indel=True)
        rev_tm = _primerx_tm(downstream_anneal, mismatch_count=0, is_indel=True)
        
        # Calculate overlap Tm (where primers anneal to each other through the mutation)
        # The overlap is the full new_seq region
        overlap_tm = _primerx_tm(new_seq, mismatch_count=0, is_indel=True) if new_seq else 0.0
        
        warnings = []
        if abs(fwd_tm - rev_tm) > 4.0:
            warnings.append(
                f"Annealing Tm mismatch: forward={fwd_tm:.1f}C, reverse={rev_tm:.1f}C"
            )
        if overlap_tm < 45.0 and len(new_seq) > 0:
            warnings.append(
                f"Low overlap Tm ({overlap_tm:.1f}C); consider longer overlap"
            )
        
        return SDMPrimerResult(
            strategy="overlapping",
            forward_primer=fwd_primer,
            reverse_primer=rev_primer,
            forward_tm=fwd_tm,
            reverse_tm=rev_tm,
            forward_anneal_seq=upstream_anneal,
            reverse_anneal_seq=downstream_anneal,
            overlap_tm=overlap_tm,
            overlap_seq=new_seq,
            warnings=warnings,
        )

    def _pick_anneal(
        self,
        seq: str,
        reverse_side: bool = False,
        tm_target: float = _ANNEAL_TM_TARGET,
    ) -> str:
        """
        Pick optimal annealing region targeting specified Tm.
        
        Args:
            seq: Sequence to pick annealing region from
            reverse_side: If True, pick from end of seq (for upstream/reverse)
            tm_target: Target Tm (default 61C for Q5)
        """
        if not seq:
            return "N" * _ANNEAL_MIN

        candidates = []
        for ln in range(_ANNEAL_MIN, min(_ANNEAL_MAX + 1, len(seq) + 1)):
            if len(seq) < ln:
                break
            
            # Pick from appropriate end
            c = seq[-ln:] if reverse_side else seq[:ln]
            
            # Calculate Tm using PrimerX formula (no mismatch for annealing region)
            tm = _primerx_tm(c, mismatch_count=0, is_indel=True)
            
            candidates.append((abs(tm - tm_target), tm, c))

        if not candidates:
            # Fallback: use minimum length or available sequence
            if reverse_side and len(seq) >= _ANNEAL_MIN:
                return seq[-_ANNEAL_MIN:]
            elif len(seq) >= _ANNEAL_MIN:
                return seq[:_ANNEAL_MIN]
            return seq if seq else "N" * _ANNEAL_MIN

        # Sort by distance from target Tm
        candidates.sort(key=lambda x: x[0])
        return candidates[0][2]

    def _compute_metrics(self, fwd_len: int, rev_len: int) -> OperatorMetrics:
        """Compute cost/time metrics for SDM."""
        metrics = OperatorMetrics()
        metrics.primer_count = 2
        metrics.pcr_count = 1
        metrics.assembly_count = 1  # KLD
        metrics.transformation_count = 1
        metrics.sequencing_count = 2  # Sequence 2 colonies

        # SDM primers can be longer due to mutation payload — bill each
        # primer at the actual length (lab.primer_cost = $0.24/bp).
        metrics.primer_cost_usd = self.lab.primer_cost(fwd_len) + self.lab.primer_cost(rev_len)

        metrics.pcr_cost_usd = self.lab.pcr_rxn_cost_usd
        # Q5 SDM Kit covers the KLD step + cells + buffers.
        metrics.assembly_cost_usd = self.lab.catalog["q5_sdm_kit"].cost_per_rxn
        metrics.transformation_cost_usd = self.lab.transformation_cost_usd
        # ONT: 1 read /construct.
        metrics.sequencing_count = self.lab.sequencing_reads_per_construct
        metrics.sequencing_cost_usd = metrics.sequencing_count * self.lab.sequencing_cost_usd

        metrics.total_cost_usd = (
            metrics.primer_cost_usd
            + metrics.pcr_cost_usd
            + metrics.assembly_cost_usd
            + metrics.transformation_cost_usd
            + metrics.sequencing_cost_usd
        )

        metrics.total_labor_hours = 3.75
        metrics.total_calendar_days = 2.0
        metrics.pcr_risk = 0.05
        metrics.assembly_risk = 0.05
        metrics.overall_risk_score = 0.05
        
        return metrics
