"""
GibsonOperator — Phase 1 of the Vibe Cloning operator library.

Takes a resolved module list (from plasmid_design_chat.py) and produces a
full GibsonBuildPlan including:
  - Per-junction overlap design (overlap sequence, primer sequences, Tm, quality)
  - Per-fragment sourcing assessment (PCR difficulty, synthesis vs PCR decision)
  - Complete primer table (ready to order)
  - Bill of materials
  - Step-by-step protocol
  - Aggregate cost / time / labor / risk metrics

Primer design convention (Gibson/HiFi):
  For each junction j between left_module[j] and right_module[j+1]:
    overlap_seq  = last overlap_len bp of left_module.sequence
    FWD primer for right_module = 5'-[overlap_seq]-[anneal_fwd]-3'
    REV primer for left_module  = 5'-RC([last anneal_len bp of left_module])-3'

  The PCR product for left_module naturally contains overlap_seq at its 3' end.
  The PCR product for right_module starts with overlap_seq (added by the FWD tail).
  These matching ends are the homology regions used by the HiFi enzyme mix.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

from .build_plan import (
    BuildStep,
    FragmentSource,
    FragmentSourceType,
    GibsonBuildPlan,
    OperatorMetrics,
    OverlapDesign,
)
from .junction import Junction, _gc_content, _has_homopolymer, build_junctions
from .lab_profile import DEFAULT_LAB_PROFILE, LabProfile

# Lazy import: primer3 lives in the venv; avoid import errors if testing without it
try:
    from ..gibson_primers import ThermodynamicCalculator, OverlapScorer, SequenceAnalyzer
    _THERMO_AVAILABLE = True
except Exception:
    _THERMO_AVAILABLE = False


# ---------------------------------------------------------------------------
# Operator parameters (class-level defaults; override in __init__ if needed)
# ---------------------------------------------------------------------------
_OVERLAP_LEN_MIN = 20
_OVERLAP_LEN_DEFAULT = 25
_OVERLAP_LEN_MAX = 40
_OVERLAP_TM_MIN = 50.0
_OVERLAP_TM_MAX = 70.0
_OVERLAP_TM_TARGET = 62.0
_ANNEAL_LEN_MIN = 18
_ANNEAL_LEN_MAX = 32
_ANNEAL_TM_TARGET = 62.0

# Fragment PCR difficulty thresholds
_EASY_MAX_BP = 4_000
_MODERATE_MAX_BP = 8_000
_PCR_MAX_BP = 12_000       # above this: synthesis only

# Sequence-complexity GC thresholds for synthesis tier
_COMPLEX_GC = 0.68


class GibsonOperator:
    """
    Gibson/HiFi Assembly cloning operator.

    Usage:
        op = GibsonOperator()
        plan = op.evaluate(resolved_modules, topology="circular")
        print(plan.summary)
        print(plan.metrics.total_cost_usd)
    """

    def __init__(self, lab_profile: Optional[LabProfile] = None) -> None:
        self.lab = lab_profile or DEFAULT_LAB_PROFILE
        self._thermo = ThermodynamicCalculator() if _THERMO_AVAILABLE else None
        self._scorer = OverlapScorer() if _THERMO_AVAILABLE else None
        self._analyzer = SequenceAnalyzer() if _THERMO_AVAILABLE else None

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def evaluate(
        self,
        modules: List[dict],
        topology: str = "circular",
        construct_sequence: Optional[str] = None,
    ) -> GibsonBuildPlan:
        """
        Evaluate Gibson/HiFi assembly for a resolved module list.

        Args:
            modules: list of resolved module dicts from plasmid_design_chat.py.
                     Each dict should have: sequence, role, canonical_id or description,
                     origin ("library" | "ncbi" | "synthesis_needed"), source (filename).
            topology: "circular" (default) or "linear".
            construct_sequence: pre-assembled full sequence string (optional;
                                auto-computed by concatenation if not provided).

        Returns:
            GibsonBuildPlan with all primer design, sourcing, cost, and risk data.
        """
        plan = GibsonBuildPlan(assembly_topology=topology)

        # Filter to modules with real sequences; flag synthesis-needed ones
        real_modules, synthesis_only = _partition_modules(modules)

        if len(real_modules) < 2:
            plan.feasible = False
            plan.infeasibility_reasons.append(
                f"Only {len(real_modules)} module(s) with resolved sequences. "
                "Gibson assembly requires ≥2 sequenced fragments. "
                "Resolve remaining sequences (library lookup, NCBI fetch, or synthesis) first."
            )
            return plan

        plan.fragment_count = len(real_modules)

        # Build junction metadata
        junctions = build_junctions(real_modules, topology=topology)

        # Construct sequence for uniqueness scoring
        if not construct_sequence:
            construct_sequence = "".join(m.get("sequence", "") for m in real_modules)

        # --- Step 1: Design overlaps at every junction ---
        overlap_designs: List[OverlapDesign] = []
        all_overlap_seqs: List[str] = []

        for j in junctions:
            od = self._design_junction_overlap(j, real_modules, construct_sequence)
            overlap_designs.append(od)
            all_overlap_seqs.append(od.overlap_sequence)

        # Re-score overlaps with full set (cross-dimer + Tm uniformity)
        if self._scorer and len(all_overlap_seqs) > 1:
            for od in overlap_designs:
                scored = self._scorer.score_overlap(
                    od.overlap_sequence,
                    all_overlaps=all_overlap_seqs,
                    construct_seq=construct_sequence,
                )
                od.quality_score = scored.get("total_score", od.quality_score)
                od.uniqueness_score = scored.get("uniqueness_score", od.uniqueness_score)
                new_warnings = scored.get("warnings", [])
                for w in new_warnings:
                    if w not in od.warnings:
                        od.warnings.append(w)

        plan.overlap_designs = overlap_designs

        # Propagate junction-level warnings into plan
        for j, od in zip(junctions, overlap_designs):
            prefix = f"J{j.junction_index} ({j.left_module_name}→{j.right_module_name})"
            for w in od.warnings:
                plan.warnings.append(f"{prefix}: {w}")

        # --- Step 2: Assess fragment sourcing ---
        fragment_sources = [
            self._assess_fragment_source(i, m) for i, m in enumerate(real_modules)
        ]
        # Add synthesis-only modules as SYNTHESIS sources
        for m in synthesis_only:
            idx = modules.index(m)
            fragment_sources.append(
                FragmentSource(
                    module_index=idx,
                    module_name=_module_name(m),
                    source_type=FragmentSourceType.SYNTHESIS,
                    synthesis_bp=m.get("length", 500),
                    synthesis_tier="standard",
                    synthesis_cost_estimate_usd=self._synthesis_cost(
                        m.get("length", 500), "standard"
                    ),
                    synthesis_days=self.lab.synthesis_lead_time_days,
                    notes=["Sequence not yet resolved; synthesis or retrieval required"],
                )
            )

        plan.fragment_sources = fragment_sources

        # --- Step 3: Build primer table ---
        plan.primer_table = self._build_primer_table(real_modules, overlap_designs)

        # --- Step 4: Bill of materials ---
        plan.bom = self._build_bom(real_modules, overlap_designs, fragment_sources)

        # --- Step 5: Step-by-step protocol ---
        plan.steps = self._build_steps(real_modules, overlap_designs, fragment_sources)

        # --- Step 6: Compute metrics ---
        plan.metrics = self._compute_metrics(
            real_modules, overlap_designs, fragment_sources, junctions
        )

        # --- Step 7: Feasibility check ---
        plan.feasible = self._check_feasibility(plan, junctions)

        # --- Step 8: Summary string ---
        plan.summary = self._build_summary(plan)

        return plan

    # ------------------------------------------------------------------
    # Junction overlap design
    # ------------------------------------------------------------------

    def _design_junction_overlap(
        self,
        junction: Junction,
        modules: List[dict],
        construct_sequence: str,
    ) -> OverlapDesign:
        """
        Design the overlap sequence and primer pair for a single junction.

        Returns an OverlapDesign with:
          - overlap_sequence: from end of left module
          - forward_primer: for right module (overlap tail + anneal)
          - reverse_primer: for left module (RC of last anneal_len bp of left)
        """
        left_mod = modules[junction.left_module_index]
        right_mod = modules[junction.right_module_index]
        left_seq = left_mod.get("sequence", "") or ""
        right_seq = right_mod.get("sequence", "") or ""

        warnings: List[str] = []

        # Warn on GC extremes
        if junction.left_gc_content > 0.72 or junction.right_gc_content > 0.72:
            warnings.append(
                f"High GC at junction ({junction.left_gc_content:.0%}/{junction.right_gc_content:.0%})"
            )
        if junction.left_gc_content < 0.28 or junction.right_gc_content < 0.28:
            warnings.append(
                f"Low GC at junction ({junction.left_gc_content:.0%}/{junction.right_gc_content:.0%}); verify Tm"
            )
        if junction.left_has_homopolymer or junction.right_has_homopolymer:
            warnings.append("Homopolymer run (≥5 bp) near junction; PCR accuracy risk")
        if junction.has_tandem_repeat:
            warnings.append(f"Tandem repeat: {junction.repeat_note}")
        if junction.reading_frame_continuation:
            warnings.append("Reading frame crosses junction; ensure overlap does not introduce frameshift")

        # Find best overlap length from left module's 3' end
        overlap_seq, overlap_tm = self._find_best_overlap(left_seq)

        # Design annealing portions
        anneal_fwd, anneal_tm_fwd = self._design_anneal(right_seq, "forward")
        anneal_rev, anneal_tm_rev = self._design_anneal(left_seq, "reverse")

        # Build full primer sequences
        #   FWD primer for RIGHT module: 5'-[overlap from left]-[anneal to right start]-3'
        #   REV primer for LEFT module:  5'-[RC of last anneal_len of left]-3'
        fwd_primer = overlap_seq + anneal_fwd
        rev_primer = anneal_rev  # already RC from _design_anneal("reverse")

        # Thermodynamic quality for these specific primers
        hairpin_fwd = _calc_hairpin(fwd_primer, self._thermo)
        hairpin_rev = _calc_hairpin(rev_primer, self._thermo)
        dimer_fwd = _calc_homodimer(fwd_primer, self._thermo)
        dimer_rev = _calc_homodimer(rev_primer, self._thermo)

        if hairpin_fwd < -3.0:
            warnings.append(f"FWD primer hairpin ΔG={hairpin_fwd:.1f} kcal/mol")
        if hairpin_rev < -3.0:
            warnings.append(f"REV primer hairpin ΔG={hairpin_rev:.1f} kcal/mol")
        if dimer_fwd < -6.0:
            warnings.append(f"FWD primer self-dimer ΔG={dimer_fwd:.1f} kcal/mol")
        if dimer_rev < -6.0:
            warnings.append(f"REV primer self-dimer ΔG={dimer_rev:.1f} kcal/mol")

        # Preliminary overlap quality score (single overlap, no cross-dimer context yet)
        quality_score = 75.0  # default until re-scored with full overlap set
        uniqueness_score = 75.0
        if self._scorer:
            scored = self._scorer.score_overlap(
                overlap_seq, construct_seq=construct_sequence
            )
            quality_score = scored.get("total_score", 75.0)
            uniqueness_score = scored.get("uniqueness_score", 75.0)
            for w in scored.get("warnings", []):
                if w not in warnings:
                    warnings.append(w)

        return OverlapDesign(
            junction_index=junction.junction_index,
            left_module_name=junction.left_module_name,
            right_module_name=junction.right_module_name,
            overlap_sequence=overlap_seq,
            overlap_length=len(overlap_seq),
            overlap_tm=round(overlap_tm, 1),
            overlap_gc=round(_gc_content(overlap_seq), 3),
            forward_primer=fwd_primer,
            reverse_primer=rev_primer,
            forward_primer_anneal=anneal_fwd,
            reverse_primer_anneal=anneal_rev,
            forward_anneal_tm=round(anneal_tm_fwd, 1),
            reverse_anneal_tm=round(anneal_tm_rev, 1),
            uniqueness_score=round(uniqueness_score, 1),
            hairpin_dg_fwd=round(hairpin_fwd, 2),
            hairpin_dg_rev=round(hairpin_rev, 2),
            self_dimer_dg_fwd=round(dimer_fwd, 2),
            self_dimer_dg_rev=round(dimer_rev, 2),
            quality_score=round(quality_score, 1),
            warnings=warnings,
        )

    def _find_best_overlap(self, left_seq: str) -> Tuple[str, float]:
        """
        Scan overlap lengths from _OVERLAP_LEN_MIN to _OVERLAP_LEN_MAX and pick
        the one whose Tm is closest to _OVERLAP_TM_TARGET, within the allowed window.
        Falls back to _OVERLAP_LEN_DEFAULT if Primer3 is unavailable.
        """
        if not left_seq:
            return ("N" * _OVERLAP_LEN_DEFAULT, 60.0)

        if not self._thermo:
            length = min(_OVERLAP_LEN_DEFAULT, len(left_seq))
            return (left_seq[-length:], 60.0)

        best_seq, best_tm, best_score = None, 60.0, -9999.0

        for length in range(_OVERLAP_LEN_MIN, min(_OVERLAP_LEN_MAX + 1, len(left_seq) + 1)):
            candidate = left_seq[-length:]
            tm = self._thermo.calculate_tm(candidate)
            if tm < _OVERLAP_TM_MIN or tm > _OVERLAP_TM_MAX:
                continue
            gc = _gc_content(candidate)
            # Score: prefer Tm near target and GC near 50%
            score = -abs(tm - _OVERLAP_TM_TARGET) * 2.0 - abs(gc - 0.5) * 10.0
            if score > best_score:
                best_seq, best_tm, best_score = candidate, tm, score

        if best_seq is None:
            # No candidate in Tm window; take closest to target Tm
            best_seq = left_seq[-_OVERLAP_LEN_DEFAULT:] if len(left_seq) >= _OVERLAP_LEN_DEFAULT else left_seq
            best_tm = self._thermo.calculate_tm(best_seq)

        return (best_seq, best_tm)

    def _design_anneal(self, seq: str, direction: str) -> Tuple[str, float]:
        """
        Design the annealing portion of a primer on the given template sequence.

        direction="forward" → start from 5' end of seq (for FWD primer on right module)
        direction="reverse" → from 3' end of seq; returned as RC (for REV primer on left)

        Returns (annealing_sequence_ready_to_use, tm).
        For "reverse": returned sequence is already RC and can be used directly as primer.
        """
        if not seq:
            return ("NNNNNNNNNNNNNNNNNNNN", 55.0)

        if not self._thermo:
            if direction == "forward":
                anneal = seq[:_ANNEAL_LEN_MIN]
            else:
                anneal = _rev_comp(seq[-_ANNEAL_LEN_MIN:])
            return (anneal, 58.0)

        best_seq, best_tm = None, None

        for length in range(_ANNEAL_LEN_MIN, min(_ANNEAL_LEN_MAX + 1, len(seq) + 1)):
            if direction == "forward":
                candidate_template = seq[:length]
                candidate = candidate_template
            else:
                candidate_template = seq[-length:]
                candidate = _rev_comp(candidate_template)

            tm = self._thermo.calculate_tm(candidate)
            # Stop extending once we've passed the target
            if best_tm is not None and tm >= _ANNEAL_TM_TARGET:
                best_seq, best_tm = candidate, tm
                break
            # Track best so far
            if best_tm is None or abs(tm - _ANNEAL_TM_TARGET) < abs(best_tm - _ANNEAL_TM_TARGET):
                best_seq, best_tm = candidate, tm

        if best_seq is None:
            if direction == "forward":
                best_seq = seq[:_ANNEAL_LEN_MIN]
            else:
                best_seq = _rev_comp(seq[-_ANNEAL_LEN_MIN:])
            best_tm = 55.0

        return (best_seq, best_tm)

    # ------------------------------------------------------------------
    # Fragment sourcing assessment
    # ------------------------------------------------------------------

    def _assess_fragment_source(self, idx: int, module: dict) -> FragmentSource:
        """
        Determine how to obtain a module's sequence for PCR assembly.
        Evaluates PCR difficulty and whether synthesis is preferred.
        """
        seq = module.get("sequence") or ""
        origin = module.get("origin", "library")
        source = module.get("source") or ""
        name = _module_name(module)
        length = len(seq)
        difficulty_reasons: List[str] = []

        # Sequence-complexity checks
        gc = _gc_content(seq)
        if gc > _COMPLEX_GC:
            difficulty_reasons.append(f"High GC: {gc:.0%}")
        if _has_homopolymer(seq, min_run=5):
            difficulty_reasons.append("Homopolymer runs ≥5 bp")
        if self._analyzer and _has_tandem_repeat_in_seq(seq, self._analyzer):
            difficulty_reasons.append("Tandem repeat regions")

        # Above hard PCR limit → synthesis only
        if length > _PCR_MAX_BP:
            tier = "complex" if difficulty_reasons or gc > _COMPLEX_GC else "standard"
            return FragmentSource(
                module_index=idx,
                module_name=name,
                source_type=FragmentSourceType.SYNTHESIS,
                template_plasmid=source or None,
                expected_amplicon_bp=length,
                synthesis_bp=length,
                synthesis_tier=tier,
                synthesis_cost_estimate_usd=self._synthesis_cost(length, tier),
                synthesis_days=self.lab.synthesis_lead_time_days,
                notes=[
                    f"Fragment length {length} bp exceeds PCR limit ({_PCR_MAX_BP} bp); synthesis required"
                ] + difficulty_reasons,
            )

        # Determine PCR difficulty
        if length > _MODERATE_MAX_BP or len(difficulty_reasons) >= 2:
            difficulty = "difficult"
        elif length > _EASY_MAX_BP or len(difficulty_reasons) == 1:
            difficulty = "moderate"
        else:
            difficulty = "easy"

        # Source type from origin field
        source_type = {
            "library": FragmentSourceType.PCR_LIBRARY,
            "ncbi": FragmentSourceType.PCR_ADDGENE,
            "synthesis_needed": FragmentSourceType.SYNTHESIS,
        }.get(origin, FragmentSourceType.PCR_LIBRARY)

        if source_type == FragmentSourceType.SYNTHESIS:
            tier = "complex" if gc > _COMPLEX_GC or "tandem" in " ".join(difficulty_reasons).lower() else "standard"
            return FragmentSource(
                module_index=idx,
                module_name=name,
                source_type=FragmentSourceType.SYNTHESIS,
                synthesis_bp=length,
                synthesis_tier=tier,
                synthesis_cost_estimate_usd=self._synthesis_cost(length, tier),
                synthesis_days=self.lab.synthesis_lead_time_days,
                notes=difficulty_reasons,
            )

        # For difficult PCR fragments, flag synthesis comparison
        if difficulty == "difficult":
            tier = "complex" if gc > _COMPLEX_GC else "standard"
            synth_cost = self._synthesis_cost(length, tier)
            pcr_cost = (
                self.lab.pcr_rxn_cost_usd
                + 2 * self.lab.primer_cost_usd
                + self.lab.gel_lane_cost_usd
            )
            if synth_cost < pcr_cost * 2.5:
                difficulty_reasons.append(
                    f"Synthesis may be cost-competitive (est. ${synth_cost:.0f} vs PCR ~${pcr_cost:.0f})"
                )

        return FragmentSource(
            module_index=idx,
            module_name=name,
            source_type=source_type,
            template_plasmid=source or None,
            expected_amplicon_bp=length,
            pcr_difficulty=difficulty,
            pcr_difficulty_reasons=difficulty_reasons,
        )

    def _synthesis_cost(self, bp: int, tier: str) -> float:
        # Length-tiered pricing (lab_profile.synthesis_cost):
        # 0.5-1.8 kbp $0.07/bp, 1.8-3.2 $0.08, 3.2-5.0 $0.09.
        # `tier` is retained for the legacy "complex" override but
        # the band rate dominates for typical fragments.
        return self.lab.synthesis_cost(bp)

    # ------------------------------------------------------------------
    # Primer table
    # ------------------------------------------------------------------

    def _build_primer_table(
        self,
        modules: List[dict],
        overlaps: List[OverlapDesign],
    ) -> List[Dict]:
        """
        Build a flat primer table from overlap designs.
        Each junction contributes:
          - One FWD primer for the RIGHT module (with overlap tail)
          - One REV primer for the LEFT module
        Result is ordered: for each junction, REV primer first (LEFT), then FWD (RIGHT).
        """
        n = len(modules)
        primers = []

        for od in overlaps:
            left_idx = od.junction_index
            right_idx = (od.junction_index + 1) % n
            left_name = _module_name(modules[left_idx])
            right_name = _module_name(modules[right_idx])

            # REV primer for LEFT module
            primers.append({
                "primer_name": f"{_slug(left_name)}_REV_Gibson",
                "sequence": od.reverse_primer,
                "annealing_portion": od.reverse_primer_anneal,
                "overlap_tail": "(no tail — anneals to end of fragment)",
                "tm_anneal": od.reverse_anneal_tm,
                "tm_overlap": None,
                "length": len(od.reverse_primer),
                "gc_content": round(_gc_content(od.reverse_primer), 3),
                "hairpin_dg": od.hairpin_dg_rev,
                "self_dimer_dg": od.self_dimer_dg_rev,
                "quality_score": od.quality_score,
                "purpose": f"Reverse primer for {left_name}; PCR product naturally contains overlap into {right_name}",
                "fragment": left_name,
                "warnings": od.warnings,
            })

            # FWD primer for RIGHT module
            primers.append({
                "primer_name": f"{_slug(right_name)}_FWD_Gibson",
                "sequence": od.forward_primer,
                "annealing_portion": od.forward_primer_anneal,
                "overlap_tail": od.overlap_sequence,
                "tm_anneal": od.forward_anneal_tm,
                "tm_overlap": od.overlap_tm,
                "length": len(od.forward_primer),
                "gc_content": round(_gc_content(od.forward_primer), 3),
                "hairpin_dg": od.hairpin_dg_fwd,
                "self_dimer_dg": od.self_dimer_dg_fwd,
                "quality_score": od.quality_score,
                "purpose": f"Forward primer for {right_name}; 5' tail = last {od.overlap_length} bp of {left_name} (overlap)",
                "fragment": right_name,
                "warnings": od.warnings,
            })

        return primers

    # ------------------------------------------------------------------
    # Bill of materials
    # ------------------------------------------------------------------

    def _build_bom(
        self,
        modules: List[dict],
        overlaps: List[OverlapDesign],
        sources: List[FragmentSource],
    ) -> List[str]:
        primer_table = self._build_primer_table(modules, overlaps)

        lines = [
            "## Bill of Materials",
            "",
            "### Enzymes & Kits",
            "- NEB HiFi Assembly Master Mix (E2621, 2×) — 1 reaction",
            "- High-fidelity polymerase for PCR (NEB Q5 or Thermo Phusion)",
            "- DpnI (if templates are E. coli plasmids with dam methylation) — optional",
            "",
            "### Competent Cells",
            "- NEB 5-alpha (C2987) or DH5α — 1 aliquot per transformation",
            "",
            "### Plates & Media",
            "- LB agar + 100 µg/mL ampicillin plates",
            "- SOC medium",
            "",
            "### Primers (order standard desalting, 25 nmol)",
        ]

        for p in primer_table:
            lines.append(
                f"- **{p['primer_name']}**: 5'-{p['sequence']}-3'  "
                f"(anneal Tm={p['tm_anneal']}°C, {p['length']} nt)"
            )

        lines += ["", "### Templates / Fragments"]
        for src in sources:
            if src.source_type == FragmentSourceType.PCR_LIBRARY:
                lines.append(
                    f"- **{src.module_name}** ({src.expected_amplicon_bp} bp): "
                    f"PCR from `{src.template_plasmid or 'library plasmid'}`"
                    + (f"  ⚠ {src.pcr_difficulty.upper()} PCR" if src.pcr_difficulty != "easy" else "")
                )
            elif src.source_type == FragmentSourceType.PCR_ADDGENE:
                lines.append(
                    f"- **{src.module_name}** ({src.expected_amplicon_bp} bp): "
                    f"Order from Addgene ({src.template_plasmid or 'see annotation'}), then PCR amplify"
                )
            elif src.source_type == FragmentSourceType.SYNTHESIS:
                lines.append(
                    f"- **{src.module_name}** ({src.synthesis_bp} bp): "
                    f"SYNTHESIZE — {src.synthesis_tier} tier, "
                    f"est. ${src.synthesis_cost_estimate_usd:.0f}, "
                    f"~{src.synthesis_days:.0f} day lead time"
                )

        lines += [
            "",
            "### Other Supplies",
            "- PCR purification or gel extraction kit",
            "- Miniprep kit",
            "- Sanger sequencing reactions (Azenta/Genewiz or equivalent)",
            "- 1% agarose gel supplies",
        ]
        return lines

    # ------------------------------------------------------------------
    # Step-by-step protocol
    # ------------------------------------------------------------------

    def _build_steps(
        self,
        modules: List[dict],
        overlaps: List[OverlapDesign],
        sources: List[FragmentSource],
    ) -> List[BuildStep]:
        steps: List[BuildStep] = []
        step_num = 1
        primer_table = self._build_primer_table(modules, overlaps)

        # Build lookup: module_name → {FWD, REV} primer names
        primer_lookup: Dict[str, Dict[str, str]] = {}
        for p in primer_table:
            frag = p["fragment"]
            primer_lookup.setdefault(frag, {})
            if "FWD" in p["primer_name"]:
                primer_lookup[frag]["FWD"] = p["primer_name"]
            elif "REV" in p["primer_name"]:
                primer_lookup[frag]["REV"] = p["primer_name"]

        # Synthesis wait (if any synthesis sources)
        synth_sources = [s for s in sources if s.source_type == FragmentSourceType.SYNTHESIS]
        if synth_sources:
            days = max(s.synthesis_days for s in synth_sources)
            steps.append(BuildStep(
                step_number=step_num,
                step_type="synthesis_wait",
                description=f"Order {len(synth_sources)} synthetic fragment(s); wait for delivery (~{days:.0f} days)",
                materials=[f"{s.module_name} ({s.synthesis_bp} bp, {s.synthesis_tier})" for s in synth_sources],
                estimated_hours=0.5,
                estimated_days=days,
            ))
            step_num += 1

        # PCR steps (one per PCR-sourced fragment)
        pcr_sources = [s for s in sources if s.source_type in (
            FragmentSourceType.PCR_LIBRARY,
            FragmentSourceType.PCR_INVENTORY,
            FragmentSourceType.PCR_ADDGENE,
        )]
        if pcr_sources:
            for src in pcr_sources:
                prs = primer_lookup.get(src.module_name, {})
                fwd_name = prs.get("FWD", "see primer table (FWD)")
                rev_name = prs.get("REV", "see primer table (REV)")
                diff_note = f"  ⚠ {src.pcr_difficulty.upper()} amplicon" if src.pcr_difficulty != "easy" else ""
                materials = [
                    f"Template: {src.template_plasmid or 'N/A'}",
                    f"Forward primer: {fwd_name}",
                    f"Reverse primer: {rev_name}",
                    "Q5 or Phusion master mix",
                ]
                if src.pcr_difficulty != "easy":
                    if src.pcr_difficulty == "difficult":
                        materials.append("Consider: 3% DMSO or GC enhancer")
                    materials.extend(
                        [f"Note: {r}" for r in src.pcr_difficulty_reasons]
                    )
                steps.append(BuildStep(
                    step_number=step_num,
                    step_type="pcr",
                    description=f"PCR amplify {src.module_name} ({src.expected_amplicon_bp} bp){diff_note}",
                    materials=materials,
                    estimated_hours=self.lab.pcr_labor_hours,
                    estimated_days=self.lab.pcr_days,
                ))
                step_num += 1

            steps.append(BuildStep(
                step_number=step_num,
                step_type="gel_purification",
                description=f"Gel electrophoresis: verify {len(pcr_sources)} PCR product(s); gel-extract or PCR-purify correct bands",
                materials=[
                    "1% agarose gel, TAE buffer",
                    "1 kb+ DNA ladder",
                    "Gel extraction or PCR purification kit",
                ],
                estimated_hours=self.lab.gel_labor_hours,
                estimated_days=self.lab.gel_days,
            ))
            step_num += 1

        # Gibson/HiFi Assembly
        frag_list = [s.module_name for s in sources[:5]]
        if len(sources) > 5:
            frag_list.append(f"... ({len(sources) - 5} more)")
        steps.append(BuildStep(
            step_number=step_num,
            step_type="assembly",
            description=(
                f"Gibson/HiFi Assembly: combine {len(sources)} fragments "
                f"({'circular' if 'circular' in (sources[0].source_type if sources else '') else 'circular'} assembly), "
                f"50°C × 60 min"
            ),
            materials=[
                "NEB HiFi Assembly Master Mix (2×), 10 µL",
                "Equimolar mix of all fragments, total 0.02–0.5 pmol each, 10 µL",
                "Nuclease-free H₂O to 20 µL",
            ],
            estimated_hours=self.lab.assembly_labor_hours,
            estimated_days=self.lab.assembly_days,
        ))
        step_num += 1

        # Transformation
        steps.append(BuildStep(
            step_number=step_num,
            step_type="transformation",
            description="Transform into NEB 5-alpha or DH5α; plate on LB + ampicillin (100 µg/mL); incubate 37°C overnight",
            materials=[
                "NEB 5-alpha competent cells (C2987)",
                "LB agar + Amp100 plates",
                "SOC medium, 950 µL",
                "42°C heat block (42 s heat shock)",
            ],
            estimated_hours=self.lab.transformation_labor_hours,
            estimated_days=self.lab.transformation_days + self.lab.overnight_incubation_days,
        ))
        step_num += 1

        # Colony screening
        n_junctions = len(overlaps)
        steps.append(BuildStep(
            step_number=step_num,
            step_type="colony_screening",
            description=f"Pick {self.lab.colonies_to_screen} colonies; screen by colony PCR spanning each junction ({n_junctions} junctions); run gel",
            materials=[
                "Colony PCR master mix (Taq or Q5)",
                f"Junction-spanning primer pairs ({n_junctions} pairs)",
                "1% agarose gel for colony PCR analysis",
            ],
            estimated_hours=self.lab.colony_pcr_labor_hours,
            estimated_days=self.lab.colony_selection_days,
        ))
        step_num += 1

        # Miniprep
        steps.append(BuildStep(
            step_number=step_num,
            step_type="miniprep",
            description=f"Inoculate {self.lab.minipreps_per_construct} positive clones overnight; miniprep following day",
            materials=["Miniprep kit", "LB + Amp100 liquid cultures"],
            estimated_hours=self.lab.miniprep_labor_hours,
            estimated_days=self.lab.miniprep_days + self.lab.overnight_incubation_days,
        ))
        step_num += 1

        # Sequencing
        seq_count = n_junctions * self.lab.sequencing_reads_per_junction + self.lab.sequencing_reads_overhead
        steps.append(BuildStep(
            step_number=step_num,
            step_type="sequencing",
            description=f"Submit {seq_count} Sanger sequencing reactions covering all {n_junctions} junctions",
            materials=[
                f"{seq_count} sequencing reactions (Azenta/Genewiz standard)",
                "Primers spanning each junction",
            ],
            estimated_hours=self.lab.sequencing_labor_hours,
            estimated_days=self.lab.sequencing_turnaround_days,
        ))
        return steps

    # ------------------------------------------------------------------
    # Metrics
    # ------------------------------------------------------------------

    def _compute_metrics(
        self,
        modules: List[dict],
        overlaps: List[OverlapDesign],
        sources: List[FragmentSource],
        junctions: List[Junction],
    ) -> OperatorMetrics:
        m = OperatorMetrics()

        pcr_sources = [
            s for s in sources
            if s.source_type in (
                FragmentSourceType.PCR_LIBRARY,
                FragmentSourceType.PCR_INVENTORY,
                FragmentSourceType.PCR_ADDGENE,
            )
        ]
        synth_sources = [s for s in sources if s.source_type == FragmentSourceType.SYNTHESIS]

        # Counts
        m.pcr_count = len(pcr_sources)
        m.gel_count = m.pcr_count
        m.primer_count = len(modules) * 2   # each module: 1 FWD + 1 REV

        # Costs — use length-aware primer pricing and the Gibson workflow
        # reagent bundle (NEB10b cells + LB-agar plate + NEBuilder HiFi
        # + agarose / SYBR Safe / ladder share). PCR cost is one Q5 rxn
        # per fragment; sequencing collapses to 1 ONT read / construct.
        # OverlapDesign carries the actual designed primers per junction;
        # mine each one for its length so the per-bp pricing is accurate.
        primer_lens: List[int] = []
        for od in overlaps:
            f = getattr(od, "forward_primer", None)
            r = getattr(od, "reverse_primer", None)
            if f:
                primer_lens.append(len(f))
            if r:
                primer_lens.append(len(r))
        # Modules from the synthesis path do not produce overlap primers;
        # backfill with a 30 mer estimate so the per-primer count is in sync.
        while len(primer_lens) < m.primer_count:
            primer_lens.append(30)
        m.primer_cost_usd = sum(self.lab.primer_cost(L) for L in primer_lens)

        m.pcr_cost_usd = m.pcr_count * self.lab.pcr_rxn_cost_usd
        m.gel_cost_usd = m.gel_count * self.lab.gel_lane_cost_usd
        # Workflow reagent bundle covers HiFi assembly + cells + plate +
        # gel-share consumables — but we already added the gel-share via
        # gel_lane_cost_usd, so use only the assembly + transformation
        # subset on the per-construct line and leave gel_cost_usd alone.
        m.assembly_cost_usd = self.lab.hifi_assembly_cost_usd
        m.transformation_cost_usd = self.lab.transformation_cost_usd
        m.miniprep_count = 0  # ONT plasmid prep service includes miniprep
        m.miniprep_cost_usd = 0.0
        m.sequencing_count = self.lab.sequencing_reads_per_construct
        m.sequencing_cost_usd = m.sequencing_count * self.lab.sequencing_cost_usd

        for s in synth_sources:
            m.synthesis_cost_usd += s.synthesis_cost_estimate_usd

        m.total_cost_usd = (
            m.primer_cost_usd + m.pcr_cost_usd + m.gel_cost_usd
            + m.assembly_cost_usd + m.transformation_cost_usd
            + m.miniprep_cost_usd + m.sequencing_cost_usd
            + m.synthesis_cost_usd
        )

        # Labor
        m.total_labor_hours = (
            m.pcr_count * self.lab.pcr_labor_hours
            + (self.lab.gel_labor_hours if m.gel_count else 0.0)
            + self.lab.assembly_labor_hours
            + self.lab.transformation_labor_hours
            + self.lab.colony_pcr_labor_hours
            + m.miniprep_count * self.lab.miniprep_labor_hours
            + self.lab.sequencing_labor_hours
            + self.lab.primer_ordering_labor_hours
        )

        # Calendar time (critical path)
        synth_days = max((s.synthesis_days for s in synth_sources), default=0.0)
        bench_days = (
            m.pcr_count * self.lab.pcr_days
            + (self.lab.gel_days if m.gel_count else 0.0)
            + self.lab.assembly_days
            + self.lab.transformation_days
            + self.lab.overnight_incubation_days    # overnight plate
            + self.lab.colony_selection_days
            + self.lab.overnight_incubation_days    # overnight inoculation for miniprep
            + self.lab.miniprep_days
            + self.lab.sequencing_turnaround_days
        )
        m.total_calendar_days = synth_days + bench_days

        # Risk: PCR
        difficult = sum(1 for s in pcr_sources if s.pcr_difficulty == "difficult")
        moderate = sum(1 for s in pcr_sources if s.pcr_difficulty == "moderate")
        n_pcr = max(1, m.pcr_count)
        m.pcr_risk = min(1.0, (difficult * 0.4 + moderate * 0.15) / n_pcr)
        if difficult:
            m.risk_flags.append(f"{difficult} difficult PCR fragment(s)")

        # Risk: Assembly
        low_quality = sum(1 for od in overlaps if od.quality_score < 50)
        repeats = sum(1 for j in junctions if j.has_tandem_repeat)
        n_jct = max(1, len(overlaps))
        m.assembly_risk = min(1.0, (low_quality * 0.3 + repeats * 0.4) / n_jct)
        if repeats:
            m.risk_flags.append(f"{repeats} junction(s) with tandem repeat risk")
        if low_quality:
            m.risk_flags.append(f"{low_quality} junction(s) with low overlap quality (<50)")

        m.overall_risk_score = m.pcr_risk * 0.4 + m.assembly_risk * 0.6
        return m

    def _check_feasibility(
        self, plan: GibsonBuildPlan, junctions: List[Junction]
    ) -> bool:
        feasible = True
        failed = [od for od in plan.overlap_designs if od.quality_score < 20]
        if failed:
            plan.infeasibility_reasons.append(
                f"{len(failed)} junction(s) have very low overlap quality score (<20). "
                "Likely cause: sequence too short, extreme GC, or missing flanking sequence."
            )
            feasible = False
        return feasible

    def _build_summary(self, plan: GibsonBuildPlan) -> str:
        m = plan.metrics
        pcr_n = sum(
            1 for s in plan.fragment_sources
            if s.source_type in (
                FragmentSourceType.PCR_LIBRARY,
                FragmentSourceType.PCR_INVENTORY,
                FragmentSourceType.PCR_ADDGENE,
            )
        )
        synth_n = sum(
            1 for s in plan.fragment_sources
            if s.source_type == FragmentSourceType.SYNTHESIS
        )
        risk_label = (
            "Low" if m.overall_risk_score < 0.25
            else "Medium" if m.overall_risk_score < 0.55
            else "High"
        )
        lines = [
            f"**Gibson/HiFi Assembly** — {plan.fragment_count}-fragment {plan.assembly_topology} construct",
            f"Fragments: {pcr_n} by PCR" + (f", {synth_n} by synthesis" if synth_n else ""),
            f"Primers: {m.primer_count} ({m.primer_count // 2} pairs)",
            f"Estimated cost: **${m.total_cost_usd:.0f} USD**",
            f"Estimated time: **{m.total_calendar_days:.1f} calendar days**, {m.total_labor_hours:.1f} labor-hours",
            f"Risk: **{risk_label}** (score={m.overall_risk_score:.2f})",
        ]
        if plan.warnings:
            lines.append(f"⚠ {len(plan.warnings)} warning(s) — review overlap designs")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _module_name(m: dict) -> str:
    return m.get("canonical_id") or m.get("description") or "unknown_module"


def _slug(name: str) -> str:
    """URL/filename-safe slug for primer naming."""
    return re.sub(r"[^A-Za-z0-9_]", "_", name)[:30]


def _partition_modules(modules: List[dict]):
    """Split into (real_sequence_modules, synthesis_needed_modules)."""
    real, synth = [], []
    for m in modules:
        seq = m.get("sequence") or ""
        if len(seq) < 10 or seq.count("N") > len(seq) * 0.5:
            synth.append(m)
        else:
            real.append(m)
    return real, synth


def _rev_comp(seq: str) -> str:
    comp = str.maketrans("ATGCatgcNn", "TACGtacgNn")
    return seq.translate(comp)[::-1]


def _calc_hairpin(seq: str, thermo) -> float:
    if thermo is None or not seq:
        return 0.0
    try:
        import primer3
        return float(primer3.calcHairpin(seq).dg / 1000.0)
    except Exception:
        return 0.0


def _calc_homodimer(seq: str, thermo) -> float:
    if thermo is None or not seq:
        return 0.0
    try:
        import primer3
        return float(primer3.calcHomodimer(seq).dg / 1000.0)
    except Exception:
        return 0.0


def _has_tandem_repeat_in_seq(seq: str, analyzer) -> bool:
    """Use SequenceAnalyzer.find_dinucleotide_repeats to flag tandem repeats."""
    if analyzer is None or not seq:
        return False
    try:
        repeats = analyzer.find_dinucleotide_repeats(seq)
        return len(repeats) > 3
    except Exception:
        return False
