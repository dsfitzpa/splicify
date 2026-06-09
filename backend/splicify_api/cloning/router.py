"""
FastAPI router: /cloning/*

Exposes cloning operator evaluation as REST endpoints.

POST /cloning/evaluate/gibson
    Evaluate Gibson/HiFi assembly for a resolved module list.
    Returns a complete GibsonBuildPlan (primers, sourcing, cost, risk).
"""
from __future__ import annotations

import base64
from pathlib import Path
from .. import _data
import re
from typing import Any, Dict, List, Literal, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from ..utils import reverse_complement
from .diff_router import DiffRouter
from .re_database import RE_DATABASE
from .gibson_operator import GibsonOperator
from .golden_gate_operator import GoldenGateOperator
from .lab_profile import LabProfile
from .restriction_operator import RestrictionOperator
from .sdm_operator import SDMOperator
from .sgrna_oligo_designer import (
    assemble_sgrna_plasmid,
    design_sgrna_oligos,
    design_sgrna_oligos_from_vector,
    find_type_iis_sites,
    load_lenticrispr_v2,
    parse_genbank_sequence,
    parse_genbank_features,
    TYPE_IIS_ENZYMES,
)
from .golden_gate_primer_designer import (
    design_multi_fragment_assembly,
    design_single_fragment_replacement,
    design_scarless_deletion,
    build_design_response,
    GoldenGatePrimerDesign,
)

router = APIRouter(prefix="/cloning", tags=["cloning"])

_DEMO_PUC19_PATH = (
    _data.data_path("Module_Library_gb")
    / "Basic Cloning Vectors"
    / "pUC19.gb"
)
_DEMO_EGFP_SEQUENCE = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG"
)
_DEMO_FORWARD_TAIL = "GGTGGTAAGCTT"
_DEMO_REVERSE_TAIL = "GGTGGTGAATTC"


# ---------------------------------------------------------------------------
# Request / Response models
# ---------------------------------------------------------------------------

class ResolvedModuleInput(BaseModel):
    """
    Mirrors a resolved module dict from plasmid_design_chat.py.
    Accepts both the existing chat.py format and a direct JSON submission.
    """
    canonical_id: Optional[str] = None
    description: Optional[str] = None
    role: str = "other"
    sequence: Optional[str] = None
    length: Optional[int] = None
    origin: str = "library"       # "library" | "ncbi" | "synthesis_needed"
    source: Optional[str] = None  # filename / Addgene accession
    # Extra fields ignored (allows passing full resolved module dicts from the frontend)

    class Config:
        extra = "allow"


class LabProfileInput(BaseModel):
    """Partial LabProfile override. Any field omitted uses the default."""
    primer_cost_usd: Optional[float] = None
    pcr_rxn_cost_usd: Optional[float] = None
    hifi_assembly_cost_usd: Optional[float] = None
    transformation_cost_usd: Optional[float] = None
    miniprep_cost_usd: Optional[float] = None
    sequencing_cost_usd: Optional[float] = None
    synthesis_cost_per_bp_simple: Optional[float] = None
    synthesis_cost_per_bp_complex: Optional[float] = None
    synthesis_threshold_bp: Optional[int] = None
    synthesis_lead_time_days: Optional[float] = None
    colonies_to_screen: Optional[int] = None


class GibsonEvaluateRequest(BaseModel):
    """
    Request body for POST /cloning/evaluate/gibson.

    Provide either:
      - `modules`: list of resolved module dicts (from plasmid_design_chat.py session)
      - `session_id`: retrieve resolved modules from an active design session (TODO)
    """
    modules: List[Dict[str, Any]] = Field(
        ...,
        description=(
            "Resolved module list from plasmid_design_chat.py. "
            "Each item needs at minimum: sequence, role, and canonical_id or description."
        ),
    )
    topology: str = Field("circular", description="'circular' or 'linear'")
    lab_profile: Optional[LabProfileInput] = Field(
        None, description="Optional cost overrides; omit to use default academic lab pricing"
    )


class RestrictionEvaluateRequest(BaseModel):
    modules: List[Dict[str, Any]] = Field(
        ...,
        description=(
            "Resolved module list. Each item needs sequence and role; "
            "canonical_id/description used for reporting."
        ),
    )
    topology: str = Field("circular", description="'circular' or 'linear'")
    lab_profile: Optional[LabProfileInput] = Field(
        None, description="Optional cost overrides"
    )


class DiffAnalyzeRequest(BaseModel):
    desired_modules: List[Dict[str, Any]] = Field(
        ...,
        description="Target module list to build.",
    )
    anchor_modules: Optional[List[Dict[str, Any]]] = Field(
        None, description="Optional explicit anchor plasmid module list."
    )
    lab_profile: Optional[LabProfileInput] = Field(
        None, description="Optional cost overrides used when SDM plan is produced"
    )


class GoldenGateEvaluateRequest(BaseModel):
    modules: List[Dict[str, Any]] = Field(
        ...,
        description="Resolved module list for Golden Gate evaluation.",
    )
    topology: str = Field("circular", description="'circular' or 'linear'")
    lab_profile: Optional[LabProfileInput] = Field(
        None, description="Optional cost overrides"
    )


class CompareRequest(BaseModel):
    modules: List[Dict[str, Any]] = Field(
        ...,
        description="Resolved module list to evaluate across operators.",
    )
    topology: str = Field("circular", description="'circular' or 'linear'")
    objective: Literal["balanced", "cheapest", "fastest", "fewest_steps", "lowest_risk"] = "balanced"
    run_diff_first: bool = True
    anchor_modules: Optional[List[Dict[str, Any]]] = None
    lab_profile: Optional[LabProfileInput] = None


class SgRNAOligoRequest(BaseModel):
    """Request model for custom sgRNA oligo design."""
    grna_sequence: str = Field(
        ...,
        description="20-25 bp gRNA target sequence (without PAM)",
        min_length=15,
        max_length=30,
    )
    vector_sequence: Optional[str] = Field(
        None,
        description="Optional accepting vector sequence (GenBank or raw DNA). If not provided, uses lentiCRISPR v2.",
    )
    enzyme: str = Field(
        "auto",
        description="Type IIS enzyme: 'BsmBI', 'BbsI', 'BsaI', or 'auto' to detect from vector",
    )


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@router.post("/evaluate/gibson", summary="Evaluate Gibson/HiFi assembly build plan")
async def evaluate_gibson(req: GibsonEvaluateRequest) -> Dict[str, Any]:
    """
    Given a list of resolved modules (from a plasmid design), evaluate Gibson/HiFi
    assembly as a construction method.

    Returns a complete GibsonBuildPlan including:
    - Per-junction overlap design (primer sequences, Tm, quality scores)
    - Per-fragment sourcing assessment (PCR difficulty, synthesis vs PCR)
    - Complete primer table (ready to order)
    - Bill of materials
    - Step-by-step protocol
    - Aggregate cost, time, labor, and risk metrics
    """
    if not req.modules:
        raise HTTPException(status_code=400, detail="modules list is empty")

    # Build LabProfile (merge overrides onto defaults)
    lab = _build_lab_profile(req.lab_profile)

    op = GibsonOperator(lab_profile=lab)

    try:
        plan = op.evaluate(
            modules=req.modules,
            topology=req.topology,
        )
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Gibson operator error: {exc}") from exc

    return {
        "ok": True,
        "method": plan.method,
        "feasible": plan.feasible,
        "infeasibility_reasons": plan.infeasibility_reasons,
        "summary": plan.summary,
        "warnings": plan.warnings,
        "fragment_count": plan.fragment_count,
        "topology": plan.assembly_topology,
        "metrics": {
            "primer_count": plan.metrics.primer_count,
            "pcr_count": plan.metrics.pcr_count,
            "total_cost_usd": round(plan.metrics.total_cost_usd, 2),
            "total_calendar_days": round(plan.metrics.total_calendar_days, 1),
            "total_labor_hours": round(plan.metrics.total_labor_hours, 1),
            "overall_risk_score": round(plan.metrics.overall_risk_score, 3),
            "risk_flags": plan.metrics.risk_flags,
            # Full cost breakdown
            "cost_breakdown": {
                "primers_usd": round(plan.metrics.primer_cost_usd, 2),
                "pcr_usd": round(plan.metrics.pcr_cost_usd, 2),
                "gel_usd": round(plan.metrics.gel_cost_usd, 2),
                "assembly_usd": round(plan.metrics.assembly_cost_usd, 2),
                "transformation_usd": round(plan.metrics.transformation_cost_usd, 2),
                "miniprep_usd": round(plan.metrics.miniprep_cost_usd, 2),
                "sequencing_usd": round(plan.metrics.sequencing_cost_usd, 2),
                "synthesis_usd": round(plan.metrics.synthesis_cost_usd, 2),
            },
        },
        "primer_table": plan.primer_table,
        "fragment_sources": [
            {
                "module_index": s.module_index,
                "module_name": s.module_name,
                "source_type": s.source_type,
                "template_plasmid": s.template_plasmid,
                "expected_amplicon_bp": s.expected_amplicon_bp,
                "pcr_difficulty": s.pcr_difficulty,
                "pcr_difficulty_reasons": s.pcr_difficulty_reasons,
                "synthesis_bp": s.synthesis_bp,
                "synthesis_tier": s.synthesis_tier,
                "synthesis_cost_estimate_usd": round(s.synthesis_cost_estimate_usd, 2),
                "synthesis_days": s.synthesis_days,
                "notes": s.notes,
            }
            for s in plan.fragment_sources
        ],
        "overlap_designs": [
            {
                "junction_index": od.junction_index,
                "left_module": od.left_module_name,
                "right_module": od.right_module_name,
                "overlap_sequence": od.overlap_sequence,
                "overlap_length": od.overlap_length,
                "overlap_tm": od.overlap_tm,
                "overlap_gc": od.overlap_gc,
                "forward_primer": od.forward_primer,
                "reverse_primer": od.reverse_primer,
                "forward_anneal_tm": od.forward_anneal_tm,
                "reverse_anneal_tm": od.reverse_anneal_tm,
                "uniqueness_score": od.uniqueness_score,
                "quality_score": od.quality_score,
                "warnings": od.warnings,
            }
            for od in plan.overlap_designs
        ],
        "bom": plan.bom,
        "steps": [
            {
                "step_number": s.step_number,
                "step_type": s.step_type,
                "description": s.description,
                "materials": s.materials,
                "estimated_hours": s.estimated_hours,
                "estimated_days": s.estimated_days,
            }
            for s in plan.steps
        ],
    }


@router.post("/evaluate/restriction", summary="Evaluate restriction/ligation cloning plan")
async def evaluate_restriction(req: RestrictionEvaluateRequest) -> Dict[str, Any]:
    if not req.modules:
        raise HTTPException(status_code=400, detail="modules list is empty")

    lab = _build_lab_profile(req.lab_profile)
    op = RestrictionOperator(lab_profile=lab)

    try:
        plan = op.evaluate(modules=req.modules, topology=req.topology)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Restriction operator error: {exc}") from exc

    return {
        "ok": True,
        "method": plan.method,
        "feasible": plan.feasible,
        "infeasibility_reasons": plan.infeasibility_reasons,
        "summary": plan.summary,
        "warnings": plan.warnings,
        "fragment_count": plan.fragment_count,
        "topology": plan.assembly_topology,
        "metrics": {
            "primer_count": plan.metrics.primer_count,
            "pcr_count": plan.metrics.pcr_count,
            "total_cost_usd": round(plan.metrics.total_cost_usd, 2),
            "total_calendar_days": round(plan.metrics.total_calendar_days, 1),
            "total_labor_hours": round(plan.metrics.total_labor_hours, 1),
            "overall_risk_score": round(plan.metrics.overall_risk_score, 3),
            "risk_flags": plan.metrics.risk_flags,
        },
        "junction_plans": [
            {
                "junction_index": jp.junction_index,
                "left_module_name": jp.left_module_name,
                "right_module_name": jp.right_module_name,
                "left_enzyme": jp.left_enzyme,
                "right_enzyme": jp.right_enzyme,
                "strategy": jp.strategy,
                "scar_sequence": jp.scar_sequence,
                "internal_conflicts": jp.internal_conflicts,
                "warnings": jp.warnings,
            }
            for jp in plan.junction_plans
        ],
        "engineered_primer_table": plan.engineered_primer_table,
        "bom": plan.bom,
        "steps": [
            {
                "step_number": s.step_number,
                "step_type": s.step_type,
                "description": s.description,
                "materials": s.materials,
                "estimated_hours": s.estimated_hours,
                "estimated_days": s.estimated_days,
            }
            for s in plan.steps
        ],
    }


@router.post("/diff/analyze", summary="Analyze design diffs and recommend cloning strategy")
async def analyze_diffs(req: DiffAnalyzeRequest) -> Dict[str, Any]:
    if not req.desired_modules:
        raise HTTPException(status_code=400, detail="desired_modules list is empty")

    diff = DiffRouter().analyze(
        desired_modules=req.desired_modules,
        anchor_modules=req.anchor_modules,
    )

    response: Dict[str, Any] = {
        "ok": True,
        "anchor_source": diff.get("anchor_source"),
        "anchor_coverage": diff.get("anchor_coverage"),
        "edit_type": diff.get("edit_type"),
        "routing_recommendation": diff.get("routing_recommendation"),
        "module_diffs": diff.get("module_diffs", []),
    }

    # If this is an SDM-class edit, attach a full SDM build plan candidate.
    routing = str(diff.get("routing_recommendation") or "")
    if routing == "sdm_q5":
        sdm_payload = _build_sdm_candidate(
            req.desired_modules,
            req.anchor_modules,
            diff,
            _build_lab_profile(req.lab_profile),
        )
        if sdm_payload:
            response["sdm_plan"] = sdm_payload

    return response


@router.post("/evaluate/golden_gate", summary="Evaluate Golden Gate assembly plan")
async def evaluate_golden_gate(req: GoldenGateEvaluateRequest) -> Dict[str, Any]:
    if not req.modules:
        raise HTTPException(status_code=400, detail="modules list is empty")

    lab = _build_lab_profile(req.lab_profile)
    op = GoldenGateOperator(lab_profile=lab)

    try:
        plan = op.evaluate(modules=req.modules, topology=req.topology)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Golden Gate operator error: {exc}") from exc

    return {
        "ok": True,
        "method": plan.method,
        "feasible": plan.feasible,
        "infeasibility_reasons": plan.infeasibility_reasons,
        "summary": plan.summary,
        "warnings": plan.warnings,
        "fragment_count": plan.fragment_count,
        "topology": plan.assembly_topology,
        "enzyme": plan.enzyme,
        "strategy": plan.strategy,
        "domestication_burden": plan.domestication_burden,
        "metrics": {
            "primer_count": plan.metrics.primer_count,
            "pcr_count": plan.metrics.pcr_count,
            "total_cost_usd": round(plan.metrics.total_cost_usd, 2),
            "total_calendar_days": round(plan.metrics.total_calendar_days, 1),
            "total_labor_hours": round(plan.metrics.total_labor_hours, 1),
            "overall_risk_score": round(plan.metrics.overall_risk_score, 3),
            "risk_flags": plan.metrics.risk_flags,
        },
        "junction_plans": [
            {
                "junction_index": jp.junction_index,
                "left_module_name": jp.left_module_name,
                "right_module_name": jp.right_module_name,
                "overhang_4bp": jp.overhang_4bp,
                "enzyme": jp.enzyme,
                "strategy": jp.strategy,
                "overhang_fidelity_score": jp.overhang_fidelity_score,
                "warnings": jp.warnings,
            }
            for jp in plan.junction_plans
        ],
        "primer_table": plan.primer_table,
        "bom": plan.bom,
        "steps": [
            {
                "step_number": s.step_number,
                "step_type": s.step_type,
                "description": s.description,
                "materials": s.materials,
                "estimated_hours": s.estimated_hours,
                "estimated_days": s.estimated_days,
            }
            for s in plan.steps
        ],
    }


@router.post("/compare", summary="Run diff-first routing and rank cloning options")
async def compare_methods(req: CompareRequest) -> Dict[str, Any]:
    if not req.modules:
        raise HTTPException(status_code=400, detail="modules list is empty")

    lab = _build_lab_profile(req.lab_profile)
    candidates: List[Dict[str, Any]] = []
    diff_result: Optional[Dict[str, Any]] = None

    if req.run_diff_first:
        diff_result = DiffRouter().analyze(
            desired_modules=req.modules,
            anchor_modules=req.anchor_modules,
        )
        routing = str(diff_result.get("routing_recommendation") or "")
        if routing == "sdm_q5":
            sdm_payload = _build_sdm_candidate(req.modules, req.anchor_modules, diff_result, lab)
            if sdm_payload:
                candidates.append(_candidate_from_plan("q5_sdm", sdm_payload))

    # Always evaluate assembly backends for comparison unless modules are empty.
    gibson = GibsonOperator(lab_profile=lab).evaluate(modules=req.modules, topology=req.topology).to_dict()
    restriction = RestrictionOperator(lab_profile=lab).evaluate(modules=req.modules, topology=req.topology).to_dict()
    golden_gate = GoldenGateOperator(lab_profile=lab).evaluate(modules=req.modules, topology=req.topology).to_dict()

    candidates.extend([
        _candidate_from_plan("gibson_hifi", gibson),
        _candidate_from_plan("restriction_cloning", restriction),
        _candidate_from_plan("golden_gate", golden_gate),
    ])

    ranked = _rank_candidates(candidates, req.objective)

    return {
        "ok": True,
        "objective": req.objective,
        "diff_result": diff_result,
        "recommended_method": ranked[0]["method"] if ranked else None,
        "ranked_options": ranked,
    }


@router.get("/health", summary="Cloning operator health check")
async def cloning_health() -> Dict[str, Any]:
    try:
        from ..gibson_primers import ThermodynamicCalculator
        tc = ThermodynamicCalculator()
        tm = tc.calculate_tm("ATGCATGCATGCATGCATGC")
        primer3_ok = tm > 0
    except Exception:
        primer3_ok = False

    return {
        "ok": True,
        "primer3_available": primer3_ok,
        "operators": ["gibson_hifi", "restriction_cloning", "golden_gate", "diff_router", "q5_sdm", "sgrna_oligo"],
        "operators_todo": [],
    }


# ---------------------------------------------------------------------------
# sgRNA Golden Gate Oligo Design Endpoints
# ---------------------------------------------------------------------------

# Demo gRNA: EMX1 gene target (commonly used in CRISPR literature)
_DEMO_GRNA_SEQUENCE = "GAGTCCGAGCAGAAGAAGAA"
_DEMO_GRNA_NAME = "EMX1"


@router.post("/design/sgrna_oligo", summary="Design sgRNA cloning oligos")
async def design_sgrna_oligo_endpoint(req: SgRNAOligoRequest) -> Dict[str, Any]:
    """
    Design oligos for cloning a custom sgRNA into a Type IIS vector.

    Provide a gRNA sequence (20-25 bp, without PAM) and optionally:
    - A vector sequence (defaults to lentiCRISPR v2)
    - An enzyme choice (defaults to auto-detect)

    Returns oligo sequences with Tm values and annealing protocol.
    """
    grna = req.grna_sequence.strip().upper()

    # Validate gRNA contains only valid bases
    invalid_chars = set(grna) - set("ACGTN")
    if invalid_chars:
        raise HTTPException(
            status_code=400,
            detail=f"gRNA contains invalid characters: {invalid_chars}. Only A, C, G, T allowed.",
        )

    # Use provided vector or default to lentiCRISPR v2
    if req.vector_sequence:
        vector_seq = req.vector_sequence
        vector_name = "custom vector"
    else:
        try:
            vector_seq = load_lenticrispr_v2()
            vector_name = "lentiCRISPR v2"
        except FileNotFoundError as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to load default vector: {exc}",
            ) from exc

    # Design oligos
    if req.enzyme == "auto" or req.vector_sequence:
        # Auto-detect from vector
        oligo_design = design_sgrna_oligos_from_vector(
            vector_sequence=vector_seq,
            grna_sequence=grna,
            enzyme=req.enzyme,
        )
    else:
        # Use specified enzyme with standard vector overhangs
        # (the oligo overhangs will be the reverse complement)
        if req.enzyme == "BsmBI":
            five_prime = "GGTG"  # Vector overhang → oligo gets CACC
            three_prime = "GTTT"  # Vector overhang → oligo gets AAAC
        elif req.enzyme == "BbsI":
            five_prime = "GTGG"  # Standard for pX330/pX335
            three_prime = "GTTT"
        elif req.enzyme == "BsaI":
            five_prime = "CATT"  # Varies by vector → oligo gets AATG
            three_prime = "CGCT"  # → oligo gets AGCG
        else:
            raise HTTPException(
                status_code=400,
                detail=f"Unknown enzyme: {req.enzyme}. Supported: BsmBI, BbsI, BsaI, auto",
            )

        oligo_design = design_sgrna_oligos(
            grna_sequence=grna,
            five_prime_overhang=five_prime,
            three_prime_overhang=three_prime,
            enzyme=req.enzyme,
        )

    # Build response
    oligo_table = f"""| Name | Sequence (5' to 3') | Length | Tm |
|------|---------------------|--------|-----|
| sgRNA_Fwd | {oligo_design.forward_oligo} | {len(oligo_design.forward_oligo)} bp | {oligo_design.forward_tm}°C |
| sgRNA_Rev | {oligo_design.reverse_oligo} | {len(oligo_design.reverse_oligo)} bp | {oligo_design.reverse_tm}°C |"""

    reply = f"""## sgRNA Oligo Design

**Target gRNA**: {oligo_design.grna_sequence}
**Vector**: {vector_name}
**Enzyme**: {oligo_design.enzyme}

### Oligos to Order

{oligo_table}

### Annealed Insert Structure

```
{oligo_design.annealed_product_display}
```

### Quick Protocol

1. Resuspend oligos to 100 µM
2. Mix 1 µL each + 8 µL annealing buffer
3. Heat 95°C 5 min, cool to 25°C slowly
4. Dilute 1:200 for ligation
5. Ligate into {oligo_design.enzyme}-digested vector
"""

    if oligo_design.warnings:
        reply += "\n### Warnings\n\n" + "\n".join(f"- {w}" for w in oligo_design.warnings)

    # CSV file
    oligo_csv = f"""Name,Sequence,Length,Tm
sgRNA_Fwd,{oligo_design.forward_oligo},{len(oligo_design.forward_oligo)},{oligo_design.forward_tm}
sgRNA_Rev,{oligo_design.reverse_oligo},{len(oligo_design.reverse_oligo)},{oligo_design.reverse_tm}
"""

    files = [
        {
            "fileName": "sgrna_oligos.csv",
            "mimeType": "text/csv",
            "dataBase64": base64.b64encode(oligo_csv.encode("utf-8")).decode("ascii"),
        },
    ]

    return {
        "ok": True,
        "reply": reply,
        "files": files,
        "oligo_design": {
            "forward_oligo": oligo_design.forward_oligo,
            "reverse_oligo": oligo_design.reverse_oligo,
            "forward_tm": oligo_design.forward_tm,
            "reverse_tm": oligo_design.reverse_tm,
            "grna_sequence": oligo_design.grna_sequence,
            "enzyme": oligo_design.enzyme,
            "five_prime_overhang": oligo_design.five_prime_overhang,
            "three_prime_overhang": oligo_design.three_prime_overhang,
            "insert_length": oligo_design.insert_length,
            "warnings": oligo_design.warnings,
        },
    }


# ---------------------------------------------------------------------------
# Golden Gate Multi-Fragment Assembly Demo
# ---------------------------------------------------------------------------

# Demo sequences for 3-fragment assembly: EF1a promoter + eGFP + bGH polyA
_DEMO_EF1A_PROMOTER = (
    "GGCTCCGGTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCCGAGAAGTTGGGGGGAGGGGTCGGCAATTGAAC"
    "CGGTGCCTAGAGAAGGTGGCGCGGGGTAAACTGGGAAAGTGATGTCGTGTACTGGCTCCGCCTTTTTCCCGAGGGTGGGG"
    "GAGAACCGTATATAAGTGCAGTAGTCGCCGTGAACGTTCTTTTTCGCAACGGGTTTGCCGCCAGAACACAGGTAAGTGCC"
    "GTGTGTGGTTCCCGCGGGCCTGGCCTCTTTACGGGTTATGGCCCTTGCGTGCCTTGAATTACTTCCACCTGGCTGCAGTA"
    "CGTGATTCTTGATCCCGAGCTTCGGGTTGGAAGTGGGTGGGAGAGTTCGAGGCCTTGCGCTTAAGGAGCCCCTTCGCCTC"
    "GTGCTTGAGTTGAGGCCTGGCCTGGGCGCTGGGGCCGCCGCGTGCGAATCTGGTGGCACCTTCGCGCCTGTCTCGCTGCT"
    "TTCGATAAGTCTCTAGCCATTTAAAATTTTTGATGACCTGCTGCGACGCTTTTTTTCTGGCAAGATAGTCTTGTAAATGC"
    "GGGCCAAGATCTGCACACTGGTATTTCGGTTTTTGGGGCCGCGGGCGGCGACGGGGCCCGTGCGTCCCAGCGCACATGTT"
    "CGGCGAGGCGGGGCCTGCGAGCGCGGCCACCGAGAATCGGACGGGGGTAGTCTCAAGCTGGCCGGCCTGCTCTGGTGCCT"
    "GGCCTCGCGCCGCCGTGTATCGCCCCGCCCTGGGCGGCAAGGCTGGCCCGGTCGGCACCAGTTGCGTGAGCGGAAAGATG"
    "GCCGCTTCCCGGCCCTGCTGCAGGGAGCTCAAAATGGAGGACGCGGCGCTCGGGAGAGCGGGCGGGTGAGTCACCCACAC"
    "AAAGGAAAAGGGCCTTTCCGTCCTCAGCCGTCGCTTCATGTGACTCCACGGAGTACCGGGCGCCGTCCAGGCACCTCGAT"
    "TAGTTCTCGAGCTTTTGGAGTACGTCGTCTTTAGGTTGGGGGGAGGGGTTTTATGCGATGGAGTTTCCCCACACTGAGTG"
    "GGTGGAGACTGAAGTTAGGCCAGCTTGGCACTTGATGTAATTCTCCTTGGAATTTGCCCTTTTTGAGTTTGGATCTTGGT"
    "TCATTCTCAAGCCTCAGACAGTGGTTCAAAGTTTTTTTCTTCCATTTCAGGTGTCGTGA"
)

_DEMO_BGH_POLYA = (
    "CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACT"
    "GTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCA"
    "GGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCATGCTGGGGATGCGGTGGGCTCTATGG"
)


@router.get("/demo/gateway_bp_real", summary="Gateway BP Real demo with actual plasmids")
async def gateway_bp_real_demo() -> Dict[str, Any]:
    """
    Demo: Clone insert from pCMV SPORT6 into pDONR221 using Gateway BP recombination.

    Uses real plasmid files from Module_Library_gb:
    - pCMV SPORT6.gb (donor plasmid with ccdB cassette)
    - pDONR221.gb (destination vector with attP sites)

    Returns:
        Complete Gateway build plan with primer sequences and visualization.
    """
    from .gateway_operator import GatewayOperator
    from Bio import SeqIO
    import os

    try:
        # Load real plasmid files
        pcmv_sport6_path = _data.data_path("Module_Library_gb", "Mammalian Expression Vectors", "pCMV SPORT6.gb")
        pdonr221_path = _data.data_path("Module_Library_gb", "Gateway Cloning Vectors", "pDONR221.gb")

        # Parse GenBank files
        pcmv_sport6_record = SeqIO.read(pcmv_sport6_path, "genbank")
        pdonr221_record = SeqIO.read(pdonr221_path, "genbank")

        # Extract sequences
        pcmv_sport6_seq = str(pcmv_sport6_record.seq)
        pdonr221_seq = str(pdonr221_record.seq)

        # Create modules for Gateway operator
        modules = [
            {
                "canonical_id": "pCMV_SPORT6",
                "sequence": pcmv_sport6_seq,
                "role": "insert",
                "description": f"pCMV SPORT6 donor plasmid ({len(pcmv_sport6_seq)} bp)"
            },
            {
                "canonical_id": "pDONR221",
                "sequence": pdonr221_seq,
                "role": "vector",
                "description": f"pDONR221 Gateway destination vector ({len(pdonr221_seq)} bp)"
            }
        ]

        # Run Gateway operator
        operator = GatewayOperator()
        plan = operator.evaluate(modules, topology="circular")

        if not plan.feasible:
            raise HTTPException(
                status_code=500,
                detail=f"Gateway design failed: {'; '.join(plan.infeasibility_reasons)}"
            )

        # Build response
        reply_parts = []

        # Header
        reply_parts.append("## Gateway BP Recombination: pCMV SPORT6 → pDONR221\n")
        reply_parts.append("This demo clones the insert from pCMV SPORT6 into pDONR221 using Gateway BP recombination.\n\n")

        # Plasmid info
        reply_parts.append("### Input Plasmids\n\n")
        reply_parts.append(f"- **pCMV SPORT6**: {len(pcmv_sport6_seq)} bp donor plasmid\n")
        reply_parts.append(f"- **pDONR221**: {len(pdonr221_seq)} bp destination vector with attP sites\n\n")

        # Reaction type
        reply_parts.append(f"**Reaction Type:** {plan.reaction_type}\n")
        reply_parts.append(f"**Topology:** {plan.assembly_topology}\n")
        reply_parts.append(f"**Fragments:** {plan.fragment_count}\n\n")

        # Junction information
        reply_parts.append("### Junctions\n\n")
        for jp in plan.junction_plans:
            reply_parts.append(f"**Junction {jp.junction_index}:** {jp.left_module_name} → {jp.right_module_name}\n")
            reply_parts.append(f"- Strategy: {jp.strategy}\n")
            reply_parts.append(f"- att sites: {jp.left_att_site} + {jp.right_att_site} → {jp.product_left_site} + {jp.product_right_site}\n")
            if jp.warnings:
                reply_parts.append(f"- Warnings: {', '.join(jp.warnings)}\n")
            reply_parts.append("\n")

        # Primer table
        if plan.primer_table:
            reply_parts.append("### PCR Primers\n\n")
            reply_parts.append("| Primer Name | Length | Tm | Purpose |\n")
            reply_parts.append("|-------------|--------|----|---------|\n")
            for primer in plan.primer_table:
                name = primer.get("primer_name", "Unknown")
                length = primer.get("length", 0)
                tm = primer.get("tm_anneal", 0)
                purpose = primer.get("purpose", "")
                reply_parts.append(f"| {name} | {length} bp | {tm:.1f}°C | {purpose} |\n")
            reply_parts.append("\n")

        # Metrics
        m = plan.metrics
        reply_parts.append("### Project Metrics\n\n")
        reply_parts.append(f"- **Cost:** ${m.total_cost_usd:.2f}\n")
        reply_parts.append(f"- **Time:** {m.total_calendar_days:.1f} days\n")
        reply_parts.append(f"- **Risk Score:** {m.overall_risk_score:.2f} (0=low, 1=high)\n")
        reply_parts.append(f"- **Primers needed:** {m.primer_count}\n")
        reply_parts.append(f"- **PCR reactions:** {m.pcr_count}\n\n")

        # Protocol steps
        reply_parts.append("### Protocol Steps\n\n")
        for step in plan.steps:
            reply_parts.append(f"{step.step_number}. **{step.description}**\n")
            if step.estimated_days > 0:
                reply_parts.append(f"   - Time: {step.estimated_days:.1f} days\n")
        reply_parts.append("\n")

        # Summary
        reply_parts.append("### Summary\n\n")
        reply_parts.append(plan.summary)

        # Build visualization data
        viz_data = {
            "type": "design",
            "sequence": plan.product_sequence,
            "annotations": plan.product_annotations,
            "topology": plan.assembly_topology,
            "title": "Gateway BP Product (Entry Clone)",
            "total_length": len(plan.product_sequence)
        }

        # Return response
        return {
            "ok": True,
            "reply": "".join(reply_parts),
            "viz": viz_data,
            "files": []
        }

    except FileNotFoundError as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Demo plasmid files not found: {exc}"
        ) from exc
    except Exception as exc:
        import traceback
        traceback.print_exc()
        raise HTTPException(
            status_code=500,
            detail=f"Gateway BP Real demo failed: {exc}"
        ) from exc


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_lab_profile(overrides: Optional[LabProfileInput]) -> LabProfile:
    lab = LabProfile()
    if overrides is None:
        return lab
    for field_name, value in overrides.model_dump(exclude_none=True).items():
        if hasattr(lab, field_name) and value is not None:
            setattr(lab, field_name, value)
    return lab


def _load_demo_puc19_record() -> Dict[str, Any]:
    text = _DEMO_PUC19_PATH.read_text(encoding="utf-8")
    sequence = ""
    features: List[Dict[str, Any]] = []
    in_origin = False
    in_features = False
    current: Optional[Dict[str, Any]] = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN"):
            in_origin = True
            in_features = False
            if current:
                features.append(current)
                current = None
            continue
        if line.startswith("//"):
            if current:
                features.append(current)
            break

        if in_origin:
            m = re.match(r"^\s*\d+\s+(.+)$", line)
            if m:
                sequence += re.sub(r"\s+", "", m.group(1)).upper()
            continue

        if in_features:
            feature_match = re.match(r"^     (\S+)\s+(.+)$", line)
            if feature_match and not line.startswith("                     /"):
                if current:
                    features.append(current)
                current = {
                    "type": feature_match.group(1),
                    "location": feature_match.group(2).strip(),
                    "qualifiers": {},
                }
                continue

            qualifier_match = re.match(r'^                     /(\w+)=?"?(.*?)(?:"?)?$', line)
            if qualifier_match and current:
                key = qualifier_match.group(1)
                value = qualifier_match.group(2).strip().strip('"')
                current["qualifiers"][key] = value

    return {"sequence": sequence, "features": features}


def _pick_best_anneal(template: str, reverse: bool, target_tm: float = 60.0) -> Tuple[str, float]:
    from ..gibson_primers import ThermodynamicCalculator

    thermo = ThermodynamicCalculator()
    seq = reverse_complement(template[-30:]) if reverse else template[:30]
    best_seq = seq[:18]
    best_tm = thermo.calculate_tm(best_seq)
    best_delta = abs(best_tm - target_tm)

    for length in range(18, min(30, len(seq)) + 1):
        candidate = seq[:length]
        tm = thermo.calculate_tm(candidate)
        delta = abs(tm - target_tm)
        if delta < best_delta:
            best_seq = candidate
            best_tm = tm
            best_delta = delta
    return best_seq, best_tm


def _remap_demo_features(
    features: List[Dict[str, Any]],
    hindiii_start_1b: int,
    replaced_end_1b: int,
    delta: int,
) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for feature in features:
        if feature.get("type") == "source":
            continue
        label = (
            feature.get("qualifiers", {}).get("label")
            or feature.get("qualifiers", {}).get("gene")
            or feature.get("qualifiers", {}).get("product")
            or feature.get("type")
            or "feature"
        )
        if str(label).lower().endswith("site"):
            continue

        ranges = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", feature.get("location", ""))]
        if not ranges:
            continue
        strand = -1 if "complement" in feature.get("location", "") else 1
        mapped_segments: List[Tuple[int, int]] = []
        for start, end in ranges:
            if end < hindiii_start_1b:
                mapped_segments.append((start, end))
            elif start > replaced_end_1b:
                mapped_segments.append((start + delta, end + delta))
            else:
                if start < hindiii_start_1b:
                    mapped_segments.append((start, hindiii_start_1b - 1))
                if end > replaced_end_1b:
                    mapped_segments.append((replaced_end_1b + delta + 1, end + delta))

        for seg_start, seg_end in mapped_segments:
            if seg_end <= seg_start:
                continue
            out.append({
                "name": str(label),
                "start": seg_start - 1,
                "end": seg_end,
                "direction": strand,
                "origin": "genbank",
                "source": "pUC19 reference",
                "color": "#86a8d9",
            })
    return out


async def _merge_plannotate_annotations(viz: Dict[str, Any]) -> Dict[str, Any]:
    seq = viz.get("sequence") or ""
    if not seq:
        return viz

    try:
        from ..plannotate_router import AnnotateSequenceRequest, annotate_sequence_endpoint
    except Exception:
        return viz

    response = await annotate_sequence_endpoint(
        AnnotateSequenceRequest(sequence=seq, circular=(viz.get("topology") != "linear"))
    )
    if not isinstance(response, dict) or not response.get("ok"):
        return viz

    merged = list(viz.get("annotations") or [])
    seen = {
        (int(a.get("start", -1)), int(a.get("end", -1)), str(a.get("name") or "").strip().lower())
        for a in merged
    }
    added: List[Dict[str, Any]] = []
    for ann in response.get("annotations") or []:
        key = (int(ann.get("start", -1)), int(ann.get("end", -1)), str(ann.get("name") or "").strip().lower())
        if key in seen:
            continue
        pl_ann = {
            "name": ann.get("name", "pLannotate feature"),
            "start": ann.get("start", 0),
            "end": ann.get("end", 0),
            "direction": ann.get("direction", 1),
            "color": "#A78BFA",
            "origin": "plannotate",
            "source": "pLannotate",
        }
        merged.append(pl_ann)
        added.append(pl_ann)
    viz["annotations"] = merged
    viz["plannotate_annotations"] = added
    viz["plannotate_annotation_count"] = len(added)
    return viz


def _dedupe_annotations(annotations: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    deduped: List[Dict[str, Any]] = []
    seen: Set[Tuple[int, int, str]] = set()
    for ann in annotations:
        key = (int(ann.get("start", -1)), int(ann.get("end", -1)), str(ann.get("name") or "").strip().lower())
        if key in seen:
            continue
        seen.add(key)
        deduped.append(ann)
    return deduped


def _build_demo_genbank(sequence: str, annotations: List[Dict[str, Any]], title: str) -> str:
    features = ["FEATURES             Location/Qualifiers", "     source          1..%d" % len(sequence)]
    features.append('                     /organism="synthetic DNA construct"')
    features.append('                     /mol_type="other DNA"')

    for ann in annotations:
        start = int(ann.get("start", 0)) + 1
        end = int(ann.get("end", 0))
        if end <= start:
            continue
        location = f"{start}..{end}"
        if int(ann.get("direction", 1)) < 0:
            location = f"complement({location})"
        features.append(f"     misc_feature    {location}")
        features.append(f'                     /label="{str(ann.get("name") or "feature").replace(chr(34), chr(39))}"')

    origin = ["ORIGIN"]
    seq_lower = sequence.lower()
    for i in range(0, len(seq_lower), 60):
        chunk = seq_lower[i:i + 60]
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        origin.append(f"{str(i + 1).rjust(9)} {groups}")

    return "\n".join([
        f"LOCUS       {title[:16].ljust(16)} {str(len(sequence)).rjust(6)} bp    DNA     circular SYN",
        "DEFINITION  pUC19 with GFP cloned into the MCS using HindIII and EcoRI.",
        "ACCESSION   DEMO0001",
        "VERSION     DEMO0001.2",
        "KEYWORDS    .",
        "SOURCE      synthetic DNA construct",
        "  ORGANISM  synthetic DNA construct",
        *features,
        *origin,
        "//",
    ])


def _build_sdm_candidate(
    desired_modules: List[Dict[str, Any]],
    anchor_modules: Optional[List[Dict[str, Any]]],
    diff: Dict[str, Any],
    lab: LabProfile,
) -> Optional[Dict[str, Any]]:
    rows = diff.get("module_diffs", [])
    modified_rows = [r for r in rows if _status_value(r.get("status")) == "MODIFIED"]
    if len(modified_rows) != 1:
        return None

    row = modified_rows[0]
    idx = int(row.get("module_index", -1))
    if idx < 0 or idx >= len(desired_modules):
        return None

    desired_mod = desired_modules[idx]
    key = _module_key(desired_mod)
    source = diff.get("anchor_source")
    resolved_anchor_modules = (
        anchor_modules if anchor_modules is not None else [m for m in desired_modules if m.get("source") == source]
    )
    anchor_lookup = {_module_key(m): m for m in resolved_anchor_modules}
    anchor_mod = anchor_lookup.get(key)
    if anchor_mod is None:
        return None

    edit = row.get("edit") or {}
    old_seq = (edit.get("old_seq") or "").upper()
    new_seq = (edit.get("new_seq") or "").upper()
    template_seq = (anchor_mod.get("sequence") or "").upper()
    if not template_seq or old_seq is None or new_seq is None:
        return None
    if old_seq == new_seq:
        return None

    template_name = (
        anchor_mod.get("canonical_id")
        or anchor_mod.get("description")
        or anchor_mod.get("source")
        or "anchor_template"
    )
    plan = SDMOperator(lab_profile=lab).evaluate(
        template_seq=template_seq,
        old_seq=old_seq,
        new_seq=new_seq,
        template_name=template_name,
    )
    return plan.to_dict()


def _status_value(status: Any) -> str:
    if hasattr(status, "value"):
        return str(status.value)
    return str(status)


def _module_key(module: Dict[str, Any]) -> str:
    return (
        module.get("canonical_id")
        or module.get("description")
        or module.get("role")
        or ""
    ).strip().lower()


def _candidate_from_plan(method: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    metrics = payload.get("metrics") or {}
    steps = payload.get("steps") or []
    return {
        "method": method,
        "feasible": bool(payload.get("feasible", True)),
        "summary": payload.get("summary", ""),
        "warnings": payload.get("warnings", []),
        "infeasibility_reasons": payload.get("infeasibility_reasons", []),
        "metrics": {
            "total_cost_usd": float(metrics.get("total_cost_usd", 1e9)),
            "total_calendar_days": float(metrics.get("total_calendar_days", 1e9)),
            "total_labor_hours": float(metrics.get("total_labor_hours", 1e9)),
            "overall_risk_score": float(metrics.get("overall_risk_score", 1.0)),
            "step_count": len(steps),
        },
        "plan": payload,
    }


def _rank_candidates(
    candidates: List[Dict[str, Any]],
    objective: str,
) -> List[Dict[str, Any]]:
    feasible = [c for c in candidates if c.get("feasible")]
    infeasible = [c for c in candidates if not c.get("feasible")]
    if not feasible:
        return candidates

    weights = _objective_weights(objective)
    metrics_keys = ["total_cost_usd", "total_calendar_days", "step_count", "overall_risk_score"]
    ranges: Dict[str, tuple[float, float]] = {}
    for key in metrics_keys:
        vals = [float(c["metrics"][key]) for c in feasible]
        ranges[key] = (min(vals), max(vals))

    ranked: List[Dict[str, Any]] = []
    for c in feasible:
        score = 0.0
        for key in metrics_keys:
            lo, hi = ranges[key]
            val = float(c["metrics"][key])
            norm = 0.0 if hi == lo else (val - lo) / (hi - lo)
            score += norm * weights[key]
        row = dict(c)
        row["weighted_score"] = round(score, 4)
        row["reason"] = _build_reason_line(row, objective)
        ranked.append(row)

    ranked.sort(key=lambda x: x["weighted_score"])
    return ranked + infeasible


def _objective_weights(objective: str) -> Dict[str, float]:
    if objective == "cheapest":
        return {"total_cost_usd": 0.6, "total_calendar_days": 0.15, "step_count": 0.1, "overall_risk_score": 0.15}
    if objective == "fastest":
        return {"total_cost_usd": 0.15, "total_calendar_days": 0.6, "step_count": 0.1, "overall_risk_score": 0.15}
    if objective == "fewest_steps":
        return {"total_cost_usd": 0.15, "total_calendar_days": 0.15, "step_count": 0.6, "overall_risk_score": 0.1}
    if objective == "lowest_risk":
        return {"total_cost_usd": 0.15, "total_calendar_days": 0.15, "step_count": 0.1, "overall_risk_score": 0.6}
    return {"total_cost_usd": 0.3, "total_calendar_days": 0.25, "step_count": 0.2, "overall_risk_score": 0.25}


def _build_reason_line(candidate: Dict[str, Any], objective: str) -> str:
    m = candidate["metrics"]
    return (
        f"Objective={objective}: cost=${m['total_cost_usd']:.2f}, "
        f"time={m['total_calendar_days']:.1f}d, steps={m['step_count']}, "
        f"risk={m['overall_risk_score']:.3f}"
    )
# Gateway Cloning Demo Endpoint
# Add this to ~/python-libraries/aiplasmiddesign_api/backend/splicify_api/cloning/router.py


# ---------------------------------------------------------------------------
# Cloning-feature scanner (direct) — applies Step 2.75 annotations without
# running the full LLM + pLannotate + rule pipeline. Used by the Plasmid
# Visualizer to apply cloning-feature annotations to an uploaded plasmid.
# ---------------------------------------------------------------------------

class ScanCloningFeaturesRequest(BaseModel):
    sequence: str = Field(..., description="DNA sequence to scan")
    enabled_sets: Optional[List[str]] = Field(
        default=None,
        description=(
            "Subset of {restriction_II, restriction_IIs, gateway, pcr}. "
            "Default: all four."
        ),
    )
    type_ii_enzymes: Optional[List[str]] = None
    type_iis_enzymes: Optional[List[str]] = None
    gateway_fuzzy_threshold: int = 0


@router.post(
    "/scan_cloning_features",
    summary="Scan cloning features (Type II/IIs sites, Gateway att, PCR warnings) — skips LLM/pLannotate",
)
async def scan_cloning_features_endpoint(
    req: ScanCloningFeaturesRequest,
) -> Dict[str, Any]:
    sequence = (req.sequence or "").strip()
    if not sequence:
        raise HTTPException(status_code=400, detail="sequence is required")

    from ..cloning_feature_annotator import (
        scan_cloning_features,
        cloning_features_to_hierarchical,
    )

    scan = scan_cloning_features(
        sequence,
        enabled_sets=req.enabled_sets,
        type_ii_enzymes=req.type_ii_enzymes,
        type_iis_enzymes=req.type_iis_enzymes,
        gateway_fuzzy_threshold=req.gateway_fuzzy_threshold,
    )
    hierarchical = cloning_features_to_hierarchical(scan.features)

    return {
        "ok": True,
        "hierarchical_annotations": hierarchical,
        "cloning_features": scan.to_dict(),
        "summary": {
            "cloning_feature_count": len(scan.features),
            "non_cutter_count": len(scan.non_cutters),
            "type_ii_count": sum(
                1 for f in scan.features if f.feature_family == "restriction_site_II"
            ),
            "type_iis_count": sum(
                1 for f in scan.features if f.feature_family == "restriction_site_IIs"
            ),
            "gateway_count": sum(
                1 for f in scan.features if f.feature_family == "gateway_att"
            ),
            "pcr_warning_count": sum(
                1 for f in scan.features if f.feature_family == "primer_design_warning"
            ),
        },
    }


