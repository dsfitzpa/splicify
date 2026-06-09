"""Plasmid design linter — deterministic biological sense validation.

Implements splicify.intent_rules.v1 adapted to the Splicify module library.

Design-level checks (whole-design):
  - R_POL2_HAS_PROMOTER      Pol II expression requires a Pol II promoter
  - R_POL2_HAS_POLYA         Pol II expression requires a polyA signal
  - R_POL2_HAS_CDS           Pol II expression requires a CDS
  - R_POL2_KOZAK_RECOMMENDED Kozak-like sequence recommended near CDS
  - R_NUCLEAR_PROTEIN_NLS    Nuclear editors (Cas9 etc.) require NLS element(s)
  - R_POL3_HAS_PROMOTER      Guide cassette requires a Pol III promoter
  - R_POL3_HAS_SCAFFOLD      Guide cassette requires gRNA scaffold
  - R_BACTERIAL_HAS_ORI      Plasmid requires bacterial origin of replication
  - R_BACTERIAL_HAS_MARKER   Plasmid requires bacterial selection marker

Sequence-level checks (per resolved CDS-only module):
  - R_CDS_HAS_START          CDS must start with ATG
  - R_CDS_HAS_STOP           CDS must end with a stop codon (unless 2A fusion)
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

# ──────────────────────────────────────────────────────────────────────────────
# Abstract role classification maps
# ──────────────────────────────────────────────────────────────────────────────

# canonical_id.lower() → list of abstract roles
_CANONICAL_ROLES: Dict[str, List[str]] = {
    # ── Pol II promoters ───────────────────────────────────────────────────────
    "promoter_cmv":               ["promoter_pol2"],
    "promoter_ef_1_alpha":        ["promoter_pol2"],
    "promoter_ef_1_alpha_core":   ["promoter_pol2"],
    "promoter_pgk":               ["promoter_pol2"],
    "promoter_cag":               ["promoter_pol2"],
    "promoter_sv40":              ["promoter_pol2"],
    "promoter_rsv":               ["promoter_pol2"],
    "promoter_sffv":              ["promoter_pol2"],
    "promoter_ubic":              ["promoter_pol2"],
    "promoter_ubc":               ["promoter_pol2"],
    "promoter_t7":                ["promoter_pol2"],
    "promoter_chicken_beta_actin":["promoter_pol2"],
    # ── Pol III promoters ──────────────────────────────────────────────────────
    "promoter_u6_human":          ["promoter_pol3"],
    "promoter_h1":                ["promoter_pol3"],
    "promoter_atu6_26":           ["promoter_pol3"],
    "promoter_ceu6":              ["promoter_pol3"],
    "promoter_du6":               ["promoter_pol3"],
    "promoter_du6_1":             ["promoter_pol3"],
    "promoter_du6_3":             ["promoter_pol3"],
    "promoter_tremod_u6":         ["promoter_pol3"],
    "promoter_u6_tet":            ["promoter_pol3"],
    # ── PolyA / Pol II terminators ─────────────────────────────────────────────
    "polya_bgh":                  ["terminator_pol2"],
    "polya_sv40":                 ["terminator_pol2"],
    "polya_generic":              ["terminator_pol2"],
    "polya_hgh":                  ["terminator_pol2"],
    "polya_hgh":                  ["terminator_pol2"],
    "polya_tk":                   ["terminator_pol2"],
    "polya_beta_globin":          ["terminator_pol2"],
    "polya_unlabeled":            ["terminator_pol2"],
    # ── Pol III terminators ────────────────────────────────────────────────────
    "terminator_poliii":          ["terminator_pol3"],
    # ── NLS elements ──────────────────────────────────────────────────────────
    "cds_nls":                    ["nls_like"],
    "cds_sv40_nls":               ["nls_like"],
    "cds_c_myc_nls":              ["nls_like"],
    "cds_egl_13_nls":             ["nls_like"],
    "cds_nucleoplasmin_nls":      ["nls_like"],
    "cds_rex_nls":                ["nls_like"],
    # ── Kozak sequences ────────────────────────────────────────────────────────
    "cds_kozak_sequence":         ["kozak_like"],
    "misc_kozak":                 ["kozak_like"],
    "cds_atg":                    ["kozak_like"],   # explicit ATG start site
    # ── gRNA scaffold ─────────────────────────────────────────────────────────
    "misc_grna_scaffold":         ["grna_scaffold_like"],
    "cds_grna_scaffold":          ["grna_scaffold_like"],
    "cds_grna_scaffold_3_portion":["grna_scaffold_like"],
    # ── 2A self-cleaving peptides ──────────────────────────────────────────────
    "cds_p2a":                    ["peptide_2a_like"],
    "cds_t2a":                    ["peptide_2a_like"],
    "cds_e2a":                    ["peptide_2a_like"],
    "cds_f2a":                    ["peptide_2a_like"],
    # ── IRES ──────────────────────────────────────────────────────────────────
    "cds_ires":                   ["ires_like"],
    "cds_ires2":                  ["ires_like"],
    # ── Bacterial origins of replication ──────────────────────────────────────
    "ori_generic":                ["bacterial_origin"],
    "ori_f1":                     ["bacterial_origin"],
    # ── Bacterial selection markers ────────────────────────────────────────────
    "marker_amp":                 ["bacterial_selection_marker"],
    "marker_kan":                 ["bacterial_selection_marker"],
    "cds_amp":                    ["bacterial_selection_marker"],
    # ── Mammalian selection markers (also CDS) ─────────────────────────────────
    "cds_puror":                  ["mammalian_selection_marker", "cds"],
    "cds_bsd":                    ["mammalian_selection_marker", "cds"],
    "cds_neor_kanr":              ["mammalian_selection_marker", "cds"],
    "cds_hygror":                 ["mammalian_selection_marker", "cds"],
    # ── Reporters (also CDS) ───────────────────────────────────────────────────
    "cds_egfp":                   ["reporter", "cds"],
    "cds_mcherry":                ["reporter", "cds"],
    "cds_tagrfp":                 ["reporter", "cds"],
    "cds_taggfp2":                ["reporter", "cds"],
    "cds_tagbfp":                 ["reporter", "cds"],
    "cds_zsgreen":                ["reporter", "cds"],
    "cds_acgfp1":                 ["reporter", "cds"],
    "cds_cfp":                    ["reporter", "cds"],
    "cds_yfp":                    ["reporter", "cds"],
    "cds_bfp":                    ["reporter", "cds"],
    "cds_citrine":                ["reporter", "cds"],
    # ── Nucleases (CDS + requires NLS) ────────────────────────────────────────
    "cds_cas9":                   ["cds", "nuclease_like"],
    "cds_dcas9":                  ["cds", "nuclease_like"],
    "cds_spcas9_d10a":            ["cds", "nuclease_like"],
    "cds_spcas9_hf1":             ["cds", "nuclease_like"],
    "cds_sacas9":                 ["cds", "nuclease_like"],
    "cds_ascpf1":                 ["cds", "nuclease_like"],
    "cds_espcas9_1_1":            ["cds", "nuclease_like"],
    "cds_nmcas9":                 ["cds", "nuclease_like"],
    "cds_icas9":                  ["cds", "nuclease_like"],
    "cds_cas9_vqr":               ["cds", "nuclease_like"],
    "cds_cas9_vrer":              ["cds", "nuclease_like"],
    "cds_cas9_h840a":             ["cds", "nuclease_like"],
    "cds_cas9m4":                 ["cds", "nuclease_like"],
    "cds_nmdcas9":                ["cds", "nuclease_like"],
    "cds_pcocas9":                ["cds", "nuclease_like"],
    "cds_st1cas9":                ["cds", "nuclease_like"],
    "cds_st1dcas9":               ["cds", "nuclease_like"],
    # ── Lentiviral elements ────────────────────────────────────────────────────
    "cds_wpre":                   ["viral_element_like"],
    "lenti_element_5_ltr":        ["viral_element_like"],
    "lenti_element_3_ltr_delta_u3": ["viral_element_like"],
    "cds_hiv_1_psi":              ["viral_element_like"],
    "cds_rre":                    ["viral_element_like"],
    "cds_cppt_cts":               ["viral_element_like"],
}

# LLM-assigned role → abstract roles (fallback when canonical_id not in map)
_LLM_ROLE_FALLBACK: Dict[str, List[str]] = {
    "promoter":           ["promoter_pol2"],
    "transgene":          ["cds"],
    "reporter":           ["reporter", "cds"],
    "polya":              ["terminator_pol2"],
    "selection_marker":   [],          # resolved below by canonical prefix
    "origin":             ["bacterial_origin"],
    # guide_cassette = complete pol3 unit: all roles implied
    "guide_cassette":     [
        "promoter_pol3", "guide_insert_region_like",
        "grna_scaffold_like", "terminator_pol3",
    ],
    "nuclease":           ["cds", "nuclease_like"],
    "lentiviral_element": ["viral_element_like"],
    "backbone":           ["bacterial_origin", "bacterial_selection_marker"],
    "other":              [],
}

# Module types that are SELF-CONTAINED expression units (include embedded
# promoter / polyA / NLS — structural checks skip these)
_SELF_CONTAINED_TYPES: Set[str] = {
    "nuclease_expression_cassette",
    "pol2_expression_cassette",
    "pol3_u6_sgrna_cassette",
    "bacterial_backbone",
    "bacterial_marker_cassette",
    "bacterial_expression_cassette",
    "sv40_neo_selection_cassette",
    "lentiviral_backbone",
    "mammalian_lentiviral_expression_module",
    "lentiviral_expression_vector",
}

# Canonical IDs whose products must reside in the nucleus
_NUCLEAR_CANONICAL: Set[str] = {
    "cds_cas9", "cds_dcas9", "cds_spcas9_d10a", "cds_spcas9_hf1",
    "cds_sacas9", "cds_ascpf1", "cds_espcas9_1_1", "cds_nmcas9",
    "cds_icas9", "cds_cas9_vqr", "cds_cas9_vrer", "cds_cas9_h840a",
    "cds_cas9m4", "cds_nmdcas9", "cds_pcocas9", "cds_st1cas9", "cds_st1dcas9",
}

_STOP_CODONS: Set[str] = {"TAA", "TAG", "TGA"}


# ──────────────────────────────────────────────────────────────────────────────
# Role helpers
# ──────────────────────────────────────────────────────────────────────────────

def _abstract_roles(mod: Dict[str, Any]) -> List[str]:
    """Return list of abstract roles for a module spec dict."""
    canonical = (mod.get("canonical_id") or "").lower()
    llm_role = mod.get("role", "other")

    if canonical in _CANONICAL_ROLES:
        return list(_CANONICAL_ROLES[canonical])

    # selection_marker: distinguish bacterial vs mammalian
    if llm_role == "selection_marker":
        if canonical.startswith("marker_") or canonical in ("cds_amp",):
            return ["bacterial_selection_marker"]
        return ["mammalian_selection_marker", "cds"]

    # Generic CDS_ patterns not in map
    if canonical.startswith("cds_"):
        return ["cds"]

    return list(_LLM_ROLE_FALLBACK.get(llm_role, []))


def _design_roles(modules: List[Dict[str, Any]]) -> Set[str]:
    """Union of all abstract roles across all modules in a design."""
    roles: Set[str] = set()
    for m in modules:
        roles.update(_abstract_roles(m))
    return roles


def _resolved_module_type(resolved_mod: Dict[str, Any]) -> str:
    """Return the DB module_type of the best resolved candidate, or ''."""
    cands = resolved_mod.get("db_candidates") or []
    return cands[0].get("module_type", "") if cands else ""


def _is_self_contained(resolved_mod: Dict[str, Any]) -> bool:
    return _resolved_module_type(resolved_mod) in _SELF_CONTAINED_TYPES


def _has_nuclear_requirement(modules: List[Dict[str, Any]]) -> bool:
    for mod in modules:
        canonical = (mod.get("canonical_id") or "").lower()
        if canonical in _NUCLEAR_CANONICAL:
            return True
        if mod.get("role") == "nuclease":
            return True
    return False


# ──────────────────────────────────────────────────────────────────────────────
# Intent detection
# ──────────────────────────────────────────────────────────────────────────────

def detect_intent(design_spec: Dict[str, Any]) -> str:
    """Map LLM design_type + module contents → intent_id from intent rules."""
    design_type = (design_spec.get("design_type") or "custom").lower()
    modules = design_spec.get("modules", [])

    has_nuclease = any(
        "nuclease_like" in _abstract_roles(m) or m.get("role") == "nuclease"
        for m in modules
    )
    has_guide = any(m.get("role") == "guide_cassette" for m in modules)
    has_lentiviral = any(m.get("role") == "lentiviral_element" for m in modules)

    if (has_nuclease or design_type == "crispr") and has_guide:
        return "crispr_cas9_mammalian_all_in_one"
    if has_lentiviral or design_type == "lentiviral":
        return "viral_vector_lentiviral_like"
    if design_type == "bacterial_expression":
        return "bacterial_protein_expression_generic"
    if design_type == "aav":
        return "viral_vector_aav_like"
    return "mammalian_protein_expression"


# ──────────────────────────────────────────────────────────────────────────────
# Violation dataclass
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class LintViolation:
    rule_id: str
    severity: str          # ERROR | WARN | INFO
    title: str
    message: str
    evidence: Dict[str, Any] = field(default_factory=dict)
    autofix_action: str = "none"    # add_module | add_nls | flag_user | none
    autofix_params: Dict[str, Any] = field(default_factory=dict)
    applied: bool = False


# ──────────────────────────────────────────────────────────────────────────────
# Structural validation (design-level, pre/post-resolution)
# ──────────────────────────────────────────────────────────────────────────────

def roles_from_canonical_ids(canonical_ids) -> Set[str]:
    """Translate a collection of canonical ID strings to their abstract roles."""
    return {role for cid in canonical_ids for role in _CANONICAL_ROLES.get(cid.lower(), [])}


def _check_structural(
    design_spec: Dict[str, Any],
    resolved_modules: Optional[List[Dict[str, Any]]],
    extra_covered_roles: Optional[Set[str]] = None,
) -> List[LintViolation]:
    """
    Check design-level biological requirements.

    When resolved_modules is provided, self-contained cassette modules
    (nuclease_expression_cassette, pol3_u6_sgrna_cassette, etc.) satisfy
    their own promoter / polyA / NLS requirements automatically.

    extra_covered_roles: abstract roles already satisfied by base plasmid modules
    (set-cover path). These suppress spurious R_BACTERIAL_HAS_ORI / _MARKER violations.
    """
    modules = design_spec.get("modules", [])
    violations: List[LintViolation] = []

    # Build per-module metadata: abstract roles + whether self-contained
    enriched: List[Dict[str, Any]] = []
    for i, mod in enumerate(modules):
        res = resolved_modules[i] if (resolved_modules and i < len(resolved_modules)) else {}
        sc = _is_self_contained(res)
        enriched.append({
            "mod": mod,
            "roles": set(_abstract_roles(mod)),
            "self_contained": sc,
            "idx": i,
        })

    # Aggregate abstract roles of non-self-contained modules
    open_roles: Set[str] = set()
    # All roles (including self-contained) — used for certain design checks
    all_roles: Set[str] = set()
    for e in enriched:
        all_roles.update(e["roles"])
        if not e["self_contained"]:
            open_roles.update(e["roles"])

    # Merge in roles already satisfied by base plasmid modules from set-cover
    if extra_covered_roles:
        all_roles.update(extra_covered_roles)

    # Detect whether there are non-self-contained CDS modules
    has_open_cds = any(
        r in ("cds", "reporter", "nuclease_like")
        for r in open_roles
    )
    has_any_cds = any(
        r in ("cds", "reporter", "nuclease_like")
        for r in all_roles
    )
    has_guide_cassette = any(
        e["roles"] & {"promoter_pol3", "guide_insert_region_like", "grna_scaffold_like"}
        or e["mod"].get("role") == "guide_cassette"
        for e in enriched
    )
    # Self-contained guide cassettes (pol3_u6_sgrna_cassette) satisfy pol3 checks
    has_sc_guide = any(
        e["self_contained"] and e["mod"].get("role") == "guide_cassette"
        for e in enriched
    )
    # Self-contained nuclease cassettes (nuclease_expression_cassette) satisfy NLS etc.
    has_sc_nuclease = any(
        e["self_contained"] and "nuclease_like" in e["roles"]
        for e in enriched
    )

    intent_id = detect_intent(design_spec)
    has_nuclear_req = _has_nuclear_requirement(modules)

    # ── R_POL2_HAS_PROMOTER ──────────────────────────────────────────────────
    if has_open_cds and "promoter_pol2" not in open_roles:
        violations.append(LintViolation(
            rule_id="R_POL2_HAS_PROMOTER",
            severity="ERROR",
            title="Missing Pol II promoter",
            message=(
                "A Pol II expression cassette contains a CDS but no Pol II promoter was found. "
                "Adding CMV promoter."
            ),
            autofix_action="add_module",
            autofix_params={
                "role": "promoter",
                "canonical_id": "PROMOTER_CMV",
                "description": "CMV promoter (auto-inferred)",
                "position": "before_first_cds",
            },
        ))

    # ── R_POL2_HAS_POLYA ─────────────────────────────────────────────────────
    if has_open_cds and "terminator_pol2" not in open_roles:
        violations.append(LintViolation(
            rule_id="R_POL2_HAS_POLYA",
            severity="ERROR",
            title="Missing polyA signal",
            message=(
                "A Pol II expression cassette contains a CDS but no polyA/termination element "
                "was found. Adding BGH polyA."
            ),
            autofix_action="add_module",
            autofix_params={
                "role": "polya",
                "canonical_id": "POLYA_BGH",
                "description": "BGH polyA signal (auto-inferred)",
                "position": "before_pol3_or_backbone",
            },
        ))

    # ── R_POL2_HAS_CDS ───────────────────────────────────────────────────────
    if "promoter_pol2" in open_roles and not has_open_cds:
        violations.append(LintViolation(
            rule_id="R_POL2_HAS_CDS",
            severity="ERROR",
            title="Pol II promoter without CDS",
            message="A Pol II promoter is present but no CDS or reporter was found in the design.",
            autofix_action="flag_user",
            autofix_params={"reason": "CDS is design-specific and cannot be auto-inferred."},
        ))

    # ── R_POL2_KOZAK_RECOMMENDED ─────────────────────────────────────────────
    if has_open_cds and "kozak_like" not in all_roles:
        violations.append(LintViolation(
            rule_id="R_POL2_KOZAK_RECOMMENDED",
            severity="WARN",
            title="Kozak sequence recommended",
            message=(
                "No Kozak-like sequence found near the CDS start. "
                "Expression efficiency may be reduced. Adding Kozak sequence."
            ),
            autofix_action="add_module",
            autofix_params={
                "role": "other",
                "canonical_id": "MISC_kozak",
                "description": "Kozak consensus sequence (auto-inferred)",
                "position": "after_promoter",
            },
        ))

    # ── R_NUCLEAR_PROTEIN_NLS_REQUIRED ───────────────────────────────────────
    if has_nuclear_req and not has_sc_nuclease:
        if "nls_like" not in all_roles:
            violations.append(LintViolation(
                rule_id="R_NUCLEAR_PROTEIN_NLS_REQUIRED",
                severity="ERROR",
                title="NLS required for nuclear protein",
                message=(
                    "Design contains a nuclear-localized editor (Cas9/dCas9/etc.) "
                    "but no NLS element is present. Adding dual NLS (N- and C-terminal)."
                ),
                autofix_action="add_nls",
                autofix_params={
                    "nls_n_term": "CDS_sv40_nls",
                    "nls_c_term": "CDS_nucleoplasmin_nls",
                    "reason": "Dual NLS (SV40 N-terminal, nucleoplasmin C-terminal) "
                              "ensures robust nuclear localization of the editor.",
                },
            ))

    # ── R_POL3_HAS_PROMOTER ──────────────────────────────────────────────────
    if has_guide_cassette and not has_sc_guide:
        if "promoter_pol3" not in all_roles:
            violations.append(LintViolation(
                rule_id="R_POL3_HAS_PROMOTER",
                severity="ERROR",
                title="Missing Pol III promoter in guide cassette",
                message=(
                    "A guide/gRNA cassette is present but no Pol III promoter (e.g., U6) was found. "
                    "Adding human U6 promoter."
                ),
                autofix_action="add_module",
                autofix_params={
                    "role": "promoter",
                    "canonical_id": "PROMOTER_U6_HUMAN",
                    "description": "Human U6 Pol III promoter (auto-inferred)",
                    "position": "before_guide",
                },
            ))

    # ── R_POL3_HAS_SCAFFOLD ──────────────────────────────────────────────────
    if has_guide_cassette and not has_sc_guide:
        if "grna_scaffold_like" not in all_roles:
            violations.append(LintViolation(
                rule_id="R_POL3_HAS_SCAFFOLD",
                severity="ERROR",
                title="Missing gRNA scaffold",
                message=(
                    "A guide cassette is present but no gRNA scaffold was found. "
                    "The scaffold is required for guide RNA function."
                ),
                autofix_action="flag_user",
                autofix_params={"reason": "gRNA scaffold is system-specific; verify the guide cassette module."},
            ))

    # ── R_BACTERIAL_HAS_ORI ──────────────────────────────────────────────────
    if "bacterial_origin" not in all_roles:
        violations.append(LintViolation(
            rule_id="R_BACTERIAL_HAS_ORI",
            severity="ERROR",
            title="Missing bacterial origin of replication",
            message=(
                "No bacterial origin of replication found. "
                "The plasmid cannot propagate in E. coli. Adding ColE1/pUC origin."
            ),
            autofix_action="add_module",
            autofix_params={
                "role": "origin",
                "canonical_id": "ORI_GENERIC",
                "description": "ColE1/pUC origin of replication (auto-inferred)",
                "position": "end",
            },
        ))

    # ── R_BACTERIAL_HAS_MARKER ───────────────────────────────────────────────
    if "bacterial_selection_marker" not in all_roles:
        violations.append(LintViolation(
            rule_id="R_BACTERIAL_HAS_MARKER",
            severity="ERROR",
            title="Missing bacterial selection marker",
            message=(
                "No bacterial selection marker (AmpR, KanR, etc.) found. "
                "Adding ampicillin resistance."
            ),
            autofix_action="add_module",
            autofix_params={
                "role": "selection_marker",
                "canonical_id": "MARKER_AMP",
                "description": "Ampicillin resistance marker (auto-inferred)",
                "position": "end",
            },
        ))

    return violations


# ──────────────────────────────────────────────────────────────────────────────
# Sequence-level validation (post-resolution)
# ──────────────────────────────────────────────────────────────────────────────

def check_sequences(
    modules: List[Dict[str, Any]],
    resolved_modules: List[Dict[str, Any]],
) -> List[LintViolation]:
    """
    Per-CDS sequence checks: start codon and stop codon.
    Only applied to cds_only module types with actual sequences.
    """
    violations: List[LintViolation] = []
    n = min(len(modules), len(resolved_modules))

    for i in range(n):
        mod = modules[i]
        resolved = resolved_modules[i]

        # Only check modules that have a CDS role
        roles = set(_abstract_roles(mod))
        if not any(r in roles for r in ("cds", "reporter", "nuclease_like")):
            continue

        seq = (resolved.get("sequence") or "").upper().strip()
        if len(seq) < 6:
            continue

        # Only check cds_only modules — full cassettes have embedded context
        mt = _resolved_module_type(resolved)
        if mt and mt != "cds_only":
            continue

        desc = mod.get("description") or resolved.get("description") or f"module {i+1}"

        # Check if next module is a 2A peptide or IRES (means no stop codon expected)
        next_mod = modules[i + 1] if i + 1 < len(modules) else {}
        next_roles = set(_abstract_roles(next_mod))
        is_fusion = bool(next_roles & {"peptide_2a_like", "ires_like"})

        # R_CDS_HAS_START
        if seq[:3] != "ATG":
            violations.append(LintViolation(
                rule_id="R_CDS_HAS_START",
                severity="ERROR",
                title="CDS missing start codon",
                message=(
                    f"'{desc}' does not begin with ATG (first codon: '{seq[:3]}'). "
                    "This CDS cannot be translated."
                ),
                evidence={"first_codon": seq[:3], "module_index": i},
                autofix_action="flag_user",
                autofix_params={
                    "reason": "Verify the CDS source; if a Kozak+ATG is embedded in the cassette "
                              "context this may be a module boundary artifact."
                },
            ))

        # R_CDS_HAS_STOP
        if seq[-3:] not in _STOP_CODONS and not is_fusion:
            violations.append(LintViolation(
                rule_id="R_CDS_HAS_STOP",
                severity="ERROR",
                title="CDS missing stop codon",
                message=(
                    f"'{desc}' does not end with a stop codon (last codon: '{seq[-3:]}'). "
                    "Translation will read through into downstream sequence."
                ),
                evidence={"last_codon": seq[-3:], "module_index": i},
                autofix_action="flag_user",
                autofix_params={
                    "reason": "If this CDS is intentionally fused to a downstream element "
                              "(tag, reporter), add a 2A/IRES module to indicate the junction."
                },
            ))
        elif seq[-3:] not in _STOP_CODONS and is_fusion:
            # Informational: stop codon absent but 2A/IRES follows — this is expected
            violations.append(LintViolation(
                rule_id="R_CDS_HAS_STOP",
                severity="INFO",
                title="CDS lacking stop codon (2A/IRES fusion — expected)",
                message=(
                    f"'{desc}' has no stop codon, but is followed by a 2A/IRES element. "
                    "This is the expected arrangement for polyprotein/bicistronic designs."
                ),
                evidence={"last_codon": seq[-3:], "module_index": i},
            ))

    return violations


# ──────────────────────────────────────────────────────────────────────────────
# Insertion position helpers
# ──────────────────────────────────────────────────────────────────────────────

def _find_backbone_start(modules: List[Dict[str, Any]]) -> int:
    """Index of first ORI or bacterial selection marker module."""
    backbone_roles = {"bacterial_origin", "bacterial_selection_marker"}
    for i, mod in enumerate(modules):
        if backbone_roles & set(_abstract_roles(mod)):
            return i
    return len(modules)


def _find_polya_insert_point(modules: List[Dict[str, Any]]) -> int:
    """
    Insertion point for a polyA signal: before the first Pol III element
    (guide_cassette, U6 promoter) or backbone element, whichever comes first.
    This ensures polyA terminates the Pol II cassette before any Pol III cassettes,
    preventing Pol II read-through into the guide RNA region.
    """
    backbone_roles = {"bacterial_origin", "bacterial_selection_marker"}
    for i, mod in enumerate(modules):
        if backbone_roles & set(_abstract_roles(mod)):
            return i
        roles = set(_abstract_roles(mod))
        if mod.get("role") == "guide_cassette" or "promoter_pol3" in roles:
            return i
    return len(modules)


def _find_first_cds_idx(modules: List[Dict[str, Any]]) -> int:
    """Index of first module with a CDS/reporter/nuclease role."""
    for i, mod in enumerate(modules):
        roles = set(_abstract_roles(mod))
        if roles & {"cds", "reporter", "nuclease_like"}:
            return i
    return -1


def _find_pol2_promoter_idx(modules: List[Dict[str, Any]]) -> int:
    """Index of first Pol II promoter module."""
    for i, mod in enumerate(modules):
        if "promoter_pol2" in _abstract_roles(mod):
            return i
    return -1


def _find_nuclease_idx(modules: List[Dict[str, Any]]) -> int:
    """Index of first nuclease-role module."""
    for i, mod in enumerate(modules):
        roles = set(_abstract_roles(mod))
        if "nuclease_like" in roles or mod.get("role") == "nuclease":
            return i
    return -1


def _find_guide_cassette_idx(modules: List[Dict[str, Any]]) -> int:
    """Index of first guide_cassette or pol3 promoter module."""
    for i, mod in enumerate(modules):
        if mod.get("role") == "guide_cassette":
            return i
        if "promoter_pol3" in _abstract_roles(mod):
            return i
    return -1


# ──────────────────────────────────────────────────────────────────────────────
# Autofix application
# ──────────────────────────────────────────────────────────────────────────────

def _make_module(
    role: str,
    canonical_id: str,
    description: str,
    reasoning: str,
) -> Dict[str, Any]:
    return {
        "role": role,
        "canonical_id": canonical_id,
        "description": description,
        "specified": False,
        "reasoning": reasoning,
        "ncbi_query": None,
    }


def apply_autofixes(
    design_spec: Dict[str, Any],
    violations: List[LintViolation],
) -> Tuple[Dict[str, Any], List[str]]:
    """
    Apply fixable autofix actions to design_spec.
    Returns (updated_spec, list_of_fix_summary_strings).
    Marks each applied violation with .applied = True.
    """
    spec = copy.deepcopy(design_spec)
    modules: List[Dict[str, Any]] = spec.setdefault("modules", [])
    applied: List[str] = []

    # Track existing canonical IDs to avoid duplicates
    existing: Set[str] = {(m.get("canonical_id") or "").lower() for m in modules}

    for v in violations:
        if v.applied or v.autofix_action not in ("add_module", "add_nls"):
            continue

        if v.autofix_action == "add_module":
            cid = v.autofix_params.get("canonical_id", "")
            if cid.lower() in existing:
                v.applied = True
                continue

            new_mod = _make_module(
                role=v.autofix_params.get("role", "other"),
                canonical_id=cid,
                description=v.autofix_params.get("description", cid),
                reasoning=v.message,
            )

            position = v.autofix_params.get("position", "before_backbone")
            backbone_i = _find_backbone_start(modules)
            promoter_i = _find_pol2_promoter_idx(modules)
            cds_i = _find_first_cds_idx(modules)
            guide_i = _find_guide_cassette_idx(modules)

            if position == "before_pol3_or_backbone":
                modules.insert(_find_polya_insert_point(modules), new_mod)
            elif position == "before_backbone":
                modules.insert(backbone_i, new_mod)
            elif position == "end":
                modules.append(new_mod)
            elif position == "before_first_cds":
                insert_at = cds_i if cds_i >= 0 else backbone_i
                modules.insert(insert_at, new_mod)
            elif position == "after_promoter":
                # Insert right after the Pol II promoter, or before the first CDS
                if promoter_i >= 0:
                    modules.insert(promoter_i + 1, new_mod)
                elif cds_i >= 0:
                    modules.insert(cds_i, new_mod)
                else:
                    modules.insert(backbone_i, new_mod)
            elif position == "before_guide":
                insert_at = guide_i if guide_i >= 0 else backbone_i
                modules.insert(insert_at, new_mod)
            else:
                modules.insert(backbone_i, new_mod)

            existing.add(cid.lower())
            v.applied = True
            applied.append(
                f"Added **{new_mod['description']}** (`{cid}`) — {v.rule_id}"
            )

        elif v.autofix_action == "add_nls":
            nls_n = v.autofix_params.get("nls_n_term", "CDS_sv40_nls")
            nls_c = v.autofix_params.get("nls_c_term", "CDS_nucleoplasmin_nls")
            nuclease_i = _find_nuclease_idx(modules)
            added_any = False

            # N-terminal NLS: insert before nuclease
            if nls_n.lower() not in existing:
                n_mod = _make_module(
                    role="other",
                    canonical_id=nls_n,
                    description=f"N-terminal NLS ({nls_n.replace('CDS_', '').replace('_', ' ').title()}) (auto-inferred)",
                    reasoning="NLS required for nuclear localization of editor protein",
                )
                insert_at = nuclease_i if nuclease_i >= 0 else _find_backbone_start(modules)
                modules.insert(insert_at, n_mod)
                # Update nuclease index after insertion
                if nuclease_i >= 0:
                    nuclease_i += 1
                existing.add(nls_n.lower())
                applied.append(f"Added **{n_mod['description']}** (`{nls_n}`) — {v.rule_id}")
                added_any = True

            # C-terminal NLS: insert after nuclease
            if nls_c.lower() not in existing:
                c_mod = _make_module(
                    role="other",
                    canonical_id=nls_c,
                    description=f"C-terminal NLS ({nls_c.replace('CDS_', '').replace('_', ' ').title()}) (auto-inferred)",
                    reasoning="Dual NLS improves nuclear import efficiency",
                )
                insert_at = (nuclease_i + 1) if nuclease_i >= 0 else _find_backbone_start(modules)
                modules.insert(insert_at, c_mod)
                existing.add(nls_c.lower())
                applied.append(f"Added **{c_mod['description']}** (`{nls_c}`) — {v.rule_id}")
                added_any = True

            if added_any or (nls_n.lower() in existing and nls_c.lower() in existing):
                v.applied = True

    spec["modules"] = modules
    return spec, applied


# ──────────────────────────────────────────────────────────────────────────────
# Main entry points
# ──────────────────────────────────────────────────────────────────────────────

def lint_and_fix(
    design_spec: Dict[str, Any],
    resolved_modules: Optional[List[Dict[str, Any]]] = None,
    extra_covered_roles: Optional[Set[str]] = None,
) -> Tuple[Dict[str, Any], List[LintViolation], List[str]]:
    """
    Run structural validation and apply autofixes to design_spec.

    Args:
        design_spec: LLM-generated design specification dict.
        resolved_modules: Optional list of resolved modules (from _resolve_modules).
                          When provided, self-contained cassette types are detected
                          and those cassettes are exempt from structural checks.
        extra_covered_roles: Abstract roles already satisfied by base plasmid modules
                             (set-cover path). Prevents spurious ORI/marker autofixes.

    Returns:
        (updated_design_spec, violations, applied_fix_summaries)
    """
    violations = _check_structural(design_spec, resolved_modules, extra_covered_roles)
    updated_spec, applied = apply_autofixes(design_spec, violations)
    return updated_spec, violations, applied


def lint_sequences(
    modules: List[Dict[str, Any]],
    resolved_modules: List[Dict[str, Any]],
) -> List[LintViolation]:
    """
    Run sequence-level checks after module resolution.

    Args:
        modules: Module specs from the (possibly updated) design_spec.
        resolved_modules: Fully resolved modules with sequences.

    Returns:
        List of LintViolation for sequence-level issues.
    """
    return check_sequences(modules, resolved_modules)


# ──────────────────────────────────────────────────────────────────────────────
# Report formatting
# ──────────────────────────────────────────────────────────────────────────────

def format_lint_report(
    violations: List[LintViolation],
    applied_fixes: List[str],
    external_errors: int = 0,
    external_warnings: int = 0,
) -> str:
    """Format linting results as a markdown section for the workflow reply."""
    errors = [v for v in violations if v.severity == "ERROR" and not v.applied]
    warns  = [v for v in violations if v.severity == "WARN"  and not v.applied]
    infos  = [v for v in violations if v.severity == "INFO"  and not v.applied]

    if not violations and not applied_fixes and external_errors == 0 and external_warnings == 0:
        return "\n\n---\n**Design Validation:** ✓ All biological sense checks passed."

    lines = ["\n\n---", "### Design Validation"]

    if applied_fixes:
        lines += ["", "**Auto-corrections applied:**"]
        for fix in applied_fixes:
            lines.append(f"  ✓ {fix}")

    if errors:
        lines += ["", "**Errors — require attention:**"]
        for v in errors:
            lines.append(f"  ✗ `{v.rule_id}` {v.title}")
            lines.append(f"    {v.message}")
            reason = v.autofix_params.get("reason", "")
            if reason:
                lines.append(f"    *Fix:* {reason}")

    if warns:
        lines += ["", "**Warnings:**"]
        for v in warns:
            lines.append(f"  ⚠ `{v.rule_id}` {v.title}")
            lines.append(f"    {v.message}")

    if infos:
        lines += ["", "**Notes:**"]
        for v in infos:
            lines.append(f"  ℹ `{v.rule_id}` {v.title}: {v.message}")

    if external_errors:
        lines += ["", f"  ✗ Analyzer reported {external_errors} blocking conflict(s)."]
    elif external_warnings:
        lines += ["", f"  ⚠ Analyzer reported {external_warnings} warning conflict(s)."]

    if not errors and not warns and external_errors == 0 and external_warnings == 0:
        lines += ["", "  ✓ All biological sense checks passed" +
                  (" (after auto-corrections)." if applied_fixes else ".")]

    return "\n".join(lines)
