"""
Plasmid design analyzer — post-assembly QA layer.

Takes the assembled module map from plasmid_design_chat.py and:
  - Groups modules into expression cassettes
  - Scans for restriction enzyme sites
  - Detects intent conflicts (duplicate nucleases, viral elements in non-viral design, etc.)
  - Generates an LLM-powered purpose summary with subtle mismatch detection
  - Formats a concise markdown analysis report appended to the workflow reply
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .plasmid_design_linter import (
    _CANONICAL_ROLES,
    _NUCLEAR_CANONICAL,
    _abstract_roles,
)
from .utils import reverse_complement


# ──────────────────────────────────────────────────────────────────────────────
# RE site databases
# ──────────────────────────────────────────────────────────────────────────────

_TYPE_II_RE_SITES: Dict[str, str] = {
    "EcoRI":  "GAATTC",
    "BamHI":  "GGATCC",
    "HindIII": "AAGCTT",
    "NheI":   "GCTAGC",
    "XhoI":   "CTCGAG",
    "XbaI":   "TCTAGA",
    "SalI":   "GTCGAC",
    "KpnI":   "GGTACC",
    "SacI":   "GAGCTC",
    "NotI":   "GCGGCCGC",
    "SpeI":   "ACTAGT",
    "ClaI":   "ATCGAT",
    "PstI":   "CTGCAG",
    "AgeI":   "ACCGGT",
    "NcoI":   "CCATGG",
    "NdeI":   "CATATG",
    "BglII":  "AGATCT",
    "AscI":   "GGCGCGCC",
    "PacI":   "TTAATTAA",
    "FseI":   "GGCCGGCC",
    "AvrII":  "CCTAGG",
    "MluI":   "ACGCGT",
    "SbfI":   "CCTGCAGG",
    "SwaI":   "ATTTAAAT",
    "MfeI":   "CAATTG",
}

_TYPE_IIS_RE_SITES: Dict[str, str] = {
    "BsaI":  "GGTCTC",   # Golden Gate (Addgene standard)
    "BsmBI": "CGTCTC",   # Golden Gate (alt)
    "BbsI":  "GAAGAC",   # CRISPR sgRNA guide insertion
    "SapI":  "GCTCTTC",  # Golden Gate (alt)
    "BtgZI": "GCGATG",
    "BfuAI": "ACCTGC",
}

_GATEWAY_SITES: Dict[str, str] = {
    "attB1": "ACAAGTTTGTACAAAAAAGCAGGCT",
    "attB2": "ACCACTTTGTACAAGAAAGCTGGGT",
}

_RECOMBINATION_SITES: Dict[str, str] = {
    "loxP":    "ATAACTTCGTATAGCATACATTATACGAAGTTAT",
    "lox2272": "ATAACTTCGTATAAAGTATATTTATACGAAGTTAT",
    "FRT":     "GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC",
}


# ──────────────────────────────────────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class ExpressionCassette:
    cassette_type: str                   # "pol2" | "pol3" | "bacterial" | "source_plasmid"
    promoter: Optional[str]              # canonical_id of promoter
    cds_list: List[str] = field(default_factory=list)
    nls_list: List[str] = field(default_factory=list)
    reporter: Optional[str] = None
    terminator: Optional[str] = None
    is_complete: bool = False
    description: str = ""
    source_plasmid: Optional[str] = None


@dataclass
class CloningFeature:
    name: str
    feature_type: str       # "type2_re" | "type2s_re" | "gateway" | "cre_lox" | "flp_frt"
    recognition_seq: str
    count: int
    positions: List[int]    # 0-based positions in assembled sequence
    note: Optional[str] = None


@dataclass
class IntentConflict:
    rule_id: str
    severity: str           # "ERROR" | "WARN" | "INFO"
    title: str
    message: str


@dataclass
class PlasmidAnalysis:
    inferred_purpose: str
    inferred_design_type: str
    expression_cassettes: List[ExpressionCassette]
    cloning_features: List[CloningFeature]
    component_summary: Dict[str, List[str]]
    conflicts: List[IntentConflict]
    llm_intent_note: Optional[str] = None


# ──────────────────────────────────────────────────────────────────────────────
# RE site scanning
# ──────────────────────────────────────────────────────────────────────────────

def _iupac_to_regex(seq: str) -> str:
    """Convert IUPAC ambiguity codes to regex character classes."""
    iupac = {
        "N": "[ACGT]", "R": "[AG]", "Y": "[CT]",
        "S": "[GC]",   "W": "[AT]", "K": "[GT]",
        "M": "[AC]",   "B": "[CGT]", "D": "[AGT]",
        "H": "[ACT]",  "V": "[ACG]",
    }
    result = ""
    for ch in seq.upper():
        result += iupac.get(ch, ch)
    return result


def scan_restriction_sites(
    sequence: str,
    assembly_strategy: str = "",
) -> List[CloningFeature]:
    """
    Scan assembled sequence for restriction enzyme recognition sites.
    Returns CloningFeature entries for all sites found (count > 0).
    """
    if not sequence:
        return []

    seq_upper = sequence.upper()
    seq_len = len(seq_upper)

    # Skip if sequence is >80% N-padded (synthesis placeholders)
    n_count = seq_upper.count("N")
    if seq_len > 0 and n_count / seq_len > 0.80:
        return []

    rc_seq = reverse_complement(seq_upper)
    results: List[CloningFeature] = []

    def _scan_sites(sites_db: Dict[str, str], feature_type: str) -> None:
        for name, recog in sites_db.items():
            pattern = _iupac_to_regex(recog)
            positions: List[int] = []
            seen: set = set()

            # Forward strand
            for m in re.finditer(pattern, seq_upper):
                pos = m.start()
                if pos not in seen:
                    positions.append(pos)
                    seen.add(pos)

            # Reverse strand — map back to forward-strand positions
            rc_pattern = _iupac_to_regex(reverse_complement(recog))
            for m in re.finditer(rc_pattern, seq_upper):
                pos = m.start()
                if pos not in seen:
                    positions.append(pos)
                    seen.add(pos)

            if not positions:
                continue

            positions.sort()
            note: Optional[str] = None

            # Assembly warning for Type IIs sites
            if feature_type == "type2s_re" and "golden_gate" in assembly_strategy.lower():
                note = (
                    f"Internal {name} site may cause unintended cleavage "
                    "during Golden Gate assembly."
                )
            # Common informational notes
            if name == "BbsI":
                note = "BbsI detected — commonly used for sgRNA guide insertion into CRISPR vectors."
            elif name in ("BsaI", "BsmBI") and "golden_gate" not in assembly_strategy.lower():
                note = f"{name} detected — Type IIs enzyme, often used in Golden Gate assemblies."

            results.append(CloningFeature(
                name=name,
                feature_type=feature_type,
                recognition_seq=recog,
                count=len(positions),
                positions=positions,
                note=note,
            ))

    _scan_sites(_TYPE_IIS_RE_SITES, "type2s_re")
    _scan_sites(_TYPE_II_RE_SITES, "type2_re")
    _scan_sites(_GATEWAY_SITES, "gateway")

    # Cre-lox and FLP-FRT together
    cre_lox_flp: Dict[str, str] = {}
    cre_lox_flp.update(_RECOMBINATION_SITES)
    for name, recog in cre_lox_flp.items():
        feature_type = "flp_frt" if name == "FRT" else "cre_lox"
        pattern = _iupac_to_regex(recog)
        positions: List[int] = []
        seen: set = set()
        for m in re.finditer(pattern, seq_upper):
            pos = m.start()
            if pos not in seen:
                positions.append(pos)
                seen.add(pos)
        rc_pat = _iupac_to_regex(reverse_complement(recog))
        for m in re.finditer(rc_pat, seq_upper):
            pos = m.start()
            if pos not in seen:
                positions.append(pos)
                seen.add(pos)
        if positions:
            positions.sort()
            results.append(CloningFeature(
                name=name,
                feature_type=feature_type,
                recognition_seq=recog,
                count=len(positions),
                positions=positions,
            ))

    return results


# ──────────────────────────────────────────────────────────────────────────────
# Component extraction
# ──────────────────────────────────────────────────────────────────────────────

def _extract_all_canonical_ids(
    resolved_modules: List[Dict[str, Any]],
) -> Dict[str, List[str]]:
    """
    Build component_summary dict categorised by abstract role.
    Handles both regular modules and base_plasmid source modules
    (which contribute via covers_ids).
    """
    summary: Dict[str, List[str]] = {
        "promoters": [],
        "nucleases": [],
        "reporters": [],
        "cds": [],
        "terminators": [],
        "nls": [],
        "bacterial_elements": [],
        "viral_elements": [],
        "cloning_elements": [],
    }

    def _add_canonical(canonical_id: str, roles: List[str]) -> None:
        cid = canonical_id.lower()
        added = False
        for role in roles:
            if role in ("promoter_pol2", "promoter_pol3"):
                if cid not in summary["promoters"]:
                    summary["promoters"].append(cid)
                added = True
            elif role == "nuclease_like":
                if cid not in summary["nucleases"]:
                    summary["nucleases"].append(cid)
                added = True
            elif role == "reporter":
                if cid not in summary["reporters"]:
                    summary["reporters"].append(cid)
                added = True
            elif role in ("terminator_pol2", "terminator_pol3"):
                if cid not in summary["terminators"]:
                    summary["terminators"].append(cid)
                added = True
            elif role == "nls_like":
                if cid not in summary["nls"]:
                    summary["nls"].append(cid)
                added = True
            elif role in ("bacterial_origin", "bacterial_selection_marker"):
                if cid not in summary["bacterial_elements"]:
                    summary["bacterial_elements"].append(cid)
                added = True
            elif role == "viral_element_like":
                if cid not in summary["viral_elements"]:
                    summary["viral_elements"].append(cid)
                added = True
            elif role in ("kozak_like", "peptide_2a_like", "ires_like"):
                if cid not in summary["cloning_elements"]:
                    summary["cloning_elements"].append(cid)
                added = True

        # CDS that isn't a nuclease or reporter
        if "cds" in roles and "nuclease_like" not in roles and "reporter" not in roles:
            if cid not in summary["cds"]:
                summary["cds"].append(cid)

    for mod in resolved_modules:
        source_type = mod.get("source_type", "")

        if source_type == "base_plasmid":
            # Extract roles from covers_ids list
            for cid in mod.get("covers_ids", []):
                cid_lower = cid.lower()
                roles = list(_CANONICAL_ROLES.get(cid_lower, []))
                if not roles:
                    # Generic cds_ pattern
                    if cid_lower.startswith("cds_"):
                        roles = ["cds"]
                if roles:
                    _add_canonical(cid_lower, roles)
        else:
            # Regular module
            canonical_id = (mod.get("canonical_id") or "").lower()
            if canonical_id:
                roles = list(_CANONICAL_ROLES.get(canonical_id, []))
                if not roles:
                    llm_role = mod.get("role", "other")
                    from .plasmid_design_linter import _LLM_ROLE_FALLBACK
                    roles = list(_LLM_ROLE_FALLBACK.get(llm_role, []))
                    if not roles and canonical_id.startswith("cds_"):
                        roles = ["cds"]
                if roles:
                    _add_canonical(canonical_id, roles)

    return summary


# ──────────────────────────────────────────────────────────────────────────────
# Cassette grouping
# ──────────────────────────────────────────────────────────────────────────────

def _build_cassette_description(cassette: ExpressionCassette) -> str:
    parts = []
    if cassette.promoter:
        parts.append(cassette.promoter.replace("promoter_", "").upper())
    parts.extend(n.replace("cds_", "").upper() for n in cassette.nls_list[:1])
    parts.extend(c.replace("cds_", "").upper() for c in cassette.cds_list)
    if cassette.reporter:
        parts.append(cassette.reporter.replace("cds_", "").upper())
    if cassette.nls_list[1:]:
        parts.extend(n.replace("cds_", "").upper() for n in cassette.nls_list[1:])
    if cassette.terminator:
        parts.append(cassette.terminator.replace("polya_", "").upper() + "-polyA")
    return " → ".join(parts) if parts else cassette.cassette_type


def _group_into_cassettes(
    resolved_modules: List[Dict[str, Any]],
) -> List[ExpressionCassette]:
    """
    State-machine pass over resolved modules to group them into
    expression cassettes.
    """
    cassettes: List[ExpressionCassette] = []
    current: Optional[ExpressionCassette] = None

    def _close(c: ExpressionCassette) -> None:
        c.is_complete = bool(c.promoter and (c.cds_list or c.reporter) and c.terminator)
        c.description = _build_cassette_description(c)
        cassettes.append(c)

    for mod in resolved_modules:
        source_type = mod.get("source_type", "")

        if source_type == "base_plasmid":
            # Treat each base-plasmid module as a self-contained cassette
            if current:
                _close(current)
                current = None
            cassette = ExpressionCassette(
                cassette_type="source_plasmid",
                promoter=None,
                source_plasmid=mod.get("source", ""),
                is_complete=True,
            )
            covers = [c.lower() for c in mod.get("covers_ids", [])]
            cassette.description = mod.get("description", "Base plasmid module")
            if covers:
                cassette.description += f" [{', '.join(covers[:5])}{'...' if len(covers) > 5 else ''}]"
            cassettes.append(cassette)
            continue

        # Determine abstract roles for this module
        roles = set(_abstract_roles(mod))

        if "promoter_pol2" in roles:
            if current and current.cassette_type == "pol2":
                _close(current)
            current = ExpressionCassette(
                cassette_type="pol2",
                promoter=(mod.get("canonical_id") or "").lower(),
            )

        elif "promoter_pol3" in roles:
            if current and current.cassette_type == "pol2":
                _close(current)
                current = None
            current = ExpressionCassette(
                cassette_type="pol3",
                promoter=(mod.get("canonical_id") or "").lower(),
            )

        elif "terminator_pol2" in roles or "terminator_pol3" in roles:
            if current:
                current.terminator = (mod.get("canonical_id") or "").lower()
                _close(current)
                current = None

        elif "nls_like" in roles:
            if current:
                current.nls_list.append((mod.get("canonical_id") or "").lower())

        elif "reporter" in roles:
            if current:
                current.reporter = (mod.get("canonical_id") or "").lower()
            else:
                # Reporter without open cassette
                current = ExpressionCassette(
                    cassette_type="pol2",
                    promoter=None,
                )
                current.reporter = (mod.get("canonical_id") or "").lower()

        elif "nuclease_like" in roles or "cds" in roles:
            cid = (mod.get("canonical_id") or "").lower()
            if current:
                if cid not in current.cds_list:
                    current.cds_list.append(cid)
            else:
                current = ExpressionCassette(
                    cassette_type="pol2",
                    promoter=None,
                )
                current.cds_list.append(cid)

    # Close any still-open cassette
    if current:
        _close(current)

    return cassettes


# ──────────────────────────────────────────────────────────────────────────────
# Conflict detection
# ──────────────────────────────────────────────────────────────────────────────

def _count_canonical_occurrences(
    resolved_modules: List[Dict[str, Any]],
) -> Dict[str, int]:
    """Count how many times each canonical_id appears across all modules."""
    counts: Dict[str, int] = {}

    for mod in resolved_modules:
        source_type = mod.get("source_type", "")

        if source_type == "base_plasmid":
            for cid in mod.get("covers_ids", []):
                key = cid.lower()
                counts[key] = counts.get(key, 0) + 1
        else:
            cid = (mod.get("canonical_id") or "").lower()
            if cid:
                counts[cid] = counts.get(cid, 0) + 1

    return counts


def _assembled_has_guide(resolved_modules: List[Dict[str, Any]]) -> bool:
    """Check if guide cassette elements appear in the assembled modules."""
    guide_roles = {"promoter_pol3", "grna_scaffold_like", "guide_insert_region_like"}
    for mod in resolved_modules:
        if mod.get("role") == "guide_cassette":
            return True
        source_type = mod.get("source_type", "")
        if source_type == "base_plasmid":
            for cid in mod.get("covers_ids", []):
                if guide_roles & set(_CANONICAL_ROLES.get(cid.lower(), [])):
                    return True
        else:
            if guide_roles & set(_abstract_roles(mod)):
                return True
    return False


def _detect_conflicts(
    resolved_modules: List[Dict[str, Any]],
    design_spec: Dict[str, Any],
    component_summary: Dict[str, List[str]],
    cloning_features: Optional[List[CloningFeature]] = None,
) -> List[IntentConflict]:
    """
    Check for conflicts between the assembled design and what was requested.

    Rules are derived from the design_spec (the LLM's interpretation of the user's
    request), not from universal biological prescriptions. This means:
      - Guide/nuclease pairing is only checked when the spec asked for BOTH.
      - A guide-only plasmid is not flagged for lacking a nuclease.
      - A nuclease-only plasmid is not flagged for lacking a guide cassette.
      - Viral backbone mismatch is checked against the requested design_type,
        not against a fixed "non-viral = no lentiviral elements" rule.
    """
    conflicts: List[IntentConflict] = []
    design_type = (design_spec.get("design_type") or "custom").lower()
    assembly_strategy = (design_spec.get("assembly_strategy") or "").lower()
    modules_spec = design_spec.get("modules", [])
    counts = _count_canonical_occurrences(resolved_modules)

    # ── What did the design_spec request? ─────────────────────────────────────
    spec_roles = {m.get("role") for m in modules_spec}
    spec_has_guide    = "guide_cassette" in spec_roles
    spec_has_nuclease = "nuclease" in spec_roles
    spec_is_lentiviral = design_type == "lentiviral" or "lentiviral_element" in spec_roles

    # ── Assembled state ───────────────────────────────────────────────────────
    nucleases = component_summary.get("nucleases", [])
    reporters = component_summary.get("reporters", [])
    viral     = component_summary.get("viral_elements", [])
    bacterial = component_summary.get("bacterial_elements", [])
    assembled_has_guide = _assembled_has_guide(resolved_modules)

    # ── Universal objective rules (always apply regardless of design intent) ──

    # R_DUPLICATE_NUCLEASE — same nuclease from multiple sources; never intentional
    for nuc in nucleases:
        if counts.get(nuc, 0) > 1:
            conflicts.append(IntentConflict(
                rule_id="R_DUPLICATE_NUCLEASE",
                severity="ERROR",
                title="Duplicate nuclease",
                message=(
                    f"`{nuc}` appears {counts[nuc]}× from different source plasmids. "
                    "Only one copy is needed."
                ),
            ))

    # R_MULTIPLE_NUCLEASE_TYPES — multiple distinct nucleases; flag as informational
    if len(nucleases) > 1:
        conflicts.append(IntentConflict(
            rule_id="R_MULTIPLE_NUCLEASE_TYPES",
            severity="INFO",
            title="Multiple nuclease types",
            message=f"Multiple nuclease types in assembly: {', '.join(nucleases)}.",
        ))

    # R_DUPLICATE_REPORTER — same reporter from multiple sources
    for rep in reporters:
        if counts.get(rep, 0) > 1:
            conflicts.append(IntentConflict(
                rule_id="R_DUPLICATE_REPORTER",
                severity="ERROR",
                title="Duplicate reporter",
                message=f"`{rep}` reporter appears {counts[rep]}× from different source plasmids.",
            ))

    # R_DUAL_BACTERIAL_ORIGIN — only fire for two HIGH-COPY, incompatible origins.
    # ori_f1 (M13/f1 phage origin) and ori_sv40 (mammalian origin) are routinely
    # combined with ColE1-type origins in phagemid and shuttle vectors without
    # causing instability.  Only ColE1-family origins conflict with each other.
    _HIGH_COPY_ORIGINS = frozenset({"ori_generic", "ori_cole1", "ori_puc", "ori_p15a"})
    _origin_set = {k.lower() for k, v in _CANONICAL_ROLES.items() if "bacterial_origin" in v}
    distinct_origins = list(dict.fromkeys(b for b in bacterial if b in _origin_set))
    distinct_high_copy = [o for o in distinct_origins if o in _HIGH_COPY_ORIGINS]
    if len(distinct_high_copy) > 1:
        conflicts.append(IntentConflict(
            rule_id="R_DUAL_BACTERIAL_ORIGIN",
            severity="ERROR",
            title="Two high-copy bacterial origins of replication",
            message=(
                f"Two ColE1-compatible origins present: {', '.join(distinct_high_copy)}. "
                "This causes replication competition and plasmid instability in E. coli."
            ),
        ))

    # R_DUPLICATE_MARKER — same selection marker from multiple sources
    for mk in [b for b in bacterial if b not in distinct_origins]:
        if counts.get(mk, 0) > 1:
            conflicts.append(IntentConflict(
                rule_id="R_DUPLICATE_MARKER",
                severity="WARN",
                title="Duplicate bacterial selection marker",
                message=f"`{mk}` bacterial marker appears {counts[mk]}×.",
            ))

    # R_TYPEIIS_INTERNAL — internal Type IIs site in a Golden Gate design
    if cloning_features and "golden_gate" in assembly_strategy:
        typeiis_found = [cf for cf in cloning_features if cf.feature_type == "type2s_re"]
        if typeiis_found:
            names = ", ".join(cf.name for cf in typeiis_found)
            conflicts.append(IntentConflict(
                rule_id="R_TYPEIIS_INTERNAL",
                severity="WARN",
                title="Type IIs site in Golden Gate assembly",
                message=(
                    f"{names} site(s) found in assembled sequence. "
                    "Internal site(s) may cause unintended cleavage during Golden Gate assembly."
                ),
            ))

    # ── Request-derived rules (only apply when relevant to the stated intent) ─

    # R_GUIDE_NUCLEASE_PAIRING — only check pairing when BOTH were requested.
    # A guide-only or nuclease-only design is a legitimate user choice.
    if spec_has_nuclease and spec_has_guide:
        if not assembled_has_guide:
            conflicts.append(IntentConflict(
                rule_id="R_GUIDE_CASSETTE_MISSING",
                severity="WARN",
                title="Guide cassette not found in assembly",
                message=(
                    "The design requested a guide cassette alongside the nuclease, "
                    "but no guide cassette elements were found in the assembled result."
                ),
            ))
        if not nucleases:
            conflicts.append(IntentConflict(
                rule_id="R_NUCLEASE_MISSING",
                severity="WARN",
                title="Nuclease not found in assembly",
                message=(
                    "The design requested a nuclease alongside a guide cassette, "
                    "but no nuclease CDS was found in the assembled result."
                ),
            ))

    # R_LENTIVIRAL_BACKBONE — check viral elements match the requested vector type.
    # Only fires when there is a clear mismatch between request and assembly.
    if spec_is_lentiviral and not viral:
        conflicts.append(IntentConflict(
            rule_id="R_LENTIVIRAL_ELEMENTS_MISSING",
            severity="WARN",
            title="Lentiviral elements not assembled",
            message=(
                "The design specified a lentiviral vector, but no lentiviral backbone "
                "elements (LTRs, Ψ, RRE, cPPT, WPRE) were found in the assembly."
            ),
        ))
    elif not spec_is_lentiviral and viral:
        # Viral elements appeared in a non-lentiviral design — likely pulled in by the
        # set-cover optimizer from a source plasmid that covers needed components.
        conflicts.append(IntentConflict(
            rule_id="R_UNEXPECTED_VIRAL_ELEMENTS",
            severity="WARN",
            title="Unexpected lentiviral elements in assembly",
            message=(
                f"Lentiviral elements ({', '.join(viral)}) were not requested but "
                "appear in the assembly, likely from a source plasmid selected by the "
                "fragment optimizer. Verify the backbone is appropriate for your application."
            ),
        ))

    return conflicts


# ──────────────────────────────────────────────────────────────────────────────
# Design type inference
# ──────────────────────────────────────────────────────────────────────────────

def _infer_design_type(
    resolved_modules: List[Dict[str, Any]],
    component_summary: Dict[str, List[str]],
) -> str:
    """Deterministically infer the design type from component summary."""
    viral = component_summary.get("viral_elements", [])
    nucleases = component_summary.get("nucleases", [])
    bacterial = component_summary.get("bacterial_elements", [])

    has_guide = any(
        mod.get("role") == "guide_cassette"
        for mod in resolved_modules
    )

    if viral:
        return "lentiviral"
    if nucleases and has_guide:
        return "crispr"
    if nucleases:
        return "nuclease_expression"
    if bacterial and not component_summary.get("promoters") and not component_summary.get("cds"):
        return "bacterial_expression"
    return "mammalian_expression"


# ──────────────────────────────────────────────────────────────────────────────
# LLM purpose summary
# ──────────────────────────────────────────────────────────────────────────────

async def _llm_purpose_summary(
    user_request: str,
    component_summary: Dict[str, List[str]],
    inferred_design_type: str,
    resolved_modules: Optional[List[Any]] = None,
) -> Dict[str, Optional[str]]:
    """
    Call gpt-4o-mini for a one-sentence purpose summary and optional
    intent mismatch note. Returns {inferred_purpose, intent_mismatch}.
    """
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        return {"inferred_purpose": None, "intent_mismatch": None}

    component_lines = "\n".join(
        f"  {k}: {', '.join(v)}"
        for k, v in component_summary.items()
        if v
    )
    
    # Build module hierarchy if available
    module_summary = ""
    if resolved_modules:
        module_summary = _build_module_hierarchy_summary(resolved_modules)
        module_summary = f"\n\nModule hierarchy:\n{module_summary}"

    prompt = (
        f"A molecular biology plasmid has been assembled with the following components:\n"
        f"{component_lines}\n"
        f"Inferred design type: {inferred_design_type}{module_summary}\n"
        f"Original user request: \"{user_request}\"\n\n"
        "Respond in JSON with exactly these keys:\n"
        "  \"inferred_purpose\": A single sentence (≤25 words) describing what this plasmid does.\n"
        "  \"intent_mismatch\": A short note (≤20 words) if the assembled components do NOT "
        "match what the user requested, or null if the design matches the request."
    )

    try:
        from openai import AsyncOpenAI
        client = AsyncOpenAI(api_key=api_key)
        response = await client.chat.completions.create(
            model="gpt-4o",
            messages=[{"role": "system", "content": MODULE_SEMANTICS_SYSTEM_PROMPT}, {"role": "user", "content": prompt}],
            response_format={"type": "json_object"},
            temperature=0.3,
            max_tokens=250,
        )
        import json
        data = json.loads(response.choices[0].message.content)
        return {
            "inferred_purpose": data.get("inferred_purpose"),
            "intent_mismatch": data.get("intent_mismatch"),
        }
    except Exception as e:
        print(f"[plasmid_analyzer] LLM call failed: {e}")
        return {"inferred_purpose": None, "intent_mismatch": None}


def _fallback_purpose(
    inferred_design_type: str,
    component_summary: Dict[str, List[str]],
) -> str:
    """Generate a description from inferred_design_type + component_summary without LLM."""
    nucleases = component_summary.get("nucleases", [])
    reporters = component_summary.get("reporters", [])
    cds = component_summary.get("cds", [])

    if inferred_design_type == "crispr":
        nuc = nucleases[0].replace("cds_", "").upper() if nucleases else "nuclease"
        return f"All-in-one CRISPR vector expressing {nuc} with guide RNA."
    if inferred_design_type == "lentiviral":
        payload = (reporters + cds)
        payload_str = reporters[0].replace("cds_", "").upper() if reporters else "transgene"
        return f"Lentiviral expression vector delivering {payload_str}."
    if reporters:
        rep = reporters[0].replace("cds_", "").upper()
        return f"Mammalian expression plasmid for {rep} reporter."
    if cds:
        gene = cds[0].replace("cds_", "").upper()
        return f"Mammalian expression plasmid for {gene}."
    return "Plasmid design assembled from module library."


# ──────────────────────────────────────────────────────────────────────────────
# Main analysis entry point
# ──────────────────────────────────────────────────────────────────────────────

async def analyze_plasmid_design(
    resolved_modules: List[Dict[str, Any]],
    design_spec: Dict[str, Any],
    session_messages: List[Dict[str, str]],
    assembled_sequence: str,
) -> PlasmidAnalysis:
    """
    Orchestrate full post-assembly analysis.
    Returns PlasmidAnalysis dataclass.
    """
    component_summary = _extract_all_canonical_ids(resolved_modules)
    cassettes = _group_into_cassettes(resolved_modules)
    inferred_design_type = _infer_design_type(resolved_modules, component_summary)

    # RE site scanning
    assembly_strategy = design_spec.get("assembly_strategy", "")
    cloning_features: List[CloningFeature] = []

    if assembled_sequence:
        seq_upper = assembled_sequence.upper()
        seq_len = len(seq_upper)
        n_count = seq_upper.count("N")
        if seq_len == 0 or n_count / seq_len <= 0.80:
            cloning_features = scan_restriction_sites(assembled_sequence, assembly_strategy)

    conflicts = _detect_conflicts(
        resolved_modules, design_spec, component_summary, cloning_features
    )

    # LLM purpose summary
    user_request = ""
    for msg in session_messages:
        if msg.get("role") == "user":
            user_request = msg.get("content", "")
            break

    llm_result = await _llm_purpose_summary(user_request, component_summary, inferred_design_type, resolved_modules)

    inferred_purpose = (
        llm_result.get("inferred_purpose")
        or _fallback_purpose(inferred_design_type, component_summary)
    )
    llm_intent_note = llm_result.get("intent_mismatch")

    return PlasmidAnalysis(
        inferred_purpose=inferred_purpose,
        inferred_design_type=inferred_design_type,
        expression_cassettes=cassettes,
        cloning_features=cloning_features,
        component_summary=component_summary,
        conflicts=conflicts,
        llm_intent_note=llm_intent_note,
    )


# ──────────────────────────────────────────────────────────────────────────────
# Report formatting
# ──────────────────────────────────────────────────────────────────────────────

def format_analysis_report(analysis: PlasmidAnalysis) -> str:
    """
    Format PlasmidAnalysis as a markdown section appended to the workflow reply.
    Shows a compact clean line when no conflicts / notable features exist.
    """
    has_conflicts = bool(analysis.conflicts)
    has_cloning_notes = any(cf.note for cf in analysis.cloning_features)
    has_note = bool(analysis.llm_intent_note)

    lines = ["\n\n---", "### Plasmid Analysis", ""]
    lines.append(f"**Purpose**: {analysis.inferred_purpose}")

    if not has_conflicts and not has_cloning_notes and not has_note:
        lines.append("")
        lines.append("**Design Validation:** ✓ All biological sense checks passed.")
        return "\n".join(lines)

    # Expression cassettes
    cassettes = [c for c in analysis.expression_cassettes if c.cassette_type != "source_plasmid"]
    if cassettes:
        lines.append("")
        lines.append("**Expression cassettes**:")
        for i, c in enumerate(cassettes, 1):
            completeness = "complete" if c.is_complete else "incomplete"
            pol = "Pol II" if c.cassette_type == "pol2" else "Pol III"
            lines.append(f"  {i}. {pol} ({completeness}): {c.description}")

    # Cloning features — only show notable flags (e.g. Type IIs), not the full RE inventory
    if has_cloning_notes:
        notable = [cf for cf in analysis.cloning_features if cf.note]
        if notable:
            lines.append("")
            for cf in notable:
                lines.append(f"  ℹ {cf.note}")

    # Conflicts
    if has_conflicts:
        lines.append("")
        lines.append("**⚠ Conflicts:**")
        for c in analysis.conflicts:
            icon = "✗" if c.severity == "ERROR" else ("⚠" if c.severity == "WARN" else "ℹ")
            lines.append(f"  {icon} `{c.rule_id}` [{c.severity}]: {c.message}")

    # LLM intent note
    if has_note:
        lines.append("")
        lines.append(f"**Intent note**: {analysis.llm_intent_note}")

    return "\n".join(lines)

# Module Semantics System Prompt
MODULE_SEMANTICS_SYSTEM_PROMPT = """Analyze plasmid with module understanding.

CDS & Proteins:
- 2A peptides (P2A/T2A/E2A/F2A) produce SEPARATE proteins
- IRES produces SEPARATE proteins  
- No separator = ONE fusion protein

Expression:
- Pol II (CMV/CAG/EF1a): Protein-coding
- Pol III (U6/H1): Guide RNA (NO protein)

Viral Packaging:
- Lentiviral: Between LTRs = packaged
- AAV: Between ITRs = packaged

Selection:
- Bacterial: AmpR/KanR (E. coli only)
- Mammalian: PuroR/NeoR/HygR

CRISPR:
- Cas9: Pol II expression (protein)
- gRNA: Pol III (RNA, not protein)
"""


def _build_module_hierarchy_summary(resolved_modules):
    """Build module summary for LLM."""
    cats = {"Pol II (Protein)": [], "Pol III (RNA)": [], "Selection": [],
            "Viral": [], "Backbone": [], "Other": []}
    for m in resolved_modules:
        r = m.get('role', '').lower()
        n = m.get('canonical_id', m.get('name', '?'))
        if 'pol3' in r or any(x in r for x in ['u6','h1','guide']):
            cats["Pol III (RNA)"].append(n)
        elif 'pol2' in r or 'promoter' in r:
            cats["Selection" if any(x in r for x in ['selection','marker']) else "Pol II (Protein)"].append(n)
        elif any(x in r for x in ['ltr','itr','lentiviral','aav']):
            cats["Viral"].append(n)
        elif 'ori' in r or 'backbone' in r:
            cats["Backbone"].append(n)
        else:
            cats["Other"].append(n)
    return "\n".join(f"{k}: {', '.join(v)}" for k, v in cats.items() if v) or "No modules"


