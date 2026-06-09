"""
Module-aware interaction graph builder.

**Design contract**: interactions are a *readout* of the module hierarchy, not
an independent inference. Every interaction is scoped to a single detected
module (rule-based module, mammalian Pol II cassette, or CDS ORF with
submodules). We do not chain across cassette boundaries, do not invent
relationships from spatial proximity on the flat feature list, and do not
emit anything for regions that aren't enclosed by a detected module.

Input:
    - rule_based_modules:        RuleBasedModuleDetector.detect_modules output
    - mammalian_pol2_cassettes:  list[MammalianPol2Cassette.to_dict()]
    - cds_submodules:            resolve_cds_submodules() output (protein,
                                 nls, tag, linker, gap submodule dicts)

Output: a list of Interaction dicts with the same shape the frontend expects:

    {
        "interaction_id": "...",
        "interaction_type": "genetic_production" | "translation" |
                             "transcription" | "inhibition" | "conversion" |
                             "cleavage" | ...,
        "sbo_term": <SBO URI>,
        "participants": [
            {"name": ..., "start": ..., "end": ..., "strand": ...,
             "role": "stimulator"|"template"|"modifier"|"inhibitor"|...,
             "sbo_role": <SBO URI>, "so_role": <SO URI>,
             "parent_module": <module_type that owns this participant>},
            ...
        ],
        "rule_id": "INT-*-NN",
        "confidence": 0.0–1.0,
        "notes": "...",
        "source_module": "<module_type>",   # which module produced it
    }
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from .so_sbo_mapping import (
    SO_PROMOTER, SO_CDS, SO_POLYA_SIGNAL, SO_TERMINATOR, SO_OPERATOR,
    SO_RBS, SO_KOZAK, SO_MCS, SO_RECOMBINATION_SITE, SO_INSULATOR,
    SO_EPITOPE_TAG, SO_NLS, SO_FLEXIBLE_LINKER, SO_POLYPEPTIDE_REGION,
    SO_RIBOSOMAL_SKIP, SO_ENGINEERED_REGION,
    SO_INVERTED_REPEAT, SO_ORIGIN_OF_REPLICATION, SO_ORIGIN_OF_TRANSFER,
    SO_ENHANCER, SO_REGULATORY_REGION,
    SBO_GENETIC_PRODUCTION, SBO_TRANSCRIPTION, SBO_TRANSLATION,
    SBO_INHIBITION, SBO_STIMULATION, SBO_CONVERSION, SBO_CLEAVAGE,
    SBO_NON_COVALENT_BINDING,
    SBO_STIMULATOR, SBO_INHIBITOR, SBO_TEMPLATE, SBO_MODIFIER, SBO_REACTANT,
    so_role_for_module_type, so_role_for_feature_type,
)

# --------------------------------------------------------------------------- #
# Participant-builder helpers
# --------------------------------------------------------------------------- #


def _participant(
    part: Dict[str, Any],
    *,
    role: str,
    sbo_role: str,
    so_role: Optional[str] = None,
    parent_module: Optional[str] = None,
) -> Dict[str, Any]:
    """Assemble a participant dict from a module/submodule/feature-ish source."""
    resolved_so = (
        so_role
        or part.get("so_role")
        or so_role_for_module_type(part.get("module_type"))
        or so_role_for_feature_type(part.get("type"))
    )
    return {
        "name": part.get("name") or part.get("module_type") or "?",
        "start": part.get("start"),
        "end": part.get("end"),
        "strand": part.get("strand") or part.get("direction") or 1,
        "role": role,
        "sbo_role": sbo_role,
        "so_role": resolved_so,
        "parent_module": parent_module,
    }


# --------------------------------------------------------------------------- #
# Pol II / lentiviral expression cassette (authoritative source)
# --------------------------------------------------------------------------- #


_GENERIC_MODULE_TYPES = {
    "protein_module", "protein_submodule", "linker_module",
    "flexible_linker_module", "tag_module", "nls_module", "gap_module",
    "cds_module",
}


def _submodule_name(sub: Dict[str, Any]) -> Optional[str]:
    """Best-effort *real* name lookup. Returns None when only the generic
    module_type is available — callers decide whether to fall through to
    feature-row overlap."""
    return (
        sub.get("name")
        or (sub.get("metadata") or {}).get("name")
        or (sub.get("metadata") or {}).get("protein_name")
        or (sub.get("metadata") or {}).get("feature_name")
        or None
    )


def _resolve_cds_name(
    cm: Dict[str, Any],
    cds_submodules: Optional[List[Dict[str, Any]]],
    feature_rows: Optional[List[Dict[str, Any]]] = None,
) -> str:
    """Find the best display name for a Pol II cassette's CDS.

    Priority order:
      1. cm['name']                              (if the Pol2 detector set one)
      2. Join of protein_module submodules that fall within the cds_module
         window, with 2A/linker names between them (e.g. "Cas9-P2A-PuroR").
      3. pLannotate feature rows of type CDS whose span overlaps the window
         — pick the one with the largest overlap.
      4. Literal "CDS".
    """
    if cm.get("name") and cm["name"].lower() not in ("cds", "cds_module"):
        return cm["name"]

    s, e = cm.get("start"), cm.get("end")
    if s is None or e is None:
        return "CDS"

    # 2) Compose from protein_module submodules
    if cds_submodules:
        proteins = sorted(
            [
                sub for sub in cds_submodules
                if sub.get("module_type") in ("protein_module", "protein_submodule")
                and sub.get("start") is not None and sub.get("end") is not None
                and sub["start"] >= s - 5 and sub["end"] <= e + 5
            ],
            key=lambda p: p.get("start", 0),
        )
        # Only compose from submodules when at least one protein submodule
        # carries a real name. Otherwise the composition would be a string of
        # module_type labels ("protein_module-flexible_linker_module-...").
        # In that case, fall through to feature-row overlap lookup below.
        protein_names = [_submodule_name(p) for p in proteins]
        has_real_names = any(n for n in protein_names)
        if proteins and has_real_names:
            linkers = [
                sub for sub in cds_submodules
                if sub.get("module_type") in ("linker_module", "flexible_linker_module")
                and sub.get("start") is not None and sub.get("end") is not None
                and sub["start"] >= s and sub["end"] <= e
            ]
            parts: List[str] = []
            for i, p in enumerate(proteins):
                if i > 0:
                    prev_end = proteins[i - 1]["end"]
                    cur_start = p["start"]
                    best_linker = None
                    for l in linkers:
                        if l["start"] >= prev_end and l["end"] <= cur_start:
                            if best_linker is None or l["start"] < best_linker["start"]:
                                best_linker = l
                    if best_linker:
                        ln = _submodule_name(best_linker)
                        if ln:
                            parts.append(ln)
                pn = _submodule_name(p)
                if pn:
                    parts.append(pn)
            joined = "-".join(parts)
            if joined:
                return joined

    # 3) Overlap against pLannotate feature rows of type CDS
    if feature_rows:
        best_name = None
        best_overlap = 0
        for row in feature_rows:
            t = (row.get("type") or "").lower()
            if t not in ("cds", "gene", "protein_generator"):
                continue
            rs = row.get("start")
            re_ = row.get("end")
            nm = row.get("name")
            if rs is None or re_ is None or not nm:
                continue
            overlap = max(0, min(e, re_) - max(s, rs))
            if overlap > best_overlap:
                best_overlap = overlap
                best_name = nm
        if best_name:
            return best_name

    return cm.get("name") or "CDS"


def _pol2_cassette_interactions(
    cas: Dict[str, Any],
    idx: int,
    cds_submodules: Optional[List[Dict[str, Any]]] = None,
    feature_rows: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """For each Pol II cassette dict (with upstream_regulatory / cds_module /
    downstream_regulatory), emit exactly one genetic_production interaction."""
    ur = cas.get("upstream_regulatory") or {}
    cm = cas.get("cds_module") or {}
    dr = cas.get("downstream_regulatory") or {}

    # Resolve promoter *name* — to_dict() emits:
    #   primary_promoter: str (name) OR None
    #   components:       List[str] (names) of every upstream feature
    # Coordinates always come from the upstream_regulatory span itself.
    promoter_name = None
    primary = ur.get("primary_promoter")
    if primary and isinstance(primary, str):
        promoter_name = primary
    else:
        for comp in ur.get("components") or []:
            nm = comp if isinstance(comp, str) else comp.get("name", "")
            if not nm:
                continue
            low = nm.lower()
            if "promoter" in low or "enhancer" in low:
                promoter_name = nm
                break
    if not promoter_name and ur.get("components"):
        # Fall back to first component name
        first = ur["components"][0]
        promoter_name = first if isinstance(first, str) else first.get("name")

    if ur.get("start") is None or cm.get("start") is None:
        return []

    promoter = {
        "name": promoter_name or "upstream regulatory",
        "start": ur.get("start"),
        "end": ur.get("end"),
        "strand": ur.get("strand", 1),
    }

    cds_name = _resolve_cds_name(cm, cds_submodules, feature_rows=feature_rows)
    cds_participant = {
        "name": cds_name,
        "start": cm.get("start"),
        "end": cm.get("end"),
        "strand": cm.get("strand", 1),
    }

    source_type = cas.get("module_type", "pol2_cassette")
    participants = [
        _participant(promoter, role="stimulator",
                     sbo_role=SBO_STIMULATOR, so_role=SO_PROMOTER,
                     parent_module=source_type),
        _participant(cds_participant,
                     role="template", sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                     parent_module=source_type),
    ]
    polya = None
    if dr and dr.get("start") is not None:
        polya_name = dr.get("polya")
        if not polya_name and dr.get("components"):
            # Prefer a polyA/terminator-looking component
            for comp in dr["components"]:
                nm = comp if isinstance(comp, str) else comp.get("name", "")
                low = (nm or "").lower()
                if "poly" in low or "terminator" in low or "pa" in low:
                    polya_name = nm
                    break
            if not polya_name:
                first = dr["components"][0]
                polya_name = first if isinstance(first, str) else first.get("name")
        polya = {
            "name": polya_name or "polyA / terminator",
            "start": dr.get("start"),
            "end": dr.get("end"),
            "strand": dr.get("strand", 1),
        }
        participants.append(
            _participant(polya, role="modifier",
                         sbo_role=SBO_MODIFIER, so_role=SO_POLYA_SIGNAL,
                         parent_module=source_type)
        )

    # Emit two peer interactions (upstream drives CDS, downstream modifies CDS)
    # instead of the single wrapper-level cassette interaction (INT-POL2-CAS-01).
    ixs: List[Dict[str, Any]] = []

    # INT-POL2-UR-01 — upstream_regulatory_module stimulates the CDS
    ixs.append({
        "interaction_id": f"pol2_ur_{idx}",
        "interaction_type": "stimulation",
        "sbo_term": SBO_STIMULATION,
        "participants": [
            _participant(promoter, role="stimulator",
                         sbo_role=SBO_STIMULATOR, so_role=SO_PROMOTER,
                         parent_module="upstream_regulatory_module"),
            _participant(cds_participant, role="template",
                         sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                         parent_module="cds_module"),
        ],
        "rule_id": "INT-POL2-UR-01",
        "confidence": cas.get("weight", 0.9),
        "notes": (
            f"{promoter['name']} (upstream regulatory) drives transcription of "
            f"{cds_name}."
        ),
        "source_module": "upstream_regulatory_module",
    })

    # INT-POL2-DR-01 — downstream_regulatory_module modifies / terminates CDS
    if polya is not None:
        ixs.append({
            "interaction_id": f"pol2_dr_{idx}",
            "interaction_type": "stimulation",
            "sbo_term": SBO_STIMULATION,
            "participants": [
                _participant(polya, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_POLYA_SIGNAL,
                             parent_module="downstream_regulatory_module"),
                _participant(cds_participant, role="template",
                             sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                             parent_module="cds_module"),
            ],
            "rule_id": "INT-POL2-DR-01",
            "confidence": cas.get("weight", 0.9),
            "notes": (
                f"{polya['name']} (downstream regulatory) terminates "
                f"transcription and stabilizes the {cds_name} mRNA."
            ),
            "source_module": "downstream_regulatory_module",
        })

    return ixs


# --------------------------------------------------------------------------- #
# lacZα blue/white screening module — rich submodule graph
# --------------------------------------------------------------------------- #


def _lac_bw_interactions(mod: Dict[str, Any], idx: int) -> List[Dict[str, Any]]:
    """Emit the within-module relationships declared by a
    lac_alpha_blue_white_module's submodule list."""
    subs_by_type: Dict[str, List[Dict[str, Any]]] = {}
    for s in mod.get("submodules") or []:
        subs_by_type.setdefault(s.get("module_type", ""), []).append(s)

    alpha = (subs_by_type.get("lac_alpha_cds") or [None])[0]
    promoter = (subs_by_type.get("lac_promoter") or [None])[0]
    operator = (subs_by_type.get("lac_operator") or [None])[0]
    mcs = (subs_by_type.get("mcs") or [None])[0]
    laci = (subs_by_type.get("lac_i_gene") or [None])[0]

    source_type = mod.get("module_type", "lac_alpha_blue_white_module")
    ixs: List[Dict[str, Any]] = []

    # lac promoter transcribes lacZα
    if promoter and alpha:
        ixs.append({
            "interaction_id": f"lac_bw_{idx}_txn",
            "interaction_type": "genetic_production",
            "sbo_term": SBO_GENETIC_PRODUCTION,
            "participants": [
                _participant(promoter, role="stimulator",
                             sbo_role=SBO_STIMULATOR, so_role=SO_PROMOTER,
                             parent_module=source_type),
                _participant(alpha, role="template",
                             sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                             parent_module=source_type),
            ],
            "rule_id": "INT-LACBW-01",
            "confidence": mod.get("weight", 0.9),
            "notes": "lac promoter drives transcription of lacZα fragment "
                     "(blue on X-gal/IPTG when intact).",
            "source_module": source_type,
        })

    # Operator inhibits promoter (LacI-mediated repression)
    if operator and promoter:
        ixs.append({
            "interaction_id": f"lac_bw_{idx}_rep",
            "interaction_type": "inhibition",
            "sbo_term": SBO_INHIBITION,
            "participants": [
                _participant(operator, role="inhibitor",
                             sbo_role=SBO_INHIBITOR, so_role=SO_OPERATOR,
                             parent_module=source_type),
                _participant(promoter, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_PROMOTER,
                             parent_module=source_type),
            ],
            "rule_id": "INT-LACBW-02",
            "confidence": mod.get("weight", 0.9),
            "notes": "lac operator is bound by LacI in absence of IPTG, "
                     "repressing the lac promoter.",
            "source_module": source_type,
        })

    # LacI binds operator (non-covalent binding)
    if laci and operator:
        ixs.append({
            "interaction_id": f"lac_bw_{idx}_bind",
            "interaction_type": "non_covalent_binding",
            "sbo_term": SBO_NON_COVALENT_BINDING,
            "participants": [
                _participant(laci, role="stimulator",
                             sbo_role=SBO_STIMULATOR, so_role=SO_CDS,
                             parent_module=source_type),
                _participant(operator, role="template",
                             sbo_role=SBO_TEMPLATE, so_role=SO_OPERATOR,
                             parent_module=source_type),
            ],
            "rule_id": "INT-LACBW-03",
            "confidence": mod.get("weight", 0.9),
            "notes": "LacI repressor binds lac operator — IPTG relieves binding.",
            "source_module": source_type,
        })

    # Insertion into MCS disrupts lacZα (blue/white gating)
    if mcs and alpha:
        ixs.append({
            "interaction_id": f"lac_bw_{idx}_disrupt",
            "interaction_type": "inhibition",
            "sbo_term": SBO_INHIBITION,
            "participants": [
                _participant(mcs, role="inhibitor",
                             sbo_role=SBO_INHIBITOR, so_role=SO_MCS,
                             parent_module=source_type),
                _participant(alpha, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_CDS,
                             parent_module=source_type),
            ],
            "rule_id": "INT-LACBW-04",
            "confidence": mod.get("weight", 0.9),
            "notes": "MCS insertion disrupts lacZα reading frame → white colony on X-gal/IPTG.",
            "source_module": source_type,
        })

    # INT-LACBW-04b: post-cloning state — an inserted CDS overlaps the
    # lacZα ORF (≥30 bp same strand). The detector populates
    # mod["inserted_cds"] when a non-MCS CDS overlaps the α window.
    inserted_cds = None
    a_strand = (alpha or {}).get("strand", 1) if alpha else 1
    for ic in mod.get("inserted_cds") or []:
        if not alpha:
            break
        a_s, a_e = int(alpha.get("start") or 0), int(alpha.get("end") or 0)
        c_s, c_e = int(ic.get("start") or 0), int(ic.get("end") or 0)
        if c_e <= a_s or a_e <= c_s:
            continue
        ov = min(c_e, a_e) - max(c_s, a_s)
        if ov >= 30 and (ic.get("strand") or 1) == a_strand:
            inserted_cds = ic
            break
    if inserted_cds and alpha:
        ixs.append({
            "interaction_id": f"lac_bw_{idx}_disrupt_post",
            "interaction_type": "inhibition",
            "sbo_term": SBO_INHIBITION,
            "participants": [
                _participant(inserted_cds, role="inhibitor",
                             sbo_role=SBO_INHIBITOR, so_role=SO_CDS,
                             parent_module=source_type),
                _participant(alpha, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_CDS,
                             parent_module=source_type),
            ],
            "rule_id": "INT-LACBW-04b",
            "confidence": mod.get("weight", 0.9),
            "notes": "Cloned insert disrupts lacZα reading frame "
                     "(post-MCS-insertion white-colony state).",
            "source_module": source_type,
        })

    return ixs


# --------------------------------------------------------------------------- #
# Recombination-flanked modules (Cre/FLP/Gateway/integrase)
# --------------------------------------------------------------------------- #

_RECOMB_MAP: Dict[str, Dict[str, str]] = {
    "floxed_region":            {"flank": "loxP",     "enzyme": "Cre recombinase",   "rule": "INT-REC-CRE-01"},
    "lsl_cassette":             {"flank": "loxP",     "enzyme": "Cre recombinase",   "rule": "INT-REC-CRE-02"},
    "floxed_cassette":          {"flank": "loxP",     "enzyme": "Cre recombinase",   "rule": "INT-REC-CRE-01"},
    "frt_flanked_cassette":     {"flank": "FRT",      "enzyme": "Flp recombinase",   "rule": "INT-REC-FLP-01"},
    "gateway_entry_cassette":   {"flank": "attL",     "enzyme": "LR clonase",        "rule": "INT-REC-GW-01"},
    "gateway_dest_cassette":    {"flank": "attR",     "enzyme": "LR clonase",        "rule": "INT-REC-GW-02"},
    "gateway_recombination":    {"flank": "att pair", "enzyme": "BP/LR clonase",     "rule": "INT-REC-GW-03"},
    "integrase_landing_pad":    {"flank": "attP/attB","enzyme": "phage integrase",   "rule": "INT-REC-INT-01"},
}


def _recombination_interaction(mod: Dict[str, Any], idx: int) -> List[Dict[str, Any]]:
    spec = _RECOMB_MAP.get(mod.get("module_type", ""))
    if not spec:
        return []

    source_type = mod.get("module_type", "recombination_cassette")
    # Try to use the module's `features`/flank metadata when available;
    # otherwise present the module as a single "cassette" participant.
    participants: List[Dict[str, Any]] = []

    # Prefer labeled flanking-feature names + spans when the detector
    # supplied them; fall back to the legacy generic 5'/3' synthesis.
    flanks = mod.get("flanking_features") or []
    if flanks:
        for f in flanks:
            participants.append({
                "name": f.get("name") or spec["flank"],
                "start": f.get("start"),
                "end": f.get("end"),
                "strand": f.get("strand", mod.get("strand", 1)),
                "role": "reactant",
                "sbo_role": SBO_STIMULATOR,
                "so_role": SO_RECOMBINATION_SITE,
                "parent_module": source_type,
            })
    else:
        flank1_name = mod.get("loxp_site_1") or mod.get("frt_site_1") or f"5' {spec['flank']}"
        flank2_name = mod.get("loxp_site_2") or mod.get("frt_site_2") or f"3' {spec['flank']}"
        participants.append({
            "name": flank1_name,
            "start": mod.get("start"), "end": mod.get("start"),
            "strand": mod.get("strand", 1),
            "role": "reactant",
            "sbo_role": SBO_STIMULATOR,
            "so_role": SO_RECOMBINATION_SITE,
            "parent_module": source_type,
        })
        participants.append({
            "name": flank2_name,
            "start": mod.get("end"), "end": mod.get("end"),
            "strand": mod.get("strand", 1),
            "role": "reactant",
            "sbo_role": SBO_STIMULATOR,
            "so_role": SO_RECOMBINATION_SITE,
            "parent_module": source_type,
        })
    participants.append({
        "name": mod.get("name") or source_type,
        "start": mod.get("start"), "end": mod.get("end"),
        "strand": mod.get("strand", 1),
        "role": "modifier",
        "sbo_role": SBO_MODIFIER,
        "so_role": so_role_for_module_type(source_type),
        "parent_module": source_type,
    })
    participants.append({
        "name": spec["enzyme"],
        "start": None, "end": None, "strand": None,
        "role": "stimulator",
        "sbo_role": SBO_STIMULATOR,
        "so_role": None,
        "parent_module": source_type,
        "external": True,
    })

    return [{
        "interaction_id": f"recomb_{idx}",
        "interaction_type": "conversion",
        "sbo_term": SBO_CONVERSION,
        "participants": participants,
        "rule_id": spec["rule"],
        "confidence": mod.get("weight", 0.9),
        "notes": (
            f"{spec['flank']}-flanked cassette — {spec['enzyme']} mediates "
            "excision / inversion / integration of the internal region."
        ),
        "source_module": source_type,
    }]


# --------------------------------------------------------------------------- #
# Insulated expression block — boundary stimulation
# --------------------------------------------------------------------------- #


def _insulator_block_interaction(
    mod: Dict[str, Any], idx: int
) -> List[Dict[str, Any]]:
    source_type = mod.get("module_type", "insulated_expression_block")
    return [{
        "interaction_id": f"insulator_block_{idx}",
        "interaction_type": "stimulation",
        "sbo_term": SBO_STIMULATION,
        "participants": [
            {
                "name": mod.get("name", "paired insulator boundary"),
                "start": mod.get("start"), "end": mod.get("end"),
                "strand": mod.get("strand", 1),
                "role": "modifier",
                "sbo_role": SBO_MODIFIER, "so_role": SO_INSULATOR,
                "parent_module": source_type,
            },
        ],
        "rule_id": "INT-INS-01",
        "confidence": mod.get("weight", 0.9),
        "notes": "Paired insulators bracket an enhancer-blocked expression boundary.",
        "source_module": source_type,
    }]


# --------------------------------------------------------------------------- #
# CDS fusion interactions (P2A/T2A self-cleavage; tag-protein fusion)
# --------------------------------------------------------------------------- #


def _is_2a_linker(sub: Dict[str, Any]) -> bool:
    """True if this linker submodule is a ribosomal-skip 2A peptide.

    Checks `name`, `payload_id`, and metadata for P2A/T2A/E2A/F2A markers.
    Accepts either `linker_module` or `flexible_linker_module` as the
    parent type (the module_extractor emits the latter for gap-detected 2As).
    """
    mt = sub.get("module_type") or ""
    if mt not in ("linker_module", "flexible_linker_module"):
        return False
    fields = [
        (sub.get("name") or ""),
        (sub.get("payload_id") or ""),
        (sub.get("metadata") or {}).get("detected_linker", "") or "",
        (sub.get("metadata") or {}).get("name", "") or "",
    ]
    blob = " ".join(str(f).lower() for f in fields)
    return any(tag in blob for tag in ("p2a", "t2a", "e2a", "f2a"))


def _2a_name(sub: Dict[str, Any]) -> str:
    """Canonical 2A peptide label for description text."""
    for field in ("name", "payload_id"):
        v = sub.get(field)
        if v:
            up = str(v).upper()
            for tag in ("P2A", "T2A", "E2A", "F2A"):
                if tag in up:
                    return tag
    return "2A peptide"


def _find_containing_pol2_cassette(
    protein_sub: Dict[str, Any],
    mammalian_pol2_cassettes: Optional[List[Dict[str, Any]]],
) -> Optional[Dict[str, Any]]:
    """Return the Pol II cassette whose cds_module span encloses the given
    protein submodule, or None if no cassette claims it."""
    if not mammalian_pol2_cassettes or protein_sub.get("start") is None:
        return None
    ps, pe = int(protein_sub["start"]), int(protein_sub["end"])
    pstrand = int(protein_sub.get("strand") or 1)
    for cas in mammalian_pol2_cassettes:
        cm = cas.get("cds_module") or {}
        if cm.get("start") is None or cm.get("end") is None:
            continue
        if int(cm.get("strand", 1)) != pstrand:
            continue
        if int(cm["start"]) <= ps and int(cm["end"]) >= pe:
            return cas
    return None


def _cassette_promoter(cas: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Extract a display promoter dict from a Pol II cassette."""
    ur = cas.get("upstream_regulatory") or {}
    primary = ur.get("primary_promoter")
    name = None
    if isinstance(primary, str) and primary:
        name = primary
    elif ur.get("components"):
        for c in ur["components"]:
            nm = c if isinstance(c, str) else c.get("name", "")
            low = (nm or "").lower()
            if "promoter" in low or "enhancer" in low:
                name = nm
                break
        if not name:
            first = ur["components"][0]
            name = first if isinstance(first, str) else first.get("name")
    if ur.get("start") is None:
        return None
    return {
        "name": name or "upstream regulatory",
        "start": ur.get("start"),
        "end": ur.get("end"),
        "strand": ur.get("strand", cas.get("strand", 1)),
    }


def _cds_fusion_interactions(
    cds_submodules: List[Dict[str, Any]],
    mammalian_pol2_cassettes: Optional[List[Dict[str, Any]]] = None,
    feature_rows: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """Emit interactions capturing 2A ribosomal-skip biology.

    For each P2A/T2A/E2A/F2A linker bridging two protein submodules inside
    the same CDS:
      1. Emit one `cleavage` interaction (the 2A peptide is the stimulator;
         the flanking proteins are products). Communicates that a single
         polypeptide is co-translationally split into two independent
         polypeptides.
      2. If the CDS is enclosed by a Pol II cassette, emit two additional
         `genetic_production` interactions (rule INT-CDS-2A-02): one from
         the cassette's promoter to each resulting protein product. This
         makes explicit that *both* proteins are driven by the same upstream
         regulatory elements.

    Submodules carry their own strand/start/end; intervening tag/NLS
    submodules between the protein and the 2A linker are skipped.
    """
    if not cds_submodules:
        return []

    # Sort by start coordinate on each strand; within that order, adjacent
    # (protein | 2A linker | protein) triples describe a ribosomal skip.
    sorted_subs = sorted(
        (s for s in cds_submodules if s.get("start") is not None),
        key=lambda s: (s.get("strand", 1), s.get("start", 0))
    )

    ixs: List[Dict[str, Any]] = []
    protein_types = ("protein_module", "protein_submodule")
    for i, cur_sub in enumerate(sorted_subs):
        if not _is_2a_linker(cur_sub):
            continue
        # Scan outward for the nearest protein submodule on each side
        # (intervening tag/NLS/linker submodules are allowed).
        prev_sub = None
        for j in range(i - 1, -1, -1):
            if sorted_subs[j].get("strand") != cur_sub.get("strand"):
                break
            if sorted_subs[j].get("module_type") in protein_types:
                prev_sub = sorted_subs[j]
                break
        next_sub = None
        for j in range(i + 1, len(sorted_subs)):
            if sorted_subs[j].get("strand") != cur_sub.get("strand"):
                break
            if sorted_subs[j].get("module_type") in protein_types:
                next_sub = sorted_subs[j]
                break
        if prev_sub is None or next_sub is None:
            continue

        # Resolve protein names: submodule payload_id/name first, feature-row overlap fallback
        def _protein_name(p):
            nm = _submodule_name(p)
            if nm:
                return nm
            if feature_rows:
                ps, pe = p.get("start"), p.get("end")
                if ps is not None and pe is not None:
                    best = None
                    best_ov = 0
                    for row in feature_rows:
                        rt = (row.get("type") or "").lower()
                        if rt not in ("cds", "gene", "protein_generator"):
                            continue
                        rs, re_, rn = row.get("start"), row.get("end"), row.get("name")
                        if rs is None or re_ is None or not rn:
                            continue
                        ov = max(0, min(pe, re_) - max(ps, rs))
                        if ov > best_ov:
                            best_ov = ov
                            best = rn
                    if best:
                        return best
            return "protein"

        prev_name = _protein_name(prev_sub)
        next_name = _protein_name(next_sub)
        tag_name = _2a_name(cur_sub)

        enclosing_cas = (_find_containing_pol2_cassette(prev_sub, mammalian_pol2_cassettes)
                         or _find_containing_pol2_cassette(next_sub, mammalian_pol2_cassettes))
        source_module = (enclosing_cas.get("module_type")
                         if enclosing_cas else "cds_orf")

        # 1) Cleavage interaction — intrinsic 2A biology
        ixs.append({
            "interaction_id": f"p2a_cleavage_{i}",
            "interaction_type": "cleavage",
            "sbo_term": SBO_CLEAVAGE,
            "participants": [
                _participant({**cur_sub, "name": tag_name},
                             role="stimulator",
                             sbo_role=SBO_STIMULATOR, so_role=SO_RIBOSOMAL_SKIP,
                             parent_module=source_module),
                _participant({**prev_sub, "name": prev_name}, role="product",
                             sbo_role=SBO_MODIFIER, so_role=SO_POLYPEPTIDE_REGION,
                             parent_module=source_module),
                _participant({**next_sub, "name": next_name}, role="product",
                             sbo_role=SBO_MODIFIER, so_role=SO_POLYPEPTIDE_REGION,
                             parent_module=source_module),
            ],
            "rule_id": "INT-CDS-2A-01",
            "confidence": 0.95,
            "notes": (
                f"{tag_name} ribosomal-skip peptide — co-translational cleavage produces "
                f"independent {prev_name} and {next_name} polypeptides from a single ORF."
            ),
            "source_module": source_module,
        })

        # 2) Cassette-coupled genetic_production — one promoter drives both proteins
        if enclosing_cas is not None:
            promoter = _cassette_promoter(enclosing_cas)
            dr = enclosing_cas.get("downstream_regulatory") or {}
            polya_label = dr.get("polya")
            if not polya_label and dr.get("components"):
                first = dr["components"][0]
                polya_label = first if isinstance(first, str) else (first or {}).get("name", "polyA")
            if promoter is not None:
                for k, (prot, prot_name) in enumerate(
                    [(prev_sub, prev_name), (next_sub, next_name)]
                ):
                    participants = [
                        _participant(promoter, role="stimulator",
                                     sbo_role=SBO_STIMULATOR, so_role=SO_PROMOTER,
                                     parent_module=source_module),
                        _participant({**prot, "name": prot_name},
                                     role="template",
                                     sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                                     parent_module=source_module),
                    ]
                    if dr.get("start") is not None:
                        participants.append(_participant(
                            {
                                "name": polya_label or "polyA",
                                "start": dr.get("start"),
                                "end": dr.get("end"),
                                "strand": dr.get("strand", 1),
                            },
                            role="modifier",
                            sbo_role=SBO_MODIFIER, so_role=SO_POLYA_SIGNAL,
                            parent_module=source_module,
                        ))
                    ixs.append({
                        "interaction_id": f"p2a_coexpr_{i}_{k}",
                        "interaction_type": "genetic_production",
                        "sbo_term": SBO_GENETIC_PRODUCTION,
                        "participants": participants,
                        "rule_id": "INT-CDS-2A-02",
                        "confidence": 0.9,
                        "notes": (
                            f"{promoter['name']} drives expression of {prot_name} — "
                            f"one of two independent proteins released by the {tag_name} "
                            f"ribosomal skip; both share the same upstream regulatory elements."
                        ),
                        "source_module": source_module,
                    })

    return ixs




def _lentiviral_three_module_interactions(
    rule_based_modules: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Emit the pair of interactions that tie the three lentiviral peer
    modules into a functional graph:
        upstream_regulatory  — STIMULATES  → payload
        downstream_regulatory — MODIFIES   → payload

    Each LTR-bounded payload gets its own pair. If multiple payloads are
    detected (rare), we pair each with the nearest regulatory modules on
    each side by position.
    """
    payloads = [m for m in rule_based_modules if m.get("module_type") == "lentiviral_payload"]
    uprs = [m for m in rule_based_modules if m.get("module_type") == "lentiviral_upstream_regulatory"]
    drrs = [m for m in rule_based_modules if m.get("module_type") == "lentiviral_downstream_regulatory"]
    if not payloads:
        return []

    ixs: List[Dict[str, Any]] = []
    for p_idx, payload in enumerate(payloads):
        p_start = int(payload.get("start", 0))
        p_end = int(payload.get("end", 0))

        # Upstream: the UR module ending closest to (and at/before) payload start
        nearest_ur = None
        for ur in uprs:
            if int(ur.get("end", 0)) <= p_end:  # must end before or at payload end
                if nearest_ur is None or ur["end"] > nearest_ur["end"]:
                    nearest_ur = ur
        # Downstream: the DR module starting closest to (and at/after) payload end
        nearest_dr = None
        for dr in drrs:
            if int(dr.get("start", 0)) >= p_start:
                if nearest_dr is None or dr["start"] < nearest_dr["start"]:
                    nearest_dr = dr

        if nearest_ur is not None:
            ixs.append({
                "interaction_id": f"lenti_ur_{p_idx}",
                "interaction_type": "stimulation",
                "sbo_term": SBO_STIMULATION,
                "participants": [
                    _participant(nearest_ur, role="stimulator",
                                 sbo_role=SBO_STIMULATOR, so_role=None,
                                 parent_module="lentiviral_upstream_regulatory"),
                    _participant(payload, role="template",
                                 sbo_role=SBO_TEMPLATE, so_role=None,
                                 parent_module="lentiviral_payload"),
                ],
                "rule_id": "INT-LENTI-UR-01",
                "confidence": min(nearest_ur.get("weight", 0.9), payload.get("weight", 0.9)),
                "notes": (
                    f"{nearest_ur.get('name', 'Lentiviral Upstream Regulatory')} drives "
                    f"transcription across the lentiviral payload "
                    f"(integration-competent LTR-to-LTR region)."
                ),
                "source_module": "lentiviral_payload",
            })

        if nearest_dr is not None:
            ixs.append({
                "interaction_id": f"lenti_dr_{p_idx}",
                "interaction_type": "stimulation",
                "sbo_term": SBO_STIMULATION,
                "participants": [
                    _participant(nearest_dr, role="modifier",
                                 sbo_role=SBO_MODIFIER, so_role=None,
                                 parent_module="lentiviral_downstream_regulatory"),
                    _participant(payload, role="template",
                                 sbo_role=SBO_TEMPLATE, so_role=None,
                                 parent_module="lentiviral_payload"),
                ],
                "rule_id": "INT-LENTI-DR-01",
                "confidence": min(nearest_dr.get("weight", 0.9), payload.get("weight", 0.9)),
                "notes": (
                    f"{nearest_dr.get('name', 'Lentiviral Downstream Regulatory')} "
                    f"post-transcriptionally regulates the lentiviral payload "
                    f"(WPRE-mediated stabilization and polyA termination)."
                ),
                "source_module": "lentiviral_payload",
            })

    return ixs

# --------------------------------------------------------------------------- #
# Dispatch table: module_type → interaction builder(s)
# --------------------------------------------------------------------------- #

def _guide_expression_cassette_interactions(mod: Dict[str, Any], idx: int) -> List[Dict[str, Any]]:
    """INT-POL3-01: Pol III promoter → sgRNA scaffold transcription.

    Emits one  interaction per guide_expression_cassette module,
    wiring the pol3_promoter submodule (stimulator) to the sgrna_scaffold
    submodule (template). Scoped strictly to within-module submodules — no
    cross-module spatial inference.
    """
    subs_by_type: Dict[str, List[Dict[str, Any]]] = {}
    for s in mod.get("submodules") or []:
        subs_by_type.setdefault(s.get("module_type", ""), []).append(s)

    promoter = (subs_by_type.get("pol3_promoter") or [None])[0]
    scaffold = (subs_by_type.get("sgrna_scaffold") or [None])[0]

    if not promoter or not scaffold:
        return []

    source_type = mod.get("module_type", "guide_expression_cassette")
    return [{
        "interaction_id": f"pol3_guide_{idx}_txn",
        "interaction_type": "transcription",
        "sbo_term": SBO_TRANSCRIPTION,
        "participants": [
            _participant(promoter, role="stimulator",
                         sbo_role=SBO_STIMULATOR, so_role=SO_PROMOTER,
                         parent_module=source_type),
            _participant(scaffold, role="template",
                         sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                         parent_module=source_type),
        ],
        "rule_id": "INT-POL3-01",
        "confidence": mod.get("weight", 0.9),
        "notes": "Pol III promoter drives transcription of guide-RNA scaffold "
                 "(U6/H1/7SK → sgRNA/tracrRNA/pegRNA scaffold).",
        "source_module": source_type,
    }]



# --------------------------------------------------------------------------- #
# Mobilizable replicon (ori + bom)  —  INT-MOB-*
# --------------------------------------------------------------------------- #


def _mobilizable_replicon_interactions(mod, idx):
    """INT-MOB-01: bom oriT paired with a replication origin.

    Detector populates submodule list with `ori` and `bom` entries; optionally
    a `mob_cds` (relaxase) entry triggers INT-MOB-02.
    """
    subs_by_type = {}
    for s in mod.get("submodules") or []:
        subs_by_type.setdefault(s.get("module_type", ""), []).append(s)
    ori = (subs_by_type.get("ori") or [None])[0]
    bom = (subs_by_type.get("bom") or [None])[0]
    mob = (subs_by_type.get("mob_cds") or [None])[0]
    source_type = mod.get("module_type", "mobilizable_replicon")
    ixs = []

    if ori and bom:
        ixs.append({
            "interaction_id": f"mob_{idx}_oriT",
            "interaction_type": "non_covalent_binding",
            "sbo_term": SBO_NON_COVALENT_BINDING,
            "participants": [
                _participant(bom, role="reactant",
                             sbo_role=SBO_REACTANT, so_role=SO_ORIGIN_OF_TRANSFER,
                             parent_module=source_type),
                _participant(ori, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_ORIGIN_OF_REPLICATION,
                             parent_module=source_type),
            ],
            "rule_id": "INT-MOB-01",
            "confidence": mod.get("weight", 0.95),
            "notes": "bom (basis of mobility) flanks ori — relaxase nick site "
                     "for conjugal transfer when MobA/TraI is supplied in trans.",
            "source_module": source_type,
        })

    if mob and bom:
        ixs.append({
            "interaction_id": f"mob_{idx}_relaxase",
            "interaction_type": "cleavage",
            "sbo_term": SBO_CLEAVAGE,
            "participants": [
                _participant(mob, role="stimulator",
                             sbo_role=SBO_STIMULATOR, so_role=SO_CDS,
                             parent_module=source_type),
                _participant(bom, role="template",
                             sbo_role=SBO_TEMPLATE, so_role=SO_ORIGIN_OF_TRANSFER,
                             parent_module=source_type),
            ],
            "rule_id": "INT-MOB-02",
            "confidence": mod.get("weight", 0.9),
            "notes": "Relaxase nicks bom oriT for single-strand transfer.",
            "source_module": source_type,
        })

    return ixs


# --------------------------------------------------------------------------- #
# Tn3 transposon — INT-TRANSP-TN3-*
# --------------------------------------------------------------------------- #


def _tn3_transposon_interactions(mod, idx):
    """INT-TRANSP-TN3-01: Tn3 inverted-repeat pair flanks transposable element.

    Detector populates `flanking_features` with the IR pair (left, right) and
    optionally `cargo_cds` with the enclosed bla / tnpA CDS.
    """
    flanks = mod.get("flanking_features") or []
    if len(flanks) < 2:
        return []
    ir_left, ir_right = flanks[0], flanks[1]
    cargo_cds_list = mod.get("cargo_cds") or []
    cargo_cds = cargo_cds_list[0] if cargo_cds_list else None
    source_type = mod.get("module_type", "tn3_transposon")
    ixs = []

    ixs.append({
        "interaction_id": f"tn3_{idx}_ir_pair",
        "interaction_type": "conversion",
        "sbo_term": SBO_CONVERSION,
        "participants": [
            _participant(ir_left, role="reactant",
                         sbo_role=SBO_REACTANT, so_role=SO_INVERTED_REPEAT,
                         parent_module=source_type),
            _participant(ir_right, role="reactant",
                         sbo_role=SBO_REACTANT, so_role=SO_INVERTED_REPEAT,
                         parent_module=source_type),
        ],
        "rule_id": "INT-TRANSP-TN3-01",
        "confidence": mod.get("weight", 0.9),
        "notes": "Tn3 inverted repeat pair flanks transposable element — "
                 "TnpA-mediated cointegrate / resolution.",
        "source_module": source_type,
    })

    if cargo_cds:
        ixs.append({
            "interaction_id": f"tn3_{idx}_cargo",
            "interaction_type": "stimulation",
            "sbo_term": SBO_STIMULATION,
            "participants": [
                _participant(cargo_cds, role="template",
                             sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                             parent_module=source_type),
            ],
            "rule_id": "INT-TRANSP-TN3-02",
            "confidence": mod.get("weight", 0.9),
            "notes": "Cargo CDS (bla / tnpA) enclosed by IR pair.",
            "source_module": source_type,
        })

    return ixs


# --------------------------------------------------------------------------- #
# Lac promoter regulatory unit — INT-LAC-OP-01, INT-CAP-01
# --------------------------------------------------------------------------- #


def _lac_promoter_regulatory_interactions(mod, idx):
    """INT-LAC-OP-01 + INT-CAP-01: lac promoter with operator/CAP modifiers
    *outside* a blue-white module.
    """
    subs_by_type = {}
    for s in mod.get("submodules") or []:
        subs_by_type.setdefault(s.get("module_type", ""), []).append(s)
    promoter = (subs_by_type.get("lac_promoter") or [None])[0]
    operator = (subs_by_type.get("lac_operator") or [None])[0]
    cap_site = (subs_by_type.get("cap_binding_site") or [None])[0]
    laci = (subs_by_type.get("lac_i_gene") or [None])[0]
    source_type = mod.get("module_type", "lac_promoter_regulatory_unit")
    ixs = []

    if operator and promoter:
        ixs.append({
            "interaction_id": f"lac_reg_{idx}_op",
            "interaction_type": "inhibition",
            "sbo_term": SBO_INHIBITION,
            "participants": [
                _participant(operator, role="inhibitor",
                             sbo_role=SBO_INHIBITOR, so_role=SO_OPERATOR,
                             parent_module=source_type),
                _participant(promoter, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_PROMOTER,
                             parent_module=source_type),
            ],
            "rule_id": "INT-LAC-OP-01",
            "confidence": mod.get("weight", 0.9),
            "notes": "lac operator → LacI-mediated repression of lac promoter.",
            "source_module": source_type,
        })

    # LacI binds the lac operator. Lactose (allolactose) or IPTG
    # inactivates LacI and relieves repression — hence the
    # "lactose-conditional" inhibition of the operator.
    if laci and operator:
        ixs.append({
            "interaction_id": f"lac_reg_{idx}_laci",
            "interaction_type": "inhibition",
            "sbo_term": SBO_INHIBITION,
            "participants": [
                _participant(laci, role="inhibitor",
                             sbo_role=SBO_INHIBITOR, so_role=SO_CDS,
                             parent_module=source_type),
                _participant(operator, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_OPERATOR,
                             parent_module=source_type),
            ],
            "rule_id": "INT-LAC-OP-02",
            "confidence": mod.get("weight", 0.9),
            "notes": ("LacI repressor binds lac operator; lactose "
                      "(allolactose) or IPTG inactivates LacI and "
                      "relieves operator binding."),
            "source_module": source_type,
        })

    if cap_site and promoter:
        ixs.append({
            "interaction_id": f"lac_reg_{idx}_cap",
            "interaction_type": "stimulation",
            "sbo_term": SBO_STIMULATION,
            "participants": [
                _participant(cap_site, role="stimulator",
                             sbo_role=SBO_STIMULATOR, so_role=SO_OPERATOR,
                             parent_module=source_type),
                _participant(promoter, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_PROMOTER,
                             parent_module=source_type),
            ],
            "rule_id": "INT-CAP-01",
            "confidence": mod.get("weight", 0.9),
            "notes": "CAP-cAMP binds upstream activator site → catabolite "
                     "activation of lac promoter (low glucose state).",
            "source_module": source_type,
        })

    return ixs


# --------------------------------------------------------------------------- #
# Baculovirus recombination cassette — INT-BAC-RECOMB-*
# --------------------------------------------------------------------------- #


def _baculovirus_recombination_interactions(mod, idx):
    """INT-BAC-RECOMB-01: ORF1629 + lef2/ORF603 homology arms drive double
    crossover into BaculoGold / bacmid in insect cells.

    The existing detector doesn't populate `flanking_features`, so we
    reconstruct the homology-arm pair from `mod["features"]` when needed by
    name-matching ORF1629 / lef2 / ORF603 substrings.
    """
    flanks = mod.get("flanking_features") or []
    if len(flanks) < 2:
        # Reconstruct from the module's feature list: name-match for the
        # known baculovirus homology-arm substrings.
        candidates = []
        for f in mod.get("features") or []:
            if not isinstance(f, dict):
                # Some detectors store feature indices, not feature dicts —
                # in that case we can't reconstruct here; bail.
                return []
            n = (f.get("name") or "").lower()
            if "orf1629" in n or "lef2" in n or "orf603" in n:
                candidates.append(f)
        if len(candidates) >= 2:
            candidates.sort(key=lambda x: x.get("start", 0))
            flanks = [candidates[0], candidates[-1]]
        else:
            return []
    flank_5p, flank_3p = flanks[0], flanks[1]
    expression_block = mod.get("expression_block")
    source_type = mod.get("module_type", "baculovirus_recombination_cassette")
    ixs = []

    ixs.append({
        "interaction_id": f"bac_recomb_{idx}_xover",
        "interaction_type": "conversion",
        "sbo_term": SBO_CONVERSION,
        "participants": [
            _participant(flank_5p, role="reactant",
                         sbo_role=SBO_REACTANT, so_role=SO_REGULATORY_REGION,
                         parent_module=source_type),
            _participant(flank_3p, role="reactant",
                         sbo_role=SBO_REACTANT, so_role=SO_REGULATORY_REGION,
                         parent_module=source_type),
        ],
        "rule_id": "INT-BAC-RECOMB-01",
        "confidence": mod.get("weight", 0.9),
        "notes": "ORF1629 + lef2/ORF603 homology arms drive double crossover "
                 "into BaculoGold / bacmid in insect cells.",
        "source_module": source_type,
    })

    if expression_block:
        ixs.append({
            "interaction_id": f"bac_recomb_{idx}_cargo",
            "interaction_type": "stimulation",
            "sbo_term": SBO_STIMULATION,
            "participants": [
                _participant(expression_block, role="template",
                             sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                             parent_module=source_type),
            ],
            "rule_id": "INT-BAC-RECOMB-02",
            "confidence": mod.get("weight", 0.9),
            "notes": "polh/p10-driven payload enclosed by the homology arms.",
            "source_module": source_type,
        })

    return ixs


# --------------------------------------------------------------------------- #
# Tet-inducible expression cassette — INT-TET-IND-*
# --------------------------------------------------------------------------- #


def _tet_inducible_interactions(mod, idx):
    """INT-TET-IND-02 / INT-TET-IND-03: Tet-responsive promoter drives
    payload CDS, with tetO as a modulator. INT-TET-IND-01 (transactivator
    binds operator) is emitted only when a separately-detected
    tet_regulator_cassette is present on the same plasmid (handled at
    cross-module pass — out of scope for this single-module builder).
    """
    subs_by_type = {}
    for s in mod.get("submodules") or []:
        subs_by_type.setdefault(s.get("module_type", ""), []).append(s)
    promoter = (subs_by_type.get("tet_responsive_promoter") or [None])[0]
    operator = (subs_by_type.get("tet_operator") or [None])[0]
    payload = (subs_by_type.get("payload_cds") or [None])[0]
    source_type = mod.get("module_type", "tet_inducible_expression_cassette")
    ixs = []

    if promoter and payload:
        ixs.append({
            "interaction_id": f"tet_ind_{idx}_txn",
            "interaction_type": "genetic_production",
            "sbo_term": SBO_GENETIC_PRODUCTION,
            "participants": [
                _participant(promoter, role="stimulator",
                             sbo_role=SBO_STIMULATOR, so_role=SO_PROMOTER,
                             parent_module=source_type),
                _participant(payload, role="template",
                             sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                             parent_module=source_type),
            ],
            "rule_id": "INT-TET-IND-02",
            "confidence": mod.get("weight", 0.9),
            "notes": "Tet-responsive promoter drives downstream CDS — "
                     "rate gated by transactivator + Dox.",
            "source_module": source_type,
        })

    if operator and promoter:
        ixs.append({
            "interaction_id": f"tet_ind_{idx}_gate",
            "interaction_type": "stimulation",
            "sbo_term": SBO_STIMULATION,
            "participants": [
                _participant(operator, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_OPERATOR,
                             parent_module=source_type),
                _participant(promoter, role="modifier",
                             sbo_role=SBO_MODIFIER, so_role=SO_PROMOTER,
                             parent_module=source_type),
            ],
            "rule_id": "INT-TET-IND-03",
            "confidence": mod.get("weight", 0.9),
            "notes": "tetO controls transcription initiation at the responsive promoter.",
            "source_module": source_type,
        })

    return ixs



# --------------------------------------------------------------------------- #
# Orientation-aware Gateway substrates — INT-REC-GW-INTER / EXC / INV
# --------------------------------------------------------------------------- #


def _gateway_orientation_interactions(mod, idx):
    """Emit one interaction per orientation-aware Gateway substrate module.

    The two att sites are reactants; the intervening cargo (when present)
    is the modifier (deletion target for excision modules, inversion
    target for inversion modules, payload for intermolecular modules).
    """
    subs_by_type = {}
    for s in mod.get("submodules") or []:
        subs_by_type.setdefault(s.get("module_type", ""), []).append(s)
    atts = subs_by_type.get("recombination_site") or []
    cargo = (subs_by_type.get("recombination_cargo") or [None])[0]
    if len(atts) < 2:
        return []
    left, right = atts[0], atts[1]
    source_type = mod.get("module_type", "gateway_recombination")

    # Pick rule_id + interaction type per module type.
    if source_type == "gateway_excision_module":
        rid = "INT-REC-GW-EXC-01"
        ix_type = "conversion"
        sbo_term = SBO_CONVERSION
        notes = ("BP / LR clonase recombines the same-strand att pair "
                 "and EXCISES the intervening cargo as a circular byproduct "
                 "(intramolecular deletion). No second plasmid required.")
    elif source_type == "gateway_inversion_module":
        rid = "INT-REC-GW-INV-01"
        ix_type = "conversion"
        sbo_term = SBO_CONVERSION
        notes = ("BP / LR clonase recombines the opposite-strand outward-"
                 "pointing att pair and INVERTS the intervening cargo "
                 "(intramolecular inversion). No second plasmid required.")
    elif source_type == "gateway_intermolecular_module":
        rid = "INT-REC-GW-INTER-01"
        ix_type = "non_covalent_binding"
        sbo_term = SBO_NON_COVALENT_BINDING
        notes = ("Compatible att pair on opposite strands pointing inward — "
                 "this plasmid is a substrate for intermolecular BP / LR "
                 "recombination with a compatible vector carrying the "
                 "matching att pair.")
    else:
        return []

    participants = [
        _participant(left, role="reactant",
                     sbo_role=SBO_REACTANT, so_role=SO_RECOMBINATION_SITE,
                     parent_module=source_type),
        _participant(right, role="reactant",
                     sbo_role=SBO_REACTANT, so_role=SO_RECOMBINATION_SITE,
                     parent_module=source_type),
    ]
    if cargo and source_type in ("gateway_excision_module",
                                  "gateway_inversion_module"):
        participants.append(
            _participant(cargo, role="modifier",
                         sbo_role=SBO_MODIFIER, so_role=SO_ENGINEERED_REGION,
                         parent_module=source_type)
        )

    return [{
        "interaction_id": f"gw_orient_{idx}",
        "interaction_type": ix_type,
        "sbo_term": sbo_term,
        "participants": participants,
        "rule_id": rid,
        "confidence": mod.get("weight", 0.93),
        "notes": notes,
        "source_module": source_type,
    }]



# --------------------------------------------------------------------------- #
# Phage RNAP (T7 / T3 / SP6) expression cassette — INT-T7-EXPR-01 etc.
# --------------------------------------------------------------------------- #


def _phage_rnap_expression_interactions(mod, idx):
    """Phage-RNAP promoter (T7/T3/SP6) drives a same-strand CDS within
    200 bp downstream (in the promoter's orientation). Emits one
    transcription interaction with the promoter as stimulator and the
    CDS as template."""
    subs_by_type = {}
    for s in mod.get("submodules") or []:
        subs_by_type.setdefault(s.get("module_type", ""), []).append(s)
    promoter = (subs_by_type.get("phage_promoter") or [None])[0]
    payload = (subs_by_type.get("payload_cds") or [None])[0]
    if not (promoter and payload):
        return []
    rid = mod.get("rule_id", "INT-T7-EXPR-01")
    sub = (mod.get("metadata") or {}).get("subtype") or "T7"
    source_type = mod.get("module_type", "phage_rnap_expression_cassette")
    return [{
        "interaction_id": f"phage_rnap_{idx}_txn",
        "interaction_type": "transcription",
        "sbo_term": SBO_TRANSCRIPTION,
        "participants": [
            _participant(promoter, role="stimulator",
                         sbo_role=SBO_STIMULATOR, so_role=SO_PROMOTER,
                         parent_module=source_type),
            _participant(payload, role="template",
                         sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                         parent_module=source_type),
        ],
        "rule_id": rid,
        "confidence": mod.get("weight", 0.92),
        "notes": (
            f"{sub.upper()} phage RNA polymerase promoter drives downstream "
            f"CDS in vitro / in T7-strain hosts."
        ),
        "source_module": source_type,
    }]


_MODULE_BUILDERS = {
    "floxed_region":          _recombination_interaction,
    "lsl_cassette":           _recombination_interaction,
    "floxed_cassette":        _recombination_interaction,
    "frt_flanked_cassette":   _recombination_interaction,
    "gateway_entry_cassette": _recombination_interaction,
    "gateway_dest_cassette":  _recombination_interaction,
    "gateway_recombination":  _recombination_interaction,
    "integrase_landing_pad":  _recombination_interaction,
    "lac_alpha_blue_white_module": _lac_bw_interactions,
    "insulated_expression_block":  _insulator_block_interaction,
    "guide_expression_cassette":   _guide_expression_cassette_interactions,
    # 2026-05-03 audit-driven additions:
    "mobilizable_replicon":              _mobilizable_replicon_interactions,
    "tn3_transposon":                    _tn3_transposon_interactions,
    "lac_promoter_regulatory_unit":      _lac_promoter_regulatory_interactions,
    "baculovirus_recombination_cassette": _baculovirus_recombination_interactions,
    # 2026-05-03 follow-up:
    "tet_inducible_expression_cassette": _tet_inducible_interactions,
    # 2026-05-05: orientation-aware Gateway substrates
    "gateway_excision_module":         _gateway_orientation_interactions,
    "gateway_inversion_module":        _gateway_orientation_interactions,
    "gateway_intermolecular_module":   _gateway_orientation_interactions,
    # 2026-05-06: T7 / T3 / SP6 phage-RNAP expression
    "phage_rnap_expression_cassette":  _phage_rnap_expression_interactions,
}




# --------------------------------------------------------------------------- #
# Bacterial selection interactions (INT-BSEL-01)
# --------------------------------------------------------------------------- #

_BACTERIAL_SELECTION_CDS = {
    "ampr", "bla", "kanr", "neo", "neor", "cmr", "camr", "chloramphenicol",
    "tcr", "tetr", "tetracycline", "zeor", "bleor", "bla(m)",
    "hygr", "purr", "puror", "puromycin", "smr", "specr", "spectinomycin",
    "gmr", "gentamicin", "blasticidin", "bsd",
}

# Bacterial promoter-name fingerprints (substring match, lowercase). These
# either repeat the selection marker name + " promoter" or match recognized
# standalone bacterial promoters (lac, T5, σ70).
_BACTERIAL_PROMOTER_FINGERPRINTS = {
    "ampr promoter", "bla promoter", "kanr promoter", "neor promoter",
    "cmr promoter", "cat promoter", "camr promoter", "chloramphenicol promoter",
    "tcr promoter", "tetr promoter", "tet promoter",
    "puror promoter", "hygr promoter", "zeor promoter", "bleor promoter",
    "em7 promoter",  # hybrid bacterial + weak mammalian
    "pem7", "lac promoter", "plac", "tac promoter", "trc promoter",
    "t5 promoter", "t7 promoter (bacterial)",
    "ampr-promoter",
}


def _is_bacterial_selection_cds(feature: Dict[str, Any]) -> bool:
    name = (feature.get("name") or "").lower()
    ftype = (feature.get("type") or "").lower()
    if ftype not in ("cds", "gene", "protein_generator", "mat_peptide"):
        return False
    # Strip trailing disambiguation suffix (e.g. AmpR-010 → ampr)
    base = name.split("-")[0].split("(")[0].strip()
    return base in _BACTERIAL_SELECTION_CDS or any(
        sel == base or base.startswith(sel) for sel in _BACTERIAL_SELECTION_CDS
    )


def _is_bacterial_promoter(feature: Dict[str, Any]) -> bool:
    name = (feature.get("name") or "").lower()
    ftype = (feature.get("type") or "").lower()
    kb_class = (feature.get("kb_class") or "").lower()
    if ftype not in ("promoter", "regulatory") and kb_class != "promoter":
        return False
    # Bacterial selection-adjacent promoter names
    return any(fp in name for fp in _BACTERIAL_PROMOTER_FINGERPRINTS)


def _bacterial_selection_interactions(
    feature_rows: Optional[List[Dict[str, Any]]],
) -> List[Dict[str, Any]]:
    """Emit INT-BSEL-01 for each bacterial promoter → bacterial-selection-CDS
    pair found on the same strand, within 500 bp, without emitting a wrapping
    module. Optionally include a terminator within 20 bp of the CDS end as a
    modifier participant."""
    if not feature_rows:
        return []

    promoters = [f for f in feature_rows if _is_bacterial_promoter(f)]
    cdses = [f for f in feature_rows if _is_bacterial_selection_cds(f)]
    terminators = [
        f for f in feature_rows
        if (f.get("type") or "").lower() == "terminator"
        or "terminator" in (f.get("name") or "").lower()
    ]

    def _prom_names_cds(prom_name: str, cds_name: str) -> bool:
        """True if the promoter name explicitly contains the selection CDS base name."""
        p = (prom_name or "").lower()
        c = (cds_name or "").lower().split("-")[0].split("(")[0].strip()
        return c and c in p

    ixs: List[Dict[str, Any]] = []
    used_pairs: set = set()
    for prom in promoters:
        p_strand = prom.get("strand") or prom.get("direction") or 1
        prom_name = prom.get("name") or ""
        # Find the best selection CDS within 500 bp. Prefer same-strand pairs;
        # fall back to opposite-strand pairs only when the promoter name
        # explicitly names the CDS (e.g. "AmpR promoter" ↔ "AmpR").
        best_cds = None
        best_gap = None
        best_same_strand = None
        for cds in cdses:
            c_strand = cds.get("strand") or cds.get("direction") or 1
            same_strand = (c_strand == p_strand)
            name_bond = _prom_names_cds(prom_name, cds.get("name") or "")
            if not same_strand and not name_bond:
                continue
            # Distance by the promoter's putative output direction
            if p_strand >= 0:
                p_end = prom.get("end", 0)
                gap = cds.get("start", 0) - p_end
            else:
                p_end = prom.get("start", 0)
                gap = p_end - cds.get("end", 0)
            # Name-bonded pairs: allow coordinate-based proximity regardless
            # of direction convention (use absolute distance).
            if name_bond and not same_strand:
                prom_midpoint = (prom.get("start", 0) + prom.get("end", 0)) // 2
                cds_midpoint = (cds.get("start", 0) + cds.get("end", 0)) // 2
                gap = abs(prom_midpoint - cds_midpoint)
                max_gap = 1500  # name-bonded can be larger (bla sig/mat peptide spreads)
            else:
                max_gap = 500
            # Allow small negative gaps for promoters whose annotation edge
            # touches/overlaps the CDS by a few bp (pLannotate often draws the
            # boundary 1 bp into the CDS, esp. on -1 strand AmpR cassettes).
            if -10 <= gap <= max_gap:
                # Prefer same-strand pairs over name-bonded opposite-strand pairs
                better = False
                if best_cds is None:
                    better = True
                elif best_same_strand is False and same_strand:
                    better = True  # upgrade from opposite-strand name-bond to same-strand
                elif best_same_strand == same_strand and gap < best_gap:
                    better = True
                if better:
                    best_cds = cds
                    best_gap = gap
                    best_same_strand = same_strand
        if best_cds is None:
            continue
        key = (prom.get("start"), prom.get("end"), best_cds.get("start"), best_cds.get("end"))
        if key in used_pairs:
            continue
        used_pairs.add(key)

        # Look for an optional terminator within 20 bp of CDS end.
        # Use the CDS's own strand (not promoter's) so name-bonded opposite-strand
        # pairs still find a terminator downstream of the CDS.
        c_strand = best_cds.get("strand") or best_cds.get("direction") or 1
        cds_end = best_cds.get("end", 0) if c_strand >= 0 else best_cds.get("start", 0)
        term_participant = None
        for term in terminators:
            t_strand = term.get("strand") or term.get("direction") or 1
            if t_strand != c_strand:
                continue
            if c_strand >= 0:
                t_gap = term.get("start", 0) - cds_end
            else:
                t_gap = cds_end - term.get("end", 0)
            if 0 <= t_gap <= 20:
                term_participant = term
                break

        # Peer interaction 1: promoter drives the selection CDS (genetic_production)
        ixs.append({
            "interaction_id": f"bsel_prom_{prom.get('start','?')}_{best_cds.get('start','?')}",
            "interaction_type": "genetic_production",
            "sbo_term": SBO_GENETIC_PRODUCTION,
            "participants": [
                _participant(prom, role="stimulator",
                             sbo_role=SBO_STIMULATOR, so_role=SO_PROMOTER,
                             parent_module="bacterial_selection_cassette"),
                _participant(best_cds, role="template",
                             sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                             parent_module="bacterial_selection_cassette"),
            ],
            "rule_id": "INT-BSEL-01",
            "confidence": 0.9,
            "notes": (
                f"{prom.get('name','bacterial promoter')} drives expression of "
                f"{best_cds.get('name','bacterial selection marker')}."
            ),
            "source_module": "bacterial_selection_cassette",
        })

        # Peer interaction 2: terminator regulates the selection CDS transcript
        # Only emitted when a terminator sits within 20 bp of the CDS end.
        if term_participant is not None:
            ixs.append({
                "interaction_id": f"bsel_term_{prom.get('start','?')}_{best_cds.get('start','?')}",
                "interaction_type": "stimulation",
                "sbo_term": SBO_STIMULATION,
                "participants": [
                    _participant(term_participant, role="modifier",
                                 sbo_role=SBO_MODIFIER, so_role=None,
                                 parent_module="bacterial_selection_cassette"),
                    _participant(best_cds, role="template",
                                 sbo_role=SBO_TEMPLATE, so_role=SO_CDS,
                                 parent_module="bacterial_selection_cassette"),
                ],
                "rule_id": "INT-BSEL-TERM-01",
                "confidence": 0.85,
                "notes": (
                    f"{term_participant.get('name','terminator')} terminates "
                    f"transcription of {best_cds.get('name','bacterial selection marker')}."
                ),
                "source_module": "bacterial_selection_cassette",
            })

    return ixs

# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #


def build_interactions(
    features: Optional[List[Dict[str, Any]]] = None,
    rule_based_modules: Optional[List[Dict[str, Any]]] = None,
    cds_submodules: Optional[List[Dict[str, Any]]] = None,
    *,
    mammalian_pol2_cassettes: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """Emit module-scoped interactions.

    Args:
        features: pLannotate feature rows (name, type, start, end, strand).
            Used only as a fallback name source for the CDS of a Pol II
            cassette when the submodule parser produced nothing. No
            interactions are inferred from these directly.
        rule_based_modules: output of RuleBasedModuleDetector.detect_modules.
        cds_submodules: output of resolve_cds_submodules().
        mammalian_pol2_cassettes: list of MammalianPol2Cassette.to_dict()
            outputs; authoritative source for Pol II expression interactions.

    Returns:
        A flat list of interaction dicts. Every interaction carries a
        `source_module` field naming the module that produced it.
    """
    interactions: List[Dict[str, Any]] = []

    # 1) Pol II / lentiviral expression cassettes — authoritative pro→CDS→polyA
    for i, cas in enumerate(mammalian_pol2_cassettes or []):
        interactions.extend(
            _pol2_cassette_interactions(
                cas, i,
                cds_submodules=cds_submodules,
                feature_rows=features,
            )
        )

    # 2) Rule-based modules — dispatch by module_type
    for i, mod in enumerate(rule_based_modules or []):
        builder = _MODULE_BUILDERS.get(mod.get("module_type", ""))
        if builder:
            interactions.extend(builder(mod, i))

    # 2b) Cross-module lentiviral interactions (upstream → payload → downstream)
    interactions.extend(_lentiviral_three_module_interactions(rule_based_modules or []))

    # 2c) Bacterial selection interactions (promoter → selection CDS, no wrapper)
    interactions.extend(_bacterial_selection_interactions(features))

    # 3) CDS fusion interactions (P2A/T2A ribosomal skip)
    interactions.extend(_cds_fusion_interactions(
        cds_submodules or [],
        mammalian_pol2_cassettes=mammalian_pol2_cassettes,
        feature_rows=features,
    ))

    return interactions


__all__ = ["build_interactions"]
