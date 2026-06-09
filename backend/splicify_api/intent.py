"""
Intent parsing for AI Plasmid Design chat endpoint.

Deterministic, LLM-free. Classifies the prompt into one of nine intents,
extracts handler-specific parameters, and pre-resolves any KB-known parts the
user named so downstream handlers do not have to re-parse the prompt.

Intents (matches chat.py dispatch):
  annotate_gb, gateway_cloning, gibson_design, plasmid_design, sdm_design,
  sgrna_golden_gate, golden_gate_primer_design, restriction_cloning, unknown

Returns the same schema the previous LLM-backed version emitted, plus a
`kb_resolved` block carrying part candidates and KB-identified features.
"""
from __future__ import annotations

import logging
import re
from typing import Any, Dict, List, Optional

logger = logging.getLogger("intent")


# ---------------------------------------------------------------------------
# Empty sub-objects
# ---------------------------------------------------------------------------
def _empty_gibson() -> Dict[str, Any]:
    return {"assembly": None, "primer_params": None}


def _empty_pcr() -> Dict[str, Any]:
    return {
        "target_start": None,
        "target_length": None,
        "product_size_min": None,
        "product_size_max": None,
        "primer_min_tm": None,
        "primer_opt_tm": None,
        "primer_max_tm": None,
        "excluded_mode": None,
        "excluded_start": None,
        "excluded_end": None,
        "excluded_flank": None,
    }


def _empty_sdm() -> Dict[str, Any]:
    return {
        "mutation_type": None,
        "target_method": None,
        "target_feature_name": None,
        "target_position_start": None,
        "target_position_end": None,
        "codon_position": None,
        "codon_from": None,
        "codon_to": None,
        "terminus": None,
        "old_sequence": None,
        "new_sequence": None,
        "description": None,
    }


def _empty_sgrna() -> Dict[str, Any]:
    return {"grna_sequence": None, "vector_name": None, "enzyme": None}


def _empty_golden_gate() -> Dict[str, Any]:
    return {"workflow_type": None, "enzyme": None}


def _empty_restriction_cloning() -> Dict[str, Any]:
    return {
        "insert_name": None,
        "insert_sequence": None,
        "vector_name": None,
        "enzymes": None,
    }


def _empty_kb_resolved() -> Dict[str, Any]:
    return {"candidates": [], "identified": []}


# ---------------------------------------------------------------------------
# Keyword tables
# ---------------------------------------------------------------------------
_GATEWAY_KEYWORDS = (
    "gateway", "bp reaction", "lr reaction",
    "attb", "attp", "attl", "attr", "pdonr", "pdest",
)

_GRNA_KEYWORDS = re.compile(r"\b(sgrna|grna|guide\s*rna)\b", re.IGNORECASE)

_GG_OLIGO_CONTEXT = (
    "golden gate", "goldengate", "oligo", "bsmbi", "bbsi", "bsai",
    "lenticrispr", "px330", "px335",
)

_TYPE2_ENZYMES = {
    "ecori": "EcoRI", "hindiii": "HindIII", "bamhi": "BamHI", "xhoi": "XhoI",
    "nhei": "NheI", "kpni": "KpnI", "saci": "SacI", "psti": "PstI",
    "sali": "SalI", "agei": "AgeI", "mlui": "MluI", "bglii": "BglII",
    "xbai": "XbaI", "spei": "SpeI", "noti": "NotI", "asci": "AscI",
    "clai": "ClaI",
}

_RESTRICTION_PHRASES = (
    "restriction cloning", "restriction digest cloning",
    "restriction digest clone", "restriction-based cloning",
    "restriction enzyme cloning", "restriction enzyme workflow",
    "restriction enzyme assembly", "restriction enzyme digest",
    "restriction enzyme", "restriction digest",
    "traditional cloning", "re cloning",
)

_SDM_KEYWORDS = (
    "mutate", "mutation", "substitute", "substitution",
    "delete the", "deletion", "remove the", "remove bp", "drop the",
    "insert ", "insertion", "change codon", "point mutation",
    "site-directed", "site directed", "sdm", "mutagenesis",
    "amino acid",
)

_AA_NOTATION = re.compile(r"\b([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])\b")

_GG_PRIMER_KEYWORDS = (
    "golden gate primer", "golden gate primers", "type iis primer",
    "bsai primer", "bsmbi primer", "modular assembly", "moclo",
    "orthogonal overhang", "multi-fragment golden gate",
    "multi fragment golden gate", "golden gate assembly primer",
    "type iis assembly", "golden gate assembly",
    "golden gate assemble", "golden gate to assemble",
)

_PLASMID_DESIGN_KEYWORDS = (
    "make a", "make me a", "build a", "build me a", "create a", "design a",
    "i need a", "i want a", "i'd like a", "id like a",
    "expression vector", "lentiviral vector", "mammalian vector",
    "bacterial vector", "shuttle vector",
    "plasmid that expresses", "plasmid expressing",
    "describe a plasmid", "describe the plasmid",
    "vector that expresses", "vector expressing",
    "swap", "replace", "add a cassette", "add puromycin", "add selection",
    "delete the", "delete a", "remove the", "remove a", "drop the", "drop a",
    "change the", "move the",
)

_GIBSON_KEYWORDS = ("gibson", "isothermal assembly")
_GIBSON_ASSEMBLY_VERBS = ("assembly", "assemble")

_ANNOTATE_KEYWORDS = ("annotate", "features", "feature", "identify")


# ---------------------------------------------------------------------------
# Helper builders
# ---------------------------------------------------------------------------
def _result(
    intent: str,
    confidence: float,
    *,
    gibson_assembly: Optional[str] = None,
    sdm_params: Optional[Dict[str, Any]] = None,
    sgrna_params: Optional[Dict[str, Any]] = None,
    golden_gate_params: Optional[Dict[str, Any]] = None,
    restriction_params: Optional[Dict[str, Any]] = None,
    errors: Optional[List[str]] = None,
) -> Dict[str, Any]:
    gib = _empty_gibson()
    if gibson_assembly:
        gib["assembly"] = gibson_assembly
    return {
        "intent": intent,
        "gibson_design": gib,
        "pcr": _empty_pcr(),
        "sdm": sdm_params or _empty_sdm(),
        "sgrna": sgrna_params or _empty_sgrna(),
        "golden_gate": golden_gate_params or _empty_golden_gate(),
        "restriction_cloning": restriction_params or _empty_restriction_cloning(),
        "kb_resolved": _empty_kb_resolved(),
        "router_notes": {"confidence": confidence, "errors": errors or []},
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
async def parse_intent(
    message: str,
    has_target: bool = False,
    has_inventory: bool = False,
    seq_count: int = 0,
    redacted_message: str = "",
) -> Dict[str, Any]:
    """Classify intent + extract params + resolve KB-known parts.

    The signature is preserved (including `redacted_message` and `async`) so
    callers in chat.py do not need to change. This implementation is fully
    synchronous and deterministic — no LLM calls.
    """
    text = redacted_message or message
    msg = text.lower()

    result = _classify(text, msg, has_target, has_inventory, seq_count)

    # Pre-resolve any KB-known parts the user named. Done here once so every
    # downstream handler reads pre-resolved hits instead of re-parsing the
    # prompt. Imports are deferred to avoid a chat ↔ intent import cycle.
    try:
        from .chat import extract_part_candidates, identify_features_from_kb
        candidates = extract_part_candidates(text)
        identified = identify_features_from_kb(candidates) if candidates else []
        result["kb_resolved"] = {
            "candidates": candidates,
            "identified": identified,
        }
    except Exception as exc:  # pragma: no cover — best effort
        logger.warning("KB pre-resolve skipped: %s", exc)

    # v13: surface a design-completeness verdict on every classified
    # intent so the dispatcher can warn early when the user request
    # lacks elements needed for a complete expression interaction
    # network (CDS without promoter, orphan tag, single recombination
    # site, etc.). Same role contract as the post-assembly verifier.
    try:
        from .target_from_inventory_router import analyze_design_intent
        result["design_completeness"] = analyze_design_intent(result, text)
    except Exception as exc:  # pragma: no cover
        logger.warning("design-intent analysis skipped: %s", exc)

    return result


# ---------------------------------------------------------------------------
# Classifier (priority-ordered, deterministic)
# ---------------------------------------------------------------------------
def _classify(
    text: str,
    msg: str,
    has_target: bool,
    has_inventory: bool,
    seq_count: int,
) -> Dict[str, Any]:
    has_plasmid_file = has_target or has_inventory

    # 1. Gateway: explicit Gateway keywords + a file present
    if has_plasmid_file and any(kw in msg for kw in _GATEWAY_KEYWORDS):
        return _result("gateway_cloning", 0.95)

    # 2. sgRNA Golden Gate (annealed-oligo into CRISPR vector)
    has_grna = bool(_GRNA_KEYWORDS.search(msg))
    has_gg_or_oligo = any(kw in msg for kw in _GG_OLIGO_CONTEXT)
    if has_grna and has_gg_or_oligo:
        return _result(
            "sgrna_golden_gate",
            0.90,
            sgrna_params=_extract_sgrna_params(text, msg),
        )

    # 3. Golden Gate multi-fragment primer design
    if any(kw in msg for kw in _GG_PRIMER_KEYWORDS):
        if _has_multi_fragment_context(msg, seq_count):
            return _result(
                "golden_gate_primer_design",
                0.90,
                golden_gate_params=_extract_golden_gate_params(msg),
            )

    # 4. Restriction cloning: explicit phrase, or "clone X into Y" + Type II enzyme
    has_restriction_phrase = any(p in msg for p in _RESTRICTION_PHRASES)
    has_type2_enzyme = any(enz in msg for enz in _TYPE2_ENZYMES)
    has_clone_verb = (
        "clone" in msg or "cloning" in msg
        or "insert " in msg or "inserting" in msg
    )
    if has_restriction_phrase or (has_clone_verb and has_type2_enzyme):
        return _result(
            "restriction_cloning",
            0.85,
            restriction_params=_extract_restriction_params(msg),
        )

    # 4b. "clone X into Y" with a target file but no other workflow signal
    # → restriction_cloning is the right default (the previous code fell
    # through to the annotate_gb has_target fallback).
    has_other_workflow_kw = (
        any(kw in msg for kw in _GIBSON_KEYWORDS)
        or any(kw in msg for kw in _GG_PRIMER_KEYWORDS)
        or any(kw in msg for kw in _GATEWAY_KEYWORDS)
        or bool(_GRNA_KEYWORDS.search(msg))
    )
    if has_clone_verb and has_plasmid_file and not has_other_workflow_kw:
        return _result(
            "restriction_cloning",
            0.70,
            restriction_params=_extract_restriction_params(msg),
        )

    # 5. SDM: a plasmid file is present + mutation language (or AA notation like Y66H)
    has_aa_notation = bool(_AA_NOTATION.search(text))
    has_sdm_keyword = any(kw in msg for kw in _SDM_KEYWORDS)
    if has_plasmid_file and (has_aa_notation or has_sdm_keyword):
        return _result(
            "sdm_design",
            0.85,
            sdm_params=_extract_sdm_params(text, msg),
        )

    # 6. plasmid_design: describe-a-plasmid language. Checked before gibson
    #    since "design" is ambiguous.
    if any(kw in msg for kw in _PLASMID_DESIGN_KEYWORDS):
        if "gibson assembly" not in msg and "pcr primer" not in msg:
            return _result("plasmid_design", 0.80)

    # 7. annotate_gb: file uploaded + annotate/feature language
    if has_target and any(kw in msg for kw in _ANNOTATE_KEYWORDS):
        return _result("annotate_gb", 0.85)

    # 8. Gibson: explicit gibson keyword, or "assemble/assembly" with multiple
    #    sequences
    has_gibson_kw = any(kw in msg for kw in _GIBSON_KEYWORDS)
    has_assemble_verb = any(kw in msg for kw in _GIBSON_ASSEMBLY_VERBS)
    if has_gibson_kw or (has_assemble_verb and seq_count >= 2):
        return _result("gibson_design", 0.90, gibson_assembly=_infer_gibson_assembly(msg))

    # 9. Multiple sequences with no other signal → likely Gibson
    if seq_count >= 2:
        return _result("gibson_design", 0.65, gibson_assembly=_infer_gibson_assembly(msg))

    # 10. Plasmid file with no clear ask → annotate
    if has_target:
        return _result("annotate_gb", 0.55)

    return _result("unknown", 0.40)


# ---------------------------------------------------------------------------
# Param extractors
# ---------------------------------------------------------------------------
def _infer_gibson_assembly(msg: str) -> Optional[str]:
    if "linear" in msg:
        return "linear"
    if any(k in msg for k in ("circular", "plasmid", "vector")):
        return "circular"
    return None


def _has_multi_fragment_context(msg: str, seq_count: int) -> bool:
    return (
        seq_count >= 2
        or any(kw in msg for kw in (
            "fragments", "multi-fragment", "3 fragment", "three fragment",
            "assemble these", "assemble",
        ))
        or msg.count("+") >= 1
        or msg.count(",") >= 2
    )


def _extract_sgrna_params(text: str, msg: str) -> Dict[str, Any]:
    out = _empty_sgrna()
    m = re.search(r"\b([ACGT]{17,30})\b", text.upper())
    if m:
        out["grna_sequence"] = m.group(1)
    if "lenticrispr" in msg:
        out["vector_name"] = "lentiCRISPR v2"
    elif "px330" in msg:
        out["vector_name"] = "pX330"
    elif "px335" in msg:
        out["vector_name"] = "pX335"
    if "bsmbi" in msg:
        out["enzyme"] = "BsmBI"
    elif "bbsi" in msg:
        out["enzyme"] = "BbsI"
    elif "bsai" in msg:
        out["enzyme"] = "BsaI"
    return out


def _extract_golden_gate_params(msg: str) -> Dict[str, Any]:
    out = _empty_golden_gate()
    if "deletion" in msg or "delete" in msg:
        out["workflow_type"] = "scarless_deletion"
    elif "replacement" in msg or "replace" in msg:
        out["workflow_type"] = "single_fragment"
    else:
        out["workflow_type"] = "multi_fragment"
    if "bsai" in msg:
        out["enzyme"] = "BsaI"
    elif "bsmbi" in msg:
        out["enzyme"] = "BsmBI"
    elif "bbsi" in msg:
        out["enzyme"] = "BbsI"
    return out


def _extract_restriction_params(msg: str) -> Dict[str, Any]:
    out = _empty_restriction_cloning()
    found: List[str] = []
    for enz_lower, enz_canonical in _TYPE2_ENZYMES.items():
        if enz_lower in msg and enz_canonical not in found:
            found.append(enz_canonical)
    if len(found) >= 2:
        out["enzymes"] = found[:2]

    m = re.search(
        r"clon(?:e|ing)\s+([a-zA-Z0-9_\- ]+?)\s+(?:in|into)\s+([a-zA-Z0-9_\-]+)",
        msg,
    )
    if m:
        out["insert_name"] = m.group(1).strip()
        out["vector_name"] = m.group(2).strip()
    return out


def _extract_sdm_params(text: str, msg: str) -> Dict[str, Any]:
    out = _empty_sdm()

    # Amino acid notation: Y66H, S65T, D10A
    aa = _AA_NOTATION.search(text)
    if aa:
        out["target_method"] = "codon"
        out["mutation_type"] = "substitution"
        out["codon_from"] = aa.group(1)
        out["codon_position"] = int(aa.group(2))
        out["codon_to"] = aa.group(3)
        # Try to scope to a feature: "Y66H on Cas9" / "Y66H in eGFP"
        scope = re.search(
            r"%s\s+(?:on|in|of)\s+([A-Za-z0-9_\-]+)" % re.escape(aa.group(0)),
            text,
        )
        if scope:
            out["target_feature_name"] = scope.group(1)
        return out

    # Mutation type from verbs
    if any(kw in msg for kw in ("delete", "deletion", "remove", "drop")):
        out["mutation_type"] = "deletion"
    elif any(kw in msg for kw in ("insert", "insertion", "add ")):
        out["mutation_type"] = "insertion"
    else:
        out["mutation_type"] = "substitution"

    # Position-based: "delete bp 100-150"
    pos = re.search(r"bp\s+(\d+)\s*[-–]\s*(\d+)", msg)
    if pos:
        out["target_method"] = "position"
        out["target_position_start"] = int(pos.group(1))
        out["target_position_end"] = int(pos.group(2))
        return out

    # Feature-based: "delete the His-tag", "remove the NLS",
    # "drop the FLAG tag from cas9" (capture "FLAG", not "FLAG tag from cas9").
    # Allow an optional descriptor suffix (tag/peptide/sequence/site/element)
    # then stop at a connector word, end-of-string, or punctuation.
    feat = re.search(
        r"(?:delete|remove|drop)\s+(?:the\s+)?([A-Za-z0-9][A-Za-z0-9_\-]*)"
        r"(?:\s+(?:tag|tags|peptide|sequence|site|element|cassette|signal|domain))?"
        r"(?:\s+(?:from|in|on|of|at|near|around|inside|out\s+of|to)\b|\s*$|[.,;:!?])",
        msg,
    )
    if feat:
        out["target_method"] = "feature"
        out["target_feature_name"] = feat.group(1).strip()

    # Terminus
    if "n-terminus" in msg or "n terminus" in msg or "n-terminal" in msg:
        out["terminus"] = "N"
    elif "c-terminus" in msg or "c terminus" in msg or "c-terminal" in msg:
        out["terminus"] = "C"

    return out
