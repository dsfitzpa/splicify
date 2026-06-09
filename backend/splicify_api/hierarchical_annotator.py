from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .Module_Library_gb.module_extractor import (
    Feature,
    Module,
    _enrich_module_metadata,
    choose_candidates_bacterial_marker_cassette,
    choose_candidates_bacterial_expression,
    choose_candidates_cds_only,
    choose_candidates_lacI_regulator,
    choose_candidates_lenti,
    choose_candidates_mammalian_lentiviral_expression,
    choose_candidates_mammalian_replication,
    choose_candidates_linker_and_targeting,
    choose_candidates_marker,
    choose_candidates_origin,
    choose_candidates_pol2,
    choose_candidates_pol3,
    choose_candidates_regulatory_elements,
    choose_candidates_shuttle_vector,
    choose_candidates_sv40_neo,
    choose_candidates_yeast_integrating,
    choose_candidates_yeast_pol2,
    contains_feature,
    extract_sequence,
    module_length,
    seq_hash,
    sha256_text,
)
from .plasmid_analyzer import analyze_plasmid_from_plannotate, KnowledgeBase
from .module_parser import MotifDetector


# =============================================================================
# KOZAK SEQUENCE DETECTION
# =============================================================================

# Kozak consensus patterns by organism class
# Format: (pattern_regex, name, organism, strength)
# Patterns are designed to match upstream context + ATG + downstream context
KOZAK_PATTERNS = [
    # Mammalian - Strong Kozak (matches at -3 A/G AND +4 G)
    (re.compile(r'GCC[AG]CCATGG', re.IGNORECASE), "Kozak (strong)", "mammalian", "strong"),
    (re.compile(r'[AG]CCATGG', re.IGNORECASE), "Kozak (adequate)", "mammalian", "adequate"),
    (re.compile(r'[ACGT]{3}ATGG', re.IGNORECASE), "Kozak (+4 G)", "mammalian", "weak"),
    (re.compile(r'[AG][ACGT]{2}ATG[ACGT]', re.IGNORECASE), "Kozak (-3 purine)", "mammalian", "weak"),

    # Insect (Drosophila) - [AC]AA[AC]ATG[AC][AC][CG]
    (re.compile(r'[AC]AA[AC]ATG[AC]', re.IGNORECASE), "Kozak (insect)", "insect", "strong"),

    # Yeast - aAaAaAATGTCt
    (re.compile(r'A[ACGT]A[ACGT]A[ACGT]ATGTC', re.IGNORECASE), "Kozak (yeast)", "yeast", "strong"),
    (re.compile(r'[ACGT]{3}ATGTC', re.IGNORECASE), "Kozak (yeast partial)", "yeast", "adequate"),

    # Plant - acAACAATGGC
    (re.compile(r'[ACGT]{2}AACAATGGC', re.IGNORECASE), "Kozak (plant)", "plant", "strong"),
    (re.compile(r'AACAATGG', re.IGNORECASE), "Kozak (plant partial)", "plant", "adequate"),

    # Plasmodium - taaAAAATGAan
    (re.compile(r'[ACGT]{3}AAAATGA[ACGT]{2}', re.IGNORECASE), "Kozak (plasmodium)", "plasmodium", "strong"),

    # Ciliates - nTaAAAATGRct (R = A or G)
    (re.compile(r'[ACGT]T[ACGT]AAAATG[AG]', re.IGNORECASE), "Kozak (ciliate)", "ciliate", "strong"),
]

# Extended search window - look for ATG and then check upstream/downstream context
KOZAK_SEARCH_WINDOW = 15  # bp upstream of CDS to search for Kozak


def detect_kozak_sequence(
    sequence: str,
    cds_start: int,  # 0-indexed start of CDS
    cds_strand: int,
    seq_len: int,
    circular: bool = True
) -> Optional[Dict[str, Any]]:
    """
    Detect Kozak sequence at the start of a CDS.

    Args:
        sequence: Full plasmid sequence
        cds_start: 0-indexed start position of CDS
        cds_strand: 1 for forward, -1 for reverse
        seq_len: Total sequence length
        circular: Whether plasmid is circular

    Returns:
        Dict with kozak info or None if not found
    """
    if not sequence or cds_start < 0:
        return None

    # Get sequence context around CDS start
    # We need ~10bp upstream and ~4bp downstream of ATG
    if cds_strand == 1:
        # Forward strand - get upstream context
        context_start = max(0, cds_start - 10)
        context_end = min(seq_len, cds_start + 7)  # ATG + 4 more

        if circular and cds_start < 10:
            # Handle wrap-around
            upstream = sequence[seq_len - (10 - cds_start):] + sequence[:cds_start]
            downstream = sequence[cds_start:context_end]
            context_seq = upstream + downstream
            context_start = cds_start - 10  # Will be negative, need to handle
        else:
            context_seq = sequence[context_start:context_end].upper()

        # Find ATG position in context
        atg_pos = context_seq.find('ATG')
        if atg_pos == -1:
            return None

    else:
        # Reverse strand - need reverse complement
        # CDS start on reverse strand means the ATG is at the complement
        context_end = min(seq_len, cds_start + 10)
        context_start = max(0, cds_start - 4)

        context_seq = sequence[context_start:context_end]
        # Reverse complement
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'N': 'N', 'n': 'n'}
        context_seq = ''.join(complement.get(b, 'N') for b in reversed(context_seq)).upper()

        atg_pos = context_seq.find('ATG')
        if atg_pos == -1:
            return None

    # Try each Kozak pattern
    best_match = None
    best_strength_rank = {'strong': 3, 'adequate': 2, 'weak': 1}

    for pattern, name, organism, strength in KOZAK_PATTERNS:
        match = pattern.search(context_seq)
        if match:
            # Check if this match includes the ATG we found
            match_start = match.start()
            match_end = match.end()

            # The ATG should be within or at the end of the match
            if match_start <= atg_pos < match_end:
                current_rank = best_strength_rank.get(strength, 0)
                if best_match is None or current_rank > best_strength_rank.get(best_match['strength'], 0):
                    # Calculate actual genomic positions
                    if cds_strand == 1:
                        kozak_start = context_start + match_start
                        kozak_end = context_start + match_end
                    else:
                        # Reverse strand - flip positions
                        kozak_end = context_end - match_start
                        kozak_start = context_end - match_end

                    best_match = {
                        'name': name,
                        'organism': organism,
                        'strength': strength,
                        'start': kozak_start,
                        'end': kozak_end,
                        'sequence': match.group(),
                        'strand': cds_strand,
                    }

    return best_match


def detect_kozak_for_cds_features(
    sequence: str,
    features: List[Dict[str, Any]],
    circular: bool = True
) -> List[Dict[str, Any]]:
    """
    Detect Kozak sequences for all CDS features.

    Returns list of Kozak annotation dictionaries.
    """
    if not sequence:
        return []

    seq_len = len(sequence)
    kozak_annotations = []

    for feat in features:
        # Check if this is a CDS-like feature
        feat_type = str(feat.get('type', feat.get('Type', ''))).lower()
        feat_role = str(feat.get('role', '')).lower()

        is_cds = (
            feat_type == 'cds' or
            'payload' in feat_role or
            'nuclease' in feat_type or
            'cas9' in str(feat.get('name', feat.get('label', ''))).lower()
        )

        if not is_cds:
            continue

        # Get CDS position
        start = feat.get('start', 0)
        end = feat.get('end', 0)
        direction = feat.get('direction', feat.get('strand', 1))
        strand = -1 if direction == -1 or direction == 'reverse' else 1

        # Detect Kozak at CDS start
        cds_start = start if strand == 1 else end

        kozak = detect_kozak_sequence(sequence, cds_start, strand, seq_len, circular)

        if kozak:
            kozak_annotations.append({
                'name': kozak['name'],
                'start': kozak['start'],
                'end': kozak['end'],
                'direction': kozak['strand'],
                'color': '#10B981',  # Emerald green for Kozak
                'layer': 'kozak',
                'source': 'hierarchical_annotator',
                'type': 'regulatory',
                'metadata': {
                    'organism': kozak['organism'],
                    'strength': kozak['strength'],
                    'sequence': kozak['sequence'],
                    'parent_cds': feat.get('name', feat.get('label', 'unknown')),
                }
            })

    return kozak_annotations


# ──────────────────────────────────────────────
# Grammar rules table
# ──────────────────────────────────────────────

_RULES_CSV = Path(__file__).parent / "annotator/feature_db/schema/hierarchical_feature_rules.csv"


@dataclass
class FeatureRule:
    sseqid: str
    feature_name: str
    feature_type: str
    functional_bucket: str
    taxonomic_scope: str
    nested_module_role: str  # hard_start | fallback_start | hard_end | internal
    starts_module_types: List[str]
    ends_module_types: List[str]
    preferred_parent_modules: List[str]
    boundary_priority: int
    classification_confidence: str


_RULES_BY_SSEQID: Optional[Dict[str, FeatureRule]] = None
_RULES_BY_NAME: Optional[Dict[str, FeatureRule]] = None


def _load_rules() -> Tuple[Dict[str, FeatureRule], Dict[str, FeatureRule]]:
    global _RULES_BY_SSEQID, _RULES_BY_NAME
    if _RULES_BY_SSEQID is not None:
        assert _RULES_BY_NAME is not None
        return _RULES_BY_SSEQID, _RULES_BY_NAME

    _RULES_BY_SSEQID = {}
    _RULES_BY_NAME = {}

    if not _RULES_CSV.exists():
        return _RULES_BY_SSEQID, _RULES_BY_NAME

    def _split(s: str) -> List[str]:
        return [x.strip() for x in (s or "").split(";") if x.strip()]

    with _RULES_CSV.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rule = FeatureRule(
                sseqid=row["sseqid"],
                feature_name=row["Feature"],
                feature_type=row["Type"],
                functional_bucket=row["functional_bucket"],
                taxonomic_scope=row["taxonomic_scope"],
                nested_module_role=row["nested_module_role"],
                starts_module_types=_split(row["starts_module_types"]),
                ends_module_types=_split(row["ends_module_types"]),
                preferred_parent_modules=_split(row["preferred_parent_modules"]),
                boundary_priority=int(row.get("boundary_priority") or 0),
                classification_confidence=row["classification_confidence"],
            )
            _RULES_BY_SSEQID[row["sseqid"]] = rule
            name_key = row["Feature"].lower().strip()
            if name_key not in _RULES_BY_NAME:
                _RULES_BY_NAME[name_key] = rule

    return _RULES_BY_SSEQID, _RULES_BY_NAME


def _lookup_rule(sseqid: str, feature_name: str) -> Optional[FeatureRule]:
    by_id, by_name = _load_rules()
    if sseqid and sseqid in by_id:
        return by_id[sseqid]
    return by_name.get((feature_name or "").lower().strip())


def _get_rule_for_feature(feat: Feature) -> Optional[FeatureRule]:
    sseqid = (feat.qualifiers or {}).get("sseqid", "") or ""
    return _lookup_rule(sseqid, feat.name or "")


# ──────────────────────────────────────────────
# Module type family compatibility
# ──────────────────────────────────────────────

# Multi-scope module types share a family; a generic end can close a specific-scope start.
_MODULE_TYPE_FAMILY: Dict[str, str] = {
    "pol2_expression_animal": "pol2_expression",
    "pol2_expression_plant": "pol2_expression",
    "pol2_expression_fungal": "pol2_expression",
    "pol2_expression_insect": "pol2_expression",
    "pol2_expression_lentiviral": "pol2_expression",
    "pol2_expression_generic": "pol2_expression",
    "pol3_expression_animal": "pol3_expression",
    "pol3_expression_fungal": "pol3_expression",
    "pol3_expression_insect": "pol3_expression",
    "pol3_u6_sgrna_cassette": "pol3_expression",
    "pol3_expression_backbone": "pol3_expression",
    "pol3_expression_generic": "pol3_expression",
}


def _types_compatible(start_type: str, end_type: str) -> bool:
    """
    True if start_type and end_type can be paired as module boundaries.
    Exact match always works; family match allows generic terminators/polyAs
    to close scope-specific starts (e.g., pol2_expression_animal ↔ pol2_expression_generic).
    """
    if start_type == end_type:
        return True
    sf = _MODULE_TYPE_FAMILY.get(start_type)
    ef = _MODULE_TYPE_FAMILY.get(end_type)
    return sf is not None and sf == ef


# ──────────────────────────────────────────────
# Canonical type mapping from functional bucket
# ──────────────────────────────────────────────

_BUCKET_TO_CTYPE: Dict[str, str] = {
    "pol2_promoter": "promoter",
    "pol3_promoter": "promoter",
    "pol2_enhancer": "enhancer",
    "bacterial_promoter": "promoter",
    "pol2_terminator": "polya",
    "pol3_terminator": "terminator",
    "generic_terminator": "terminator",
    "bacterial_terminator": "terminator",
    "replication_origin": "origin",
    "agrobacterium_replication": "origin",
    "lentiviral_boundary": "lenti_element",
    "lentiviral_packaging_signal": "lenti_element",
    "lentiviral_internal_element": "lenti_element",
    "tdna_border": "misc",
    "guide_scaffold": "scaffold",
    "guide_payload": "misc",
    "selectable_marker_cds": "marker",
    "coding_sequence": "cds",
    "translation_initiation": "rbs",
    "recombination_site": "misc",
    "binding_site": "misc",
}

# Grammar module types produced by boundary pairing
GRAMMAR_MODULE_PRIORITIES: Dict[str, int] = {
    "lentiviral_payload": 100,
    "tdna_payload": 98,
    "pol2_expression_animal": 87,
    "pol2_expression_plant": 86,
    "pol2_expression_fungal": 85,
    "pol2_expression_insect": 84,
    "pol2_expression_generic": 83,
    "pol3_expression_animal": 82,
    "pol3_expression_generic": 80,
    "pol3_expression_insect": 80,
    "bacterial_expression": 78,
}

_GRAMMAR_FAMILY: Dict[str, str] = {
    "pol2_expression_animal": "pol2_expression",
    "pol2_expression_plant": "pol2_expression",
    "pol2_expression_fungal": "pol2_expression",
    "pol2_expression_insect": "pol2_expression",
    "pol2_expression_generic": "pol2_expression",
    "pol3_expression_animal": "pol3_expression",
    "pol3_expression_generic": "pol3_expression",
    "pol3_expression_insect": "pol3_expression",
    "pol3_u6_sgrna_cassette": "pol3_expression",
    "pol3_expression_backbone": "pol3_expression",
    "lentiviral_payload": "lentiviral_backbone",
    "mammalian_replication_module": "mammalian_replication",
}

MODULE_COLORS = {
    "pol2_expression": "#D97706",
    "pol3_expression": "#B45309",
    "lentiviral_backbone": "#1D4ED8",
    "coding_part": "#BE185D",
    "regulatory_element": "#7C2D12",
    "mammalian_replication": "#A855F7",  # violet/purple for mammalian episomal origins
    "misc": "#4B5563",
}

MODULE_NESTING_RULES: Dict[str, List[str]] = {
    "shuttle_vector_backbone": ["yeast_pol2_expression_cassette"],
    "yeast_pol2_expression_cassette": ["cds_module"],
}

MODULE_PRIORITY: Dict[str, int] = {
    "mammalian_lentiviral_expression": 100,
    "mammalian_replication_module": 95,  # Higher than bacterial backbone, lower than lentiviral
    "lentiviral_expression_vector": 100,
    "shuttle_vector_backbone": 98,
    "yeast_integrating_marker": 90,
    "lentiviral_backbone": 90,
    "nuclease_expression_cassette": 85,
    "reporter_expression_cassette": 84,
    "pol2_expression_cassette": 83,
    "pol3_u6_sgrna_cassette": 82,
    "sv40_neo_selection_cassette": 81,
    "bacterial_expression_cassette": 80,
    "yeast_pol2_expression_cassette": 78,
    "bacterial_marker_cassette": 75,
    "lacI_regulator_module": 74,
    "pol3_expression_backbone": 73,
    "bacterial_origin": 60,
    **GRAMMAR_MODULE_PRIORITIES,
}
TOP_LEVEL_PRIORITIES = MODULE_PRIORITY


# ──────────────────────────────────────────────
# Canonical type / ID inference
# ──────────────────────────────────────────────

def _clean_text(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]+", " ", str(value or "").upper()).strip()


def _contains_any(text: str, needles: Iterable[str]) -> bool:
    return any(needle in text for needle in needles)


def _slug(text: str) -> str:
    norm = re.sub(r"[^A-Z0-9]+", "_", text).strip("_")
    return norm or "UNKNOWN"


def infer_canonical_type(feature_name: str, feature_type: str, description: str) -> Optional[str]:
    # Try rules lookup first (by feature name)
    rule = _lookup_rule("", feature_name)
    if rule and rule.functional_bucket in _BUCKET_TO_CTYPE and _BUCKET_TO_CTYPE[rule.functional_bucket] != "misc":
        return _BUCKET_TO_CTYPE[rule.functional_bucket]

    # Heuristic fallback
    text = " ".join(filter(None, [_clean_text(feature_name), _clean_text(feature_type), _clean_text(description)]))

    if "PROMOTER" in text:
        return "promoter"
    if _contains_any(text, ["POLYA", "POLY A", "POLYADENYLATION", " BGH PA", "SV40 PA", "HGH PA", "TK PA"]):
        return "polya"
    if _contains_any(text, ["SGRNA SCAFFOLD", "GRNA SCAFFOLD", "TRACRRNA", "SCAFFOLD"]):
        return "scaffold"
    if _contains_any(text, ["TERMINATOR", "POLY T", "POLYT"]):
        return "terminator"
    if _contains_any(text, [" ORI", "ORIGIN", "COLE1", "PBR322", "P15A", "SC101", "F1 ORIGIN"]):
        return "origin"
    if _contains_any(text, ["LTR", "WPRE", "RRE", "CPPT", "PSI PACKAGING", "PACKAGING SIGNAL"]):
        return "lenti_element"
    if "RBS" in text:
        return "rbs"
    if _contains_any(text, ["OPERATOR", "LACO", "LAC O"]):
        return "operator"
    if _contains_any(text, ["AMPICILLIN", "AMP R", "BLA", "KANAMYCIN", "NEOMYCIN", "PUROMYCIN",
                             "HYGROMYCIN", "ZEOCIN", "SPECTINOMYCIN", "CHLORAMPHENICOL",
                             "RESISTANCE", "SELECTABLE MARKER"]):
        return "marker"
    if _contains_any(text, ["EGFP", "GFP", "MCHERRY", "DSRED", "BFP", "RFP", "LUC2", "LUCIFERASE", "ZSGREEN"]):
        return "reporter"
    if _contains_any(text, ["CDS", "PROTEIN", "GENE", "CAS9", "CAS12", "CPF1"]):
        return "cds"
    return None


def infer_canonical_id(feature_name: str, feature_type: str, description: str, canonical_type: Optional[str]) -> Optional[str]:
    text = " ".join(filter(None, [_clean_text(feature_name), _clean_text(feature_type), _clean_text(description)]))

    if canonical_type == "promoter":
        if "DU6" in text or "DROSOPHILA" in text and "U6" in text:
            return "PROMOTER_DU6"
        if "U6" in text:
            return "PROMOTER_U6_HUMAN"
        if re.search(r"\bH1\b", text):
            return "PROMOTER_H1"
        if "7SK" in text:
            return "PROMOTER_7SK"
        if "CMV" in text:
            return "PROMOTER_CMV"
        if "EF1A" in text or "EF 1A" in text:
            return "PROMOTER_EF1A"
        if "PGK" in text:
            return "PROMOTER_PGK"
        if "CAG" in text:
            return "PROMOTER_CAG"
        if "SV40" in text:
            return "PROMOTER_SV40"
        if "TRE" in text or "TET RESPONSE" in text:
            return "PROMOTER_TRE"
        if "T7" in text:
            return "PROMOTER_T7"
        if "LACI" in text:
            return "PROMOTER_LACI"
        if "LAC" in text:
            return "PROMOTER_LAC"
        if _contains_any(text, ["AMP R PROMOTER", "BLA PROMOTER", "AMPICILLIN PROMOTER"]):
            return "PROMOTER_AMPR"
        return f"PROMOTER_{_slug(feature_name or description)}"

    if canonical_type == "polya":
        if "BGH" in text:
            return "POLYA_BGH"
        if "SV40" in text:
            return "POLYA_SV40"
        if "HGH" in text:
            return "POLYA_HGH"
        if re.search(r"\bTK\b", text):
            return "POLYA_TK"
        return f"POLYA_{_slug(feature_name or description)}"

    if canonical_type == "scaffold":
        return "GRNA_SCAFFOLD"

    if canonical_type == "terminator":
        return f"TERMINATOR_{_slug(feature_name or description)}"

    if canonical_type == "origin":
        if "PUC" in text:
            return "ORI_PUC"
        if "COLE1" in text:
            return "ORI_COLE1"
        if "PBR322" in text:
            return "ORI_PBR322"
        if "P15A" in text:
            return "ORI_P15A"
        if "SC101" in text:
            return "ORI_SC101"
        # Phage origins
        if re.search(r"\bM13\b", text) or "M13 ORI" in text:
            return "ORI_M13"
        if re.search(r"\bN15\b", text) or "N15 ORI" in text:
            return "ORI_N15"
        if re.search(r"\bF1\b", text) or "F1 ORI" in text:
            return "ORI_F1"
        return f"ORI_{_slug(feature_name or description)}"

    if canonical_type == "marker":
        if _contains_any(text, ["AMPICILLIN", "AMP R", "BLA"]):
            return "MARKER_AMP"
        if _contains_any(text, ["KANAMYCIN", "KAN R", "NEOMYCIN", "NEO R", "G418"]):
            return "MARKER_NEO"
        if "PUROMYCIN" in text:
            return "MARKER_PURO"
        if "HYGROMYCIN" in text:
            return "MARKER_HYGRO"
        if "ZEOCIN" in text:
            return "MARKER_ZEO"
        if "SPECTINOMYCIN" in text:
            return "MARKER_SPEC"
        if "CHLORAMPHENICOL" in text:
            return "MARKER_CM"
        return f"MARKER_{_slug(feature_name or description)}"

    if canonical_type == "lenti_element":
        if "WPRE" in text:
            return "LENTI_ELEMENT_WPRE"
        if "RRE" in text:
            return "LENTI_ELEMENT_RRE"
        if "CPPT" in text or "CTS" in text:
            return "LENTI_ELEMENT_CPPT"
        if "PSI" in text or "PACKAGING SIGNAL" in text:
            return "LENTI_ELEMENT_PSI"
        if "LTR" in text:
            return "LENTI_ELEMENT_LTR"
        return f"LENTI_ELEMENT_{_slug(feature_name or description)}"

    if "KOZAK" in text:
        return "MISC_KOZAK"
    if "P2A" in text:
        return "CDS_P2A"
    if "T2A" in text:
        return "CDS_T2A"
    if "E2A" in text:
        return "CDS_E2A"
    if "F2A" in text:
        return "CDS_F2A"
    if "IRES" in text:
        return "MISC_IRES"
    if _contains_any(text, ["SV40 NLS", "NUCLEAR LOCALIZATION SIGNAL", "NLS"]):
        if "NUCLEOPLASMIN" in text:
            return "CDS_NUCLEOPLASMIN_NLS"
        return "CDS_SV40_NLS"

    if canonical_type == "reporter":
        if "EGFP" in text or re.search(r"\bGFP\b", text):
            return "CDS_EGFP"
        if "MCHERRY" in text:
            return "CDS_MCHERRY"
        if "DSRED" in text:
            return "REPORTER_DSRED"
        if "BFP" in text:
            return "REPORTER_BFP"
        if "LUC2" in text or "LUCIFERASE" in text:
            return "REPORTER_LUC2"
        return f"REPORTER_{_slug(feature_name or description)}"

    if canonical_type == "cds":
        if "CAS9" in text:
            return "CDS_CAS9"
        if "DCAS9" in text:
            return "CDS_DCAS9"
        if "CAS12" in text or "CPF1" in text:
            return "CDS_CPF1"
        if re.search(r"\bLACI\b", text):
            return "CDS_LACI"
        if re.search(r"\bROP\b", text):
            return "GENE_ROP"
        return f"CDS_{_slug(feature_name or description)}"

    if canonical_type == "rbs":
        return "RBS"
    if canonical_type == "operator":
        return "LAC_OPERATOR"
    return None


# ──────────────────────────────────────────────
# Feature construction from pLannotate rows
# ──────────────────────────────────────────────

def _feature_from_plannotate_row(row: Dict[str, Any]) -> Optional[Feature]:
    try:
        qstart = int(row.get("qstart", 1))
        qend = int(row.get("qend", 0))
    except Exception:
        return None

    if qend <= 0 or qend <= qstart - 1:
        return None

    raw_frame = row.get("sframe", 1)
    try:
        strand = 1 if int(raw_frame) > 0 else -1
    except Exception:
        strand = 1

    sseqid = str(row.get("Feature", "") or "")  # pLannotate's Feature column is the sseqid
    feature_type = str(row.get("Type", "misc_feature") or "misc_feature")
    description = str(row.get("Description", "") or "")

    # KB lookup for better feature naming
    kb_info = KnowledgeBase.lookup(sseqid, description)
    if kb_info:
        # Use KB feature_name for better readability
        feature_name = kb_info.get('feature_name', sseqid)
        # Get KB descriptions if available
        if 'source' in kb_info and 'descriptions' in kb_info['source']:
            kb_descriptions = kb_info['source']['descriptions']
            if kb_descriptions:
                description = '; '.join(kb_descriptions)
    else:
        feature_name = sseqid

    # Prefer rule-based canonical type when sseqid matches
    rule = _lookup_rule(sseqid if sseqid else feature_name, feature_name)
    if rule and rule.functional_bucket in _BUCKET_TO_CTYPE and _BUCKET_TO_CTYPE[rule.functional_bucket] != "misc":
        canonical_type = _BUCKET_TO_CTYPE[rule.functional_bucket]
    else:
        canonical_type = infer_canonical_type(feature_name, feature_type, description)

    canonical_id = infer_canonical_id(feature_name, feature_type, description, canonical_type)
    normalized_text = " ".join(filter(None, [feature_name.strip().lower(), description.strip().lower()])) or None

    return Feature(
        type=feature_type,
        name=feature_name or None,
        start=max(0, qstart - 1),
        end=max(0, qend),
        strand=strand,
        canonical_id=canonical_id,
        canonical_type=canonical_type,
        normalized_text=normalized_text,
        qualifiers={
            "label": feature_name,
            "note": description,
            "source": "pLannotate",
            "sseqid": sseqid,
        },
    )


def infer_features_from_plannotate(rows: Sequence[Dict[str, Any]]) -> List[Feature]:
    features: List[Feature] = []
    for row in rows:
        feat = _feature_from_plannotate_row(row)
        if feat is not None:
            features.append(feat)
    return features


# ──────────────────────────────────────────────
# Grammar-driven module construction
# ──────────────────────────────────────────────

def _downstream_distance(
    from_feat: Feature,
    to_feat: Feature,
    strand: int,
    seq_len: int,
    circular: bool,
) -> Optional[int]:
    """
    Distance from from_feat to to_feat going in `strand` direction.
    Returns a non-negative integer, or None if to_feat is not downstream.
    """
    if strand >= 0:  # +1 strand: downstream = increasing coords
        direct = to_feat.start - from_feat.end
        if direct >= 0:
            return direct
        if circular:
            return seq_len - from_feat.end + to_feat.start
        return None
    else:  # -1 strand: downstream = decreasing coords
        direct = from_feat.start - to_feat.end
        if direct >= 0:
            return direct
        if circular:
            return from_feat.start + (seq_len - to_feat.end)
        return None


def _compute_module_span(
    start_feat: Feature,
    end_feat: Feature,
    strand: int,
    seq_len: int,
    circular: bool,
) -> Optional[Tuple[int, int, bool]]:
    """
    Returns (mod_start, mod_end, wraps) for the module bounded by start and end features.
    mod_start > mod_end iff wraps=True (module crosses position 0).
    Returns None if the boundary order is incompatible.
    """
    if strand >= 0:
        # +1 strand: start_feat is upstream (lower coords), end_feat is downstream (higher coords)
        if end_feat.start >= start_feat.end:
            return start_feat.start, end_feat.end, False
        elif circular and end_feat.end <= start_feat.start:
            # end_feat is to the left; module wraps around origin
            return start_feat.start, end_feat.end, True  # start > end → wraps
        return None
    else:
        # -1 strand: start_feat is upstream (higher coords), end_feat is downstream (lower coords)
        if end_feat.end <= start_feat.start:
            # Check whether wrapping via the origin is shorter than going directly left
            left_dist = start_feat.start - end_feat.end
            if circular:
                wrap_dist = start_feat.start + (seq_len - end_feat.end)
                if wrap_dist < left_dist:
                    # Wrap is shorter: module goes left from start_feat past origin to end_feat
                    return start_feat.start, end_feat.end, True
            return end_feat.start, start_feat.end, False
        elif circular and end_feat.start >= start_feat.end:
            # end_feat is to the right of start_feat; module wraps the other way
            return end_feat.start, start_feat.end, True
        return None


def _infer_payload_id(features: List[Feature]) -> Optional[str]:
    """Pick a payload CDS canonical_id from the features in a module."""
    preferred = [
        f for f in features
        if f.canonical_type == "reporter" and f.canonical_id
    ]
    if preferred:
        return preferred[0].canonical_id
    payload = [
        f for f in features
        if f.canonical_type == "cds"
        and f.canonical_id
        and not (f.canonical_id or "").startswith("MARKER_")
    ]
    if payload:
        return payload[0].canonical_id
    return None


def _build_grammar_modules(
    sequence: str,
    circular: bool,
    features: Sequence[Feature],
    plasmid_id: str,
) -> List[Module]:
    """
    Construct modules by pairing hard_start / fallback_start features with compatible
    hard_end features, using the hierarchical_feature_rules grammar sheet.
    """
    seq_len = len(sequence)
    features = list(features)

    # Attach rules
    enriched: List[Tuple[Feature, Optional[FeatureRule]]] = [
        (f, _get_rule_for_feature(f)) for f in features
    ]

    # Partition into starts and ends
    starts: List[Tuple[Feature, FeatureRule]] = [
        (f, r) for f, r in enriched
        if r
        and r.nested_module_role in ("hard_start", "fallback_start")
        and r.starts_module_types
    ]
    ends: List[Tuple[Feature, FeatureRule]] = [
        (f, r) for f, r in enriched
        if r
        and r.nested_module_role == "hard_end"
        and r.ends_module_types
    ]

    if not starts or not ends:
        return []

    # Sort starts: hard_start before fallback_start, then higher priority first, then position
    starts_sorted = sorted(
        starts,
        key=lambda x: (
            0 if x[1].nested_module_role == "hard_start" else 1,
            -x[1].boundary_priority,
            x[0].start,
        ),
    )

    modules: List[Module] = []
    # Track which (start_feat, module_type) pairs have already been consumed
    used: set = set()

    for start_feat, start_rule in starts_sorted:
        strand = start_feat.strand if start_feat.strand != 0 else 1

        for module_type in start_rule.starts_module_types:
            key = (id(start_feat), module_type)
            if key in used:
                continue

            # Compatible ends: matching module type (exact or same family), same strand
            compatible_ends = [
                (f, r) for f, r in ends
                if any(_types_compatible(module_type, et) for et in r.ends_module_types)
                and (f.strand == strand or f.strand == 0 or strand == 0)
            ]
            if not compatible_ends:
                continue

            # Find nearest downstream end
            best: Optional[Tuple[Feature, FeatureRule]] = None
            best_dist: float = float("inf")
            for ef, er in compatible_ends:
                dist = _downstream_distance(start_feat, ef, strand, seq_len, circular)
                if dist is not None and 0 < dist < best_dist:
                    best_dist = dist
                    best = (ef, er)

            if best is None:
                continue

            end_feat, _ = best
            span = _compute_module_span(start_feat, end_feat, strand, seq_len, circular)
            if span is None:
                continue

            mod_start, mod_end, wraps = span
            mod_len = module_length(seq_len, mod_start, mod_end, wraps)

            # Reject implausibly large or tiny modules
            if mod_len > 0.95 * seq_len or mod_len < 50:
                continue

            module_seq = extract_sequence(sequence, mod_start, mod_end, wraps)
            enclosed = [f for f in features if contains_feature(seq_len, mod_start, mod_end, wraps, f)]
            payload_id = _infer_payload_id(enclosed)
            family = _GRAMMAR_FAMILY.get(module_type, "misc")

            mod_id = sha256_text(
                f"{plasmid_id}|{module_type}|{mod_start}|{mod_end}|{payload_id or ''}"
            )[:24]

            mod = Module(
                id=mod_id,
                plasmid_id=plasmid_id,
                module_type=module_type,
                payload_id=payload_id,
                start=mod_start,
                end=mod_end,
                wraps=wraps,
                length=mod_len,
                sequence=module_seq,
                seq_hash=seq_hash(module_seq),
                end_inferred=False,
                metadata={
                    "module_family": family,
                    "taxonomic_scope": start_rule.taxonomic_scope,
                    "boundary_start": start_feat.qualifiers.get("sseqid", start_feat.name or ""),
                    "boundary_end": end_feat.qualifiers.get("sseqid", end_feat.name or ""),
                    "source": "grammar",
                },
                features=enclosed,
            )
            modules.append(mod)
            used.add(key)
            break  # one module per (start_feat, strand) pass

    return modules


# ──────────────────────────────────────────────
# Module extraction: grammar-first, legacy fallback
# ──────────────────────────────────────────────

def extract_modules_from_features(
    sequence: str,
    circular: bool,
    features: Sequence[Feature],
    plasmid_id: Optional[str] = None
) -> List[Module]:
    seq = (sequence or "").upper()
    seq_len = len(seq)
    topology = "circular" if circular else "linear"
    pid = plasmid_id or sha256_text(seq)[:24]

    # ── 1. Grammar-driven modules (pol2, pol3, lentiviral, T-DNA) ──
    grammar_modules = _build_grammar_modules(seq, circular, features, pid)

    # Post-process grammar modules: truncate pol2/pol3 modules that extend past lentiviral_payload
    # Find all lentiviral_payload modules
    lentiviral_payloads = [m for m in grammar_modules if m.module_type == "lentiviral_payload"]
    
    if lentiviral_payloads:
        # Check each pol2/pol3 module
        truncated_modules = []
        for mod in grammar_modules:
            if mod.module_type in ("pol2_expression_animal", "pol2_expression_generic", "pol2_expression_plant",
                                  "pol3_expression_backbone", "pol3_u6_sgrna_cassette"):
                # Check if this module starts inside any payload but extends past it
                for payload in lentiviral_payloads:
                    if payload.wraps or mod.wraps:
                        # Skip wrapping case for now
                        truncated_modules.append(mod)
                        break
                    elif mod.start >= payload.start and mod.start < payload.end and mod.end > payload.end:
                        # Module starts inside payload but extends past it - truncate to payload end
                        # Create a new module with truncated end
                        from .Module_Library_gb.module_extractor import extract_sequence
                        new_end = payload.end
                        new_seq = extract_sequence(seq, mod.start, new_end, False)
                        truncated_mod = Module(
                            id=sha256_text(f"{pid}|{mod.module_type}|{mod.start}|{new_end}")[:24],
                            plasmid_id=pid,
                            module_type=mod.module_type,
                            payload_id=mod.payload_id,
                            start=mod.start,
                            end=new_end,
                            wraps=False,
                            length=len(new_seq),
                            sequence=new_seq,
                            seq_hash=mod.seq_hash,
                            end_inferred=True,
                            metadata=mod.metadata.copy(),
                            features=[f for f in mod.features if f.start >= mod.start and f.end <= new_end],
                        )
                        truncated_modules.append(truncated_mod)
                        break
                else:
                    # No truncation needed
                    truncated_modules.append(mod)
            else:
                # Not a pol2/pol3 module, keep as is
                truncated_modules.append(mod)
        
        grammar_modules = truncated_modules

    # Filter out pol2_expression_generic modules that encounter an LTR before finding a CDS
    # (These are likely misdetected - the promoter is being used for lentiviral expression, not generic expression)
    filtered_modules = []
    for mod in grammar_modules:
        if mod.module_type == "pol2_expression_generic":
            # Sort features by position to check order correctly
            sorted_features = sorted(mod.features, key=lambda f: f.start)
            
            # Check if there's an LTR feature before the first CDS
            has_ltr_before_cds = False
            has_cds = False
            
            for f in sorted_features:
                # Check if this is an LTR - look for multiple indicators
                # LTRs may have canonical_type "ltr", "lenti_element", or just "ltr" in the name/sseqid
                sseqid = (f.qualifiers or {}).get("sseqid", "")
                name_lower = (f.name or "").lower()
                canonical_id_lower = (f.canonical_id or "").lower()
                
                is_ltr = (
                    (f.canonical_type == "ltr") or
                    (f.canonical_type == "lenti_element" and "ltr" in name_lower) or
                    ("ltr" in name_lower) or
                    ("ltr" in canonical_id_lower) or
                    ("ltr" in sseqid.lower())
                )
                
                # Check if this is a CDS
                is_cds = (
                    (f.canonical_type == "cds") or
                    (f.canonical_type == "protein") or
                    (f.canonical_type == "reporter") or
                    ("cds" in canonical_id_lower)
                )
                
                if is_cds:
                    has_cds = True
                
                # If we find an LTR before finding a CDS, mark as invalid
                if is_ltr and not has_cds:
                    has_ltr_before_cds = True
                    break
            
            # Skip this module if it has LTR before CDS
            if has_ltr_before_cds:
                continue
        
        filtered_modules.append(mod)
    
    grammar_modules = filtered_modules

    # ── 2. Legacy modules for backbone / origin / marker / misc ──
    legacy_modules: List[Module] = []

    # Module detection: Always use OLD pipeline for now
    # (New ModuleDetectionPipeline needs more work for production use)
    if False:  # Disabled: use_new_pipeline for module detection
        # NEW PIPELINE - Use ModuleDetectionPipeline
        from .annotation_pipeline import (
            ModuleDetectionPipeline,
            OriginDetector,
            MarkerDetector,
            BackboneDetector
        )

        # Convert features to dict format for new pipeline
        feature_dicts = [{'id': f.canonical_id or f.name or str(i), 
                         'canonical_type': f.canonical_type or f.type,
                         'name': f.name or '', 
                         'start': f.start, 
                         'end': f.end,
                         'strand': f.strand} for i, f in enumerate(features)]

        pipeline = ModuleDetectionPipeline()
        pipeline.register_detector(OriginDetector())
        pipeline.register_detector(MarkerDetector())
        pipeline.register_detector(BackboneDetector())

        module_candidates = pipeline.detect_all(feature_dicts, seq)

        # Convert module dicts back to Module objects
        # (sha256_text and seq_hash are already imported at top of file)
        for i, mod_dict in enumerate(module_candidates):
            mod_start = mod_dict['start']
            mod_end = mod_dict['end']
            module_seq = seq[mod_start:mod_end] if mod_end <= len(seq) else seq[mod_start:] + seq[:mod_end - len(seq)]
            mod_id = f"{pid}_new_pipeline_mod_{i}"
            
            mod = Module(
                id=mod_id,
                plasmid_id=pid,
                module_type=mod_dict['module_type'],
                payload_id=None,
                start=mod_start,
                end=mod_end,
                wraps=False,
                length=mod_end - mod_start,
                sequence=module_seq,
                seq_hash=seq_hash(module_seq),
                end_inferred=False,
                features=[],  # Features will be populated by subsequent processing
                metadata={
                    **mod_dict.get('metadata', {}),
                    'strand': mod_dict.get('strand', 1)  # Store strand in metadata
                }
            )
            legacy_modules.append(mod)
    else:
        # OLD PIPELINE - Use choose_candidates_* functions
        legacy_modules.extend(choose_candidates_regulatory_elements(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_linker_and_targeting(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_marker(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_origin(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_bacterial_expression(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_bacterial_marker_cassette(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_yeast_pol2(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_shuttle_vector(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_yeast_integrating(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_sv40_neo(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_mammalian_replication(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_lacI_regulator(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_mammalian_lentiviral_expression(list(features), seq_len, topology, pid, seq))
        legacy_modules.extend(choose_candidates_cds_only(list(features), seq_len, topology, pid, seq))

    # ── 3. Lentiviral wrapper: grammar lentiviral_payload + legacy pol2 cassettes ──
    # Build grammar lentiviral wrappers from lentiviral_payload + pol2 grammar modules
    lenti_grammar = [m for m in grammar_modules if m.module_type == "lentiviral_payload"]
    pol2_grammar = [m for m in grammar_modules if m.module_type.startswith("pol2_")]
    # Also keep legacy pol2 for lenti wrapping if grammar pol2 is empty
    if not pol2_grammar:
        pol2_grammar = choose_candidates_pol2(list(features), seq_len, topology, pid, seq)

    # Skip creating lentiviral_expression_vector if mammalian_lentiviral_expression exists
    has_mammalian_lenti = any(m.module_type == "mammalian_lentiviral_expression" for m in legacy_modules)
    if has_mammalian_lenti:
        lenti_wrappers = []
    else:
        lenti_wrappers: List[Module] = []

        for lenti in lenti_grammar:
            for pol2 in pol2_grammar:
                start = min(lenti.start, pol2.start)
                end = max(lenti.end, pol2.end)
                wraps = start > end
                if wraps and not circular:
                    continue
                span = module_length(seq_len, start, end, wraps)
                if span > 0.9 * seq_len:
                    continue
                module_seq = extract_sequence(seq, start, end, wraps)
                lenti_wrappers.append(
                    Module(
                        id=sha256_text(f"{pid}|lentiviral_expression_vector|{start}|{end}|{pol2.payload_id or ''}")[:24],
                        plasmid_id=pid,
                        module_type="lentiviral_expression_vector",
                        payload_id=pol2.payload_id,
                        start=start,
                        end=end,
                        wraps=wraps,
                        length=len(module_seq),
                        sequence=module_seq,
                        seq_hash=seq_hash(module_seq),
                        end_inferred=False,
                        metadata={
                        "module_family": "lentiviral_backbone",
                        "lenti_module_id": lenti.id,
                        "expression_module_id": pol2.id,
                        "payload_family": pol2.metadata.get("payload_family"),
                    },
                        features=[f for f in features if contains_feature(seq_len, start, end, wraps, f)],
                    )
                )

    # ── 4. Merge: grammar takes precedence; legacy fills gaps ──
    all_modules = grammar_modules + lenti_wrappers + legacy_modules

    # Get expression cassette positions to filter overlapping cds_only modules
    expression_types = {
        "pol2_expression_cassette", "nuclease_expression_cassette",
        "reporter_expression_cassette", "pol2_expression_animal",
        "pol2_expression_generic", "pol2_expression_plant",
        "pol2_expression_fungal", "pol2_expression_insect",
        "bacterial_expression_cassette", "yeast_pol2_expression_cassette",
    }
    expression_ranges = []
    for mod in all_modules:
        if mod.module_type in expression_types:
            if mod.wraps:
                expression_ranges.append((mod.start, seq_len))
                expression_ranges.append((0, mod.end))
            else:
                expression_ranges.append((mod.start, mod.end))

    def is_cds_contained_in_cassette(cds_mod):
        """Check if a cds_only module is contained within an expression cassette."""
        if cds_mod.module_type != "cds_only":
            return False
        for exp_start, exp_end in expression_ranges:
            if cds_mod.start >= exp_start and cds_mod.end <= exp_end:
                return True
        return False

    # Filter out cds_only modules contained within expression cassettes
    # (they will be detected as cds_modules instead)
    filtered_modules = [
        mod for mod in all_modules
        if not is_cds_contained_in_cassette(mod)
    ]

    unique: Dict[str, Module] = {}
    for mod in filtered_modules:
        if not mod.features:
            mod.features = [f for f in features if contains_feature(seq_len, mod.start, mod.end, mod.wraps, f)]
        _enrich_module_metadata(mod)
        if mod.id not in unique:
            unique[mod.id] = mod

    return list(unique.values())


# ──────────────────────────────────────────────
# Annotation helpers
# ──────────────────────────────────────────────

def _module_segments(module: Module, seq_len: int) -> List[Tuple[int, int]]:
    if not module.wraps:
        return [(module.start, module.end)]
    return [(module.start, seq_len), (0, module.end)]


def _dominant_direction(module: Module) -> int:
    score = sum(f.strand for f in module.features if f.strand in (-1, 1))
    if score < 0:
        return -1
    if score > 0:
        return 1
    return 0


def _module_label(module: Module) -> str:
    payload = module.payload_id or ""
    base = module.module_type.replace("_", " ")
    if payload:
        return f"{base} ({payload})"
    return base


def _module_to_annotation(module: Module, seq_len: int) -> List[Dict[str, Any]]:
    family = module.metadata.get("module_family", "misc")
    color = _module_color(module)
    direction = _dominant_direction(module)
    annotations: List[Dict[str, Any]] = []
    for start, end in _module_segments(module, seq_len):
        if end <= start:
            continue
        annotations.append(
            {
                "name": _module_label(module),
                "start": start,
                "end": end,
                "direction": direction,
                "color": color,
                "source": "hierarchical_annotator",
                "layer": "module",
                "module_type": module.module_type,
                "module_family": family,
                "payload_id": module.payload_id,
                "metadata": module.metadata,
            }
        )
    return annotations


def _module_priority(module: Module) -> Tuple[int, int]:
    return (MODULE_PRIORITY.get(module.module_type, 0), module.length)




def _motif_color(motif_type: str) -> str:
    """Return color for different motif types"""
    colors = {
        "start_codon": "#22c55e",    # Green
        "stop_codon": "#ef4444",     # Red
        "2a_peptide": "#f59e0b",     # Amber
        "kozak": "#3b82f6",          # Blue
        "internal_stop": "#dc2626",  # Dark red
        "intron": "#8b5cf6",         # Purple
    }
    return colors.get(motif_type, "#6b7280")  # Gray default

def _module_color(module: Module) -> str:
    type_colors = {
        "shuttle_vector_backbone": "#6366F1",
        "yeast_integrating_marker": "#B45309",
        "yeast_pol2_expression_cassette": "#CA8A04",
    }
    if module.module_type in type_colors:
        return type_colors[module.module_type]
    family = module.metadata.get("module_family", "misc")
    return MODULE_COLORS.get(family, MODULE_COLORS["misc"])



# Container modules that can fully overlap with their children
CONTAINER_MODULES = {
    "mammalian_lentiviral_expression",
    "mammalian_lentiviral_expression_module",
    "lentiviral_expression_vector",
    "shuttle_vector_backbone",
}

def _calculate_overlap_fraction(mod1, mod2, seq_len):
    """Calculate the fraction of overlap between two modules (0.0 to 1.0)."""
    if mod1.wraps or mod2.wraps:
        return 0.0
    overlap_start = max(mod1.start, mod2.start)
    overlap_end = min(mod1.end, mod2.end)
    if overlap_end <= overlap_start:
        return 0.0
    overlap_len = overlap_end - overlap_start
    smaller_len = min(mod1.end - mod1.start, mod2.end - mod2.start)
    return overlap_len / smaller_len if smaller_len > 0 else 0.0

def select_top_level_modules(modules: Sequence[Module], seq_len: int) -> List[Module]:
    """
    Select non-overlapping top-level modules, preferring higher priority modules.
    Uses position-based overlap detection for proper circular handling.
    """
    selected: List[Module] = []
    occupied_positions: set = set()

    candidates = [
        mod
        for mod in modules
        if MODULE_PRIORITY.get(mod.module_type, 0) > 0
    ]
    candidates.sort(key=_module_priority, reverse=True)

    for mod in candidates:
        segments = _module_segments(mod, seq_len)
        mod_positions = set()
        for seg_start, seg_end in segments:
            if seg_end > seg_start:
                for pos in range(seg_start, seg_end):
                    mod_positions.add(pos)
        
        # Check overlap - allow small overlaps (< 10% of module)
        overlap = mod_positions & occupied_positions
        overlap_fraction = len(overlap) / len(mod_positions) if mod_positions else 0
        
        # Container modules are allowed to fully overlap with their children
        if mod.module_type in CONTAINER_MODULES:
            pass
        elif overlap_fraction > 0.1:  # More than 10% overlap
            continue
            continue
        
        selected.append(mod)
        occupied_positions |= mod_positions

    return sorted(selected, key=lambda mod: mod.start)


def gap_annotations_for_modules(modules: Sequence[Module], seq_len: int) -> List[Dict[str, Any]]:
    """
    Generate inter-module gap annotations, properly handling origin-crossing modules.
    """
    # Build a set of all occupied positions for accurate gap detection
    occupied_positions = set()
    for mod in modules:
        for seg_start, seg_end in _module_segments(mod, seq_len):
            if seg_end > seg_start:
                for pos in range(seg_start, seg_end):
                    occupied_positions.add(pos)
    
    # Find contiguous gaps
    gaps: List[Dict[str, Any]] = []
    gap_start = None
    gap_index = 1
    
    for pos in range(seq_len):
        if pos not in occupied_positions:
            if gap_start is None:
                gap_start = pos
        else:
            if gap_start is not None:
                # End of gap
                gaps.append(
                    {
                        "name": f"inter-module gap {gap_index}",
                        "start": gap_start,
                        "end": pos,
                        "direction": 0,
                        "color": "#9CA3AF",
                        "source": "hierarchical_annotator",
                        "layer": "gap",
                        "module_type": "inter_module_gap",
                    }
                )
                gap_index += 1
                gap_start = None
    
    # Handle trailing gap
    if gap_start is not None:
        gaps.append(
            {
                "name": f"inter-module gap {gap_index}",
                "start": gap_start,
                "end": seq_len,
                "direction": 0,
                "color": "#9CA3AF",
                "source": "hierarchical_annotator",
                "layer": "gap",
                "module_type": "inter_module_gap",
            }
        )

    return gaps


# ──────────────────────────────────────────────
# Public entry point
# ──────────────────────────────────────────────

def _filter_motifs_by_modules(
    motifs: List[Dict[str, Any]],
    modules: List[Dict[str, Any]],
    features: List[Dict[str, Any]],
    sequence: str
) -> List[Dict[str, Any]]:
    """Filter motifs to only those at actual CDS boundaries."""
    filtered = []
    
    # Get CDS features
    cds_features = [
        f for f in features
        if f.get("feature_type", "").lower() == "cds" or
           f.get("role", "") in ("expression_payload", "editing_payload", 
                                  "reporter_payload", "selection_payload", "bacterial_marker")
    ]
    
    for motif in motifs:
        motif_type = motif.get("motif_type", "")
        motif_start = motif.get("start", 0)
        motif_end = motif.get("end", 0)
        strand = motif.get("direction", 1)
        
        keep = False
        
        if motif_type == "start_codon":
            # Only keep if at actual start of a CDS
            for cds in cds_features:
                cds_start = cds.get("start", 0)
                cds_strand = cds.get("strand", 1)
                if cds_strand == 1 and abs(motif_start - cds_start) <= 3 and strand == cds_strand:
                    keep = True
                    break
                elif cds_strand == -1:
                    cds_end = cds.get("end", 0)
                    if abs(motif_end - cds_end) <= 3 and strand == cds_strand:
                        keep = True
                        break
        
        elif motif_type == "stop_codon":
            for cds in cds_features:
                cds_end = cds.get("end", 0)
                cds_strand = cds.get("strand", 1)
                if cds_strand == 1 and abs(motif_end - cds_end) <= 3 and strand == cds_strand:
                    keep = True
                    break
                elif cds_strand == -1:
                    cds_start = cds.get("start", 0)
                    if abs(motif_start - cds_start) <= 3 and strand == cds_strand:
                        keep = True
                        break
        
        elif motif_type == "kozak":
            for cds in cds_features:
                role = cds.get("role", "")
                if role == "bacterial_marker":
                    continue
                cds_start = cds.get("start", 0)
                cds_strand = cds.get("strand", 1)
                if cds_strand == 1 and abs(motif_start - (cds_start - 6)) <= 5:
                    keep = True
                    break
                elif cds_strand == -1:
                    cds_end = cds.get("end", 0)
                    if abs(motif_end - (cds_end + 6)) <= 5:
                        keep = True
                        break
        
        else:
            keep = True
        
        if keep:
            filtered.append(motif)
    
    return filtered



def annotate_hierarchy_from_plannotate(
    sequence: str,
    circular: bool,
    plannotate_rows: Sequence[Dict[str, Any]],
) -> Dict[str, Any]:
    features = infer_features_from_plannotate(plannotate_rows)
    modules = extract_modules_from_features(sequence=sequence, circular=circular, features=features)
    top_level = select_top_level_modules(modules, len(sequence))
    module_annotations: List[Dict[str, Any]] = []
    for mod in top_level:
        module_annotations.extend(_module_to_annotation(mod, len(sequence)))
    gap_annotations = gap_annotations_for_modules(top_level, len(sequence))

    # Detect motifs (start/stop codons, 2A peptides, Kozak sequences)
    motif_annotations: List[Dict[str, Any]] = []
    try:
        # Convert Feature objects to dicts with role for motif detection
        feature_dicts = []
        feature_instances = analyzer.get("feature_instances", [])
        for fi in feature_instances:
            role = fi.get("role", "") if isinstance(fi, dict) else getattr(fi, "role", "")
            ftype = fi.get("feature_type", "") if isinstance(fi, dict) else getattr(fi, "feature_type", "")
            start = fi.get("start", 0) if isinstance(fi, dict) else getattr(fi, "start", 0)
            end = fi.get("end", 0) if isinstance(fi, dict) else getattr(fi, "end", 0)
            strand = fi.get("strand", 1) if isinstance(fi, dict) else getattr(fi, "strand", 1)
            fname = fi.get("feature_name", "") if isinstance(fi, dict) else getattr(fi, "feature_name", "")
            feature_dicts.append({
                "role": role,
                "feature_type": ftype,
                "start": start,
                "end": end,
                "strand": strand,
                "feature_name": fname,
            })
        motif_detector = MotifDetector(sequence, circular)
        detected_motifs = motif_detector.detect_all_motifs(feature_dicts)
        for motif in detected_motifs:
            motif_annotations.append({
                "name": motif.name,
                "start": motif.start,
                "end": motif.end,
                "direction": motif.strand,
                "color": _motif_color(motif.motif_type),
                "source": "motif_detector",
                "layer": "motif",
                "type": "motif",
                "motif_type": motif.motif_type,
                "description": motif.description,
                "sequence": motif.sequence,
            })
    except Exception as e:
        print(f"[WARN] Motif detection failed: {e}")

    # Filter motifs to only those at actual CDS boundaries
    feature_dicts_for_filter = [{"feature_type": getattr(f, "type", "") or "", "role": getattr(f, "role", "") if hasattr(f, "role") else "", "start": getattr(f, "start", 0), "end": getattr(f, "end", 0), "strand": getattr(f, "strand", 1)} for f in features]
    module_dicts_for_filter = [{"module_type": m.module_type, "start": m.start, "end": m.end} for m in top_level]
    motif_annotations = _filter_motifs_by_modules(motif_annotations, module_dicts_for_filter, feature_dicts_for_filter, sequence)

    return {
        "features_inferred": len(features),
        "motif_annotations": motif_annotations,
        "module_count": len(top_level),
        "modules": [
            {
                "id": mod.id,
                "module_type": mod.module_type,
                "start": mod.start,
                "end": mod.end,
                "wraps": mod.wraps,
                "payload_id": mod.payload_id,
                "metadata": mod.metadata,
            }
            for mod in top_level
        ],
        "module_annotations": module_annotations,
        "feature_adjustments": feature_adjustments,
        "deleted_features": deleted_features,
        "gap_annotations": gap_annotations,
    }


# ──────────────────────────────────────────────
# Enhanced imports for new functionality
# ──────────────────────────────────────────────

from .Module_Library_gb.module_extractor import (
    choose_candidates_pol3_with_terminator_detection,
    detect_cds_modules,
    detect_bacterial_cds_module,
    is_pol3_promoter,
    adjust_feature_boundaries_to_submodules,
)


# ──────────────────────────────────────────────
# Enhanced Pol III promoter recognition
# ──────────────────────────────────────────────

def _is_pol3_promoter_id(canonical_id: str, feature_name: str) -> bool:
    """Check if this is a Pol III promoter based on ID or name."""
    cid = (canonical_id or "").upper()
    name = (feature_name or "").upper()
    
    pol3_markers = ["U6", "H1", "7SK", "DU6"]
    for marker in pol3_markers:
        if marker in cid or marker in name:
            return True
    return False


# ──────────────────────────────────────────────
# Enhanced module extraction with CDS detection
# ──────────────────────────────────────────────

def detect_lentiviral_cis_module(
    features: Sequence[Feature],
    seq_len: int,
    topology: str,
    plasmid_id: str
) -> List[Module]:
    """
    Detect lentiviral cis-acting element modules.
    Groups all lentiviral_cis_rna features into a single module spanning from
    first to last cis element (excluding LTRs).
    """
    modules = []
    
    # Find all lentiviral cis-acting RNA elements (excluding LTRs)
    cis_elements = [
        f for f in features 
        if f.canonical_type == 'regulatory' 
        and any(x in (f.name or '').lower() for x in [
            'psi', 'rre', 'cppt', 'cts', 'wpre', 'tar', 'rnai', 'dis', 'sd', 'sa'
        ])
        and 'ltr' not in (f.name or '').lower()
    ]
    
    if len(cis_elements) >= 2:
        # Sort by position
        sorted_cis = sorted(cis_elements, key=lambda f: f.start)
        first_element = sorted_cis[0]
        last_element = sorted_cis[-1]
        
        # Create module spanning from first to last cis element
        module = Module(
            module_id=f"mod_lvcis_{plasmid_id[:8]}_{first_element.start}",
            module_type="lentiviral_cis_module",
            start=first_element.start,
            end=last_element.end,
            features=cis_elements,
            feature_ids=[f.canonical_id for f in cis_elements if f.canonical_id],
            strand=first_element.strand,
            wraps=False,
            circular=topology == "circular",
            metadata={
                "cis_element_count": len(cis_elements),
                "elements": [f.name for f in sorted_cis],
                "function": "Lentiviral packaging and reverse transcription signals"
            }
        )
        modules.append(module)
    
    return modules


def extract_modules_with_cds_detection(
    sequence: str,
    circular: bool,
    features: Sequence[Feature],
    plasmid_id: Optional[str] = None
) -> Tuple[List[Module], List[Module]]:
    """
    Extract modules and detect CDS modules within expression cassettes.
    
    Returns:
        (all_modules, cds_modules)
    """
    # Get base modules
    modules = extract_modules_from_features(
        sequence=sequence,
        circular=circular,
        features=features,
        plasmid_id=plasmid_id
    )
    
    # Replace pol3 modules with enhanced detection
    non_pol3_modules = [m for m in modules if not m.module_type.startswith("pol3_")]
    
    seq_len = len(sequence)
    topology = "circular" if circular else "linear"
    pid = plasmid_id or sha256_text(sequence)[:24]
    
    pol3_modules = choose_candidates_pol3_with_terminator_detection(
        list(features), seq_len, topology, pid, sequence
    )
    
    # Detect lentiviral cis-acting element modules
    lentiviral_cis_modules = detect_lentiviral_cis_module(
        features=list(features),
        seq_len=seq_len,
        topology=topology,
        plasmid_id=pid
    )
    
    # Combine modules
    all_modules = non_pol3_modules + pol3_modules + lentiviral_cis_modules
    
    # Detect CDS modules within pol2 cassettes
    cds_modules = []
    pol2_types = {
        "pol2_expression_cassette",
        "nuclease_expression_cassette", 
        "reporter_expression_cassette",
        "yeast_pol2_expression_cassette",
        "pol2_expression_animal",
        "pol2_expression_generic",
        "pol2_expression_plant",
        "pol2_expression_fungal",
        "pol2_expression_insect",
    }
    
    for mod in all_modules:
        if mod.module_type in pol2_types:
            detected_cds = detect_cds_modules(mod, sequence, list(features), pid)
            cds_modules.extend(detected_cds)
    
    # Detect bacterial CDS modules (from RBS or promoter)
    bacterial_cds = detect_bacterial_cds_module(
        features=list(features),
        seq_len=len(sequence),
        topology="circular" if circular else "linear",
        plasmid_id=pid,
        seq=sequence,
    )
    cds_modules.extend(bacterial_cds)
    
    # Filter out pol2_promoter_only and pol3_promoter_only modules that are part of larger expression cassettes
    expression_cassette_types = {
        "pol2_expression_cassette", "nuclease_expression_cassette",
        "reporter_expression_cassette", "pol2_expression_animal",
        "pol2_expression_generic", "pol2_expression_plant",
        "pol2_expression_fungal", "pol2_expression_insect",
        "pol3_u6_sgrna_cassette", "pol3_expression_backbone",
        "bacterial_expression_cassette", "yeast_pol2_expression_cassette",
        "mammalian_lentiviral_expression",
    }
    
    # Find all expression cassettes
    cassettes = [m for m in all_modules if m.module_type in expression_cassette_types]
    
    # Filter promoter-only modules
    filtered_modules = []
    for mod in all_modules:
        if mod.module_type in ("pol2_promoter_only", "pol3_promoter_only"):
            # Check if this promoter is inside any expression cassette
            is_in_cassette = False
            for cassette in cassettes:
                if cassette.wraps or mod.wraps:
                    # Simplified check for wrapping
                    is_in_cassette = True
                    break
                elif mod.start >= cassette.start and mod.end <= cassette.end:
                    # Promoter is fully contained in cassette
                    is_in_cassette = True
                    break
            
            # Skip promoter-only if it's part of a cassette
            if is_in_cassette:
                continue
        
        filtered_modules.append(mod)
    
    all_modules = filtered_modules

    return all_modules, cds_modules


# ──────────────────────────────────────────────
# Updated public entry point
# ──────────────────────────────────────────────

def annotate_hierarchy_from_plannotate_v2(
    sequence: str,
    circular: bool,
    plannotate_rows: Sequence[Dict[str, Any]]
) -> Dict[str, Any]:
    """
    Enhanced hierarchical annotation with:
    - Pol III terminator detection via poly-T pattern
    - CDS module detection within expression cassettes
    - Support for dU6 and other Pol III promoter variants
    """
    features = infer_features_from_plannotate(plannotate_rows)
    analyzer = analyze_plasmid_from_plannotate(
        sequence=sequence,
        circular=circular,
        plannotate_rows=plannotate_rows
    )
    
    # Use enhanced extraction
    modules, cds_modules = extract_modules_with_cds_detection(
        sequence=sequence,
        circular=circular,
        features=features
    )
    
    # Combine all modules
    all_modules = modules + cds_modules
    
    # Select top level and generate annotations
    top_level = select_top_level_modules(modules, len(sequence))
    
    module_annotations: List[Dict[str, Any]] = []
    for mod in top_level:
        module_annotations.extend(_module_to_annotation(mod, len(sequence)))
    
    # Add child modules with multi-level nesting support
    # Level 1: Modules inside containers (e.g., lentiviral_payload inside mammalian_lentiviral_expression)
    # Level 2: Modules inside payloads (e.g., pol2_expression_animal inside lentiviral_payload)
    # Level 3: CDS modules inside expression cassettes
    
    child_modules = []
    
    # First, identify ALL modules inside containers (potential level1 children)
    potential_level1 = []
    for mod in modules:
        if mod in top_level or mod.module_type in ("cds_module", "cds_only"):
            continue
        
        for container in top_level:
            if container.module_type not in CONTAINER_MODULES:
                continue
            if mod.wraps or container.wraps:
                potential_level1.append((mod, container))
                break
            elif mod.start >= container.start and mod.end <= container.end:
                potential_level1.append((mod, container))
                break
    
    # Separate into level1 (direct children of containers) and level2 (children of level1)
    # A module is level2 if it's inside another potential_level1 module (especially lentiviral_payload)
    level1_children = []
    level2_children = []
    
    for mod, container in potential_level1:
        # Check if this module is inside any other potential_level1 module
        is_level2 = False
        for other_mod, other_container in potential_level1:
            if mod == other_mod:
                continue
            # Check if mod is inside other_mod
            if other_mod.module_type in ("lentiviral_payload", "shuttle_vector_backbone"):
                if mod.start >= other_mod.start and mod.end <= other_mod.end:
                    # This module is nested inside a payload/backbone, so it's level2
                    level2_children.append(mod)
                    is_level2 = True
                    break
        
        if not is_level2:
            level1_children.append(mod)
    
    
    # Now find modules nested inside level1_children (e.g., pol2/pol3 inside lentiviral_payload)
    level2_children = []
    for mod in modules:
        if mod in top_level or mod in level1_children or mod.module_type in ("cds_module", "cds_only"):
            continue
        
        # Check if this module is inside any level1 child (especially lentiviral_payload)
        for parent in level1_children:
            # For lentiviral_payload parents, pol2/pol3 modules must be fully contained
            if parent.module_type == "lentiviral_payload":
                # pol2/pol3 expression modules that start inside payload must end inside it
                if mod.module_type in ("pol2_expression_animal", "pol2_expression_generic", "pol2_expression_plant",
                                      "pol3_expression_backbone", "pol3_u6_sgrna_cassette"):
                    if mod.start >= parent.start and mod.end <= parent.end:
                        level2_children.append(mod)
                        break
            else:
                # For other parent types, use simple containment
                if mod.start >= parent.start and mod.end <= parent.end:
                    level2_children.append(mod)
                    break
    
    # Deduplicate each nesting level separately to preserve proper hierarchy
    # Level 1 deduplication
    deduplicated_level1 = []
    by_type_l1 = {}
    for mod in level1_children:
        by_type_l1.setdefault(mod.module_type, []).append(mod)
    
    for mod_type, mods in by_type_l1.items():
        if len(mods) == 1:
            deduplicated_level1.extend(mods)
            continue
        to_keep, skipped = [], set()
        for i, mod1 in enumerate(mods):
            if i in skipped:
                continue
            should_keep = True
            for j, mod2 in enumerate(mods):
                if i == j or j in skipped:
                    continue
                if _calculate_overlap_fraction(mod1, mod2, len(sequence)) > 0.9:
                    if mod2.length > mod1.length or (mod2.length == mod1.length and mod2.start < mod1.start):
                        should_keep = False
                        break
                    skipped.add(j)
            if should_keep:
                to_keep.append(mod1)
        deduplicated_level1.extend(to_keep)
    
    # Level 2 deduplication (separate from level 1)
    deduplicated_level2 = []
    by_type_l2 = {}
    for mod in level2_children:
        by_type_l2.setdefault(mod.module_type, []).append(mod)
    
    for mod_type, mods in by_type_l2.items():
        if len(mods) == 1:
            deduplicated_level2.extend(mods)
            continue
        to_keep, skipped = [], set()
        for i, mod1 in enumerate(mods):
            if i in skipped:
                continue
            should_keep = True
            for j, mod2 in enumerate(mods):
                if i == j or j in skipped:
                    continue
                if _calculate_overlap_fraction(mod1, mod2, len(sequence)) > 0.9:
                    if mod2.length > mod1.length or (mod2.length == mod1.length and mod2.start < mod1.start):
                        should_keep = False
                        break
                    skipped.add(j)
            if should_keep:
                to_keep.append(mod1)
        deduplicated_level2.extend(to_keep)
    
    # Cross-level deduplication: if same module type appears in both levels with >90% overlap,
    # prefer the one from the deeper level (level2)
    final_level1 = []
    for l1_mod in deduplicated_level1:
        should_keep = True
        for l2_mod in deduplicated_level2:
            if l1_mod.module_type == l2_mod.module_type:
                overlap_frac = _calculate_overlap_fraction(l1_mod, l2_mod, len(sequence))
                if overlap_frac > 0.9:
                    # Same type with >90% overlap - skip level1, keep level2
                    should_keep = False
                    break
        if should_keep:
            final_level1.append(l1_mod)
    
    # Add level1 annotations
    for mod in final_level1:
        module_annotations.extend(_module_to_annotation(mod, len(sequence)))
    
    # Cross-type preference: prefer pol3_u6_sgrna_cassette over pol3_expression_backbone
    # when they overlap >90%
    final_level2 = []
    for mod in deduplicated_level2:
        should_keep = True
        if mod.module_type == "pol3_expression_backbone":
            # Check if there's a pol3_u6_sgrna_cassette with >90% overlap
            for other in deduplicated_level2:
                if other.module_type == "pol3_u6_sgrna_cassette":
                    overlap_frac = _calculate_overlap_fraction(mod, other, len(sequence))
                    if overlap_frac > 0.9:
                        # Skip pol3_expression_backbone, keep pol3_u6_sgrna_cassette
                        should_keep = False
                        break
        if should_keep:
            final_level2.append(mod)
    
    # Add level2 annotations  
    for mod in final_level2:
        module_annotations.extend(_module_to_annotation(mod, len(sequence)))

    # Add CDS module annotations (these are children, not top-level)
    # First, remove any linker_module/nls_module annotations that fall inside CDS modules
    # These will be replaced by codon-aligned versions from CDS submodule resolution
    cds_ranges = [(cds_mod.start, cds_mod.end) for cds_mod in cds_modules]

    def _is_inside_cds(ann):
        ann_start = ann.get("start", 0)
        ann_end = ann.get("end", 0)
        for cds_start, cds_end in cds_ranges:
            if ann_start >= cds_start - 5 and ann_end <= cds_end + 5:
                return True
        return False

    # Filter out linker/nls modules inside CDS (they'll be re-added with correct boundaries)
    module_annotations = [
        ann for ann in module_annotations
        if not (ann.get("module_type") in ("linker_module", "nls_module") and _is_inside_cds(ann))
    ]

    # Deduplicate by position to avoid multiple annotations for same CDS
    seen_cds_positions = set()
    for cds_mod in cds_modules:
        pos_key = (cds_mod.start, cds_mod.end, cds_mod.payload_id)
        if pos_key not in seen_cds_positions:
            seen_cds_positions.add(pos_key)
            module_annotations.extend(_module_to_annotation(cds_mod, len(sequence)))
    
    # Resolve CDS sub-modules (protein modules, NLS, linkers, tags, gaps)
    from .Module_Library_gb.module_extractor import resolve_cds_submodules, Module, sha256_text, seq_hash
    plasmid_id = sha256_text(sequence)[:24]
    cds_submodules = []
    all_filtered_features = []  # Features removed during CDS resolution (>90% covered)
    all_boundary_corrections = []  # Boundary corrections for resolved features
    for cds_mod in cds_modules:
        try:
            result = resolve_cds_submodules(
                cds_mod, sequence, list(features), plasmid_id,
                Module, sha256_text, seq_hash
            )
            cds_submodules.extend(result["submodules"])
            all_filtered_features.extend(result.get("filtered_features", []))
            all_boundary_corrections.extend(result.get("boundary_corrections", []))
        except Exception as e:
            print(f"[WARN] CDS submodule resolution failed for {cds_mod.id}: {e}")
            import traceback
            traceback.print_exc()
    
    # Add CDS submodule annotations (deduplicated)
    seen_submod_positions = set()
    for submod in cds_submodules:
        pos_key = (submod.start, submod.end, submod.module_type)
        if pos_key not in seen_submod_positions:
            seen_submod_positions.add(pos_key)
            module_annotations.extend(_module_to_annotation(submod, len(sequence)))

    # Filter out mammalian_selection_cassette modules that overlap with CDS protein submodules
    # These selection markers are now resolved as protein_module within the CDS
    protein_submodule_ranges = [
        (ann.get('start', 0), ann.get('end', 0))
        for ann in module_annotations
        if ann.get('module_type') == 'protein_module'
        and ann.get('metadata', {}).get('module_family') == 'cds_submodule'
    ]

    def _overlaps_protein_submodule(ann):
        if ann.get('module_type') != 'mammalian_selection_cassette':
            return False
        ann_start = ann.get('start', 0)
        ann_end = ann.get('end', 0)
        ann_len = ann_end - ann_start
        if ann_len <= 0:
            return False
        for ps_start, ps_end in protein_submodule_ranges:
            overlap_start = max(ann_start, ps_start)
            overlap_end = min(ann_end, ps_end)
            overlap = max(0, overlap_end - overlap_start)
            if overlap / ann_len > 0.8:  # >80% overlap
                return True
        return False

    module_annotations = [
        ann for ann in module_annotations
        if not _overlaps_protein_submodule(ann)
    ]
    
    # Adjust feature boundaries to match submodules
    adjusted_features, features_to_delete = adjust_feature_boundaries_to_submodules(
        list(features), cds_modules, cds_submodules, sequence
    )
    
    # Build feature adjustments list for output
    feature_adjustments = []
    for orig in features:
        for adj in adjusted_features:
            if adj.name == orig.name and adj.canonical_id == orig.canonical_id:
                if adj.start != orig.start or adj.end != orig.end:
                    feature_adjustments.append({
                        "feature_name": orig.name,
                        "original_start": orig.start,
                        "original_end": orig.end,
                        "adjusted_start": adj.start,
                        "adjusted_end": adj.end,
                    })
                break
    
    # Build list of features to delete
    deleted_features = [
        {"feature_name": f.name, "start": f.start, "end": f.end, "reason": "no_submodule_pair"}
        for f in features_to_delete
    ]
    
    # Update analyzer's feature_instances with adjusted features
    # This ensures the GenBank output uses corrected feature boundaries
    feature_name_to_adjusted = {}
    for adj_feat in adjusted_features:
        key = (adj_feat.name, adj_feat.canonical_id or "")
        if key not in feature_name_to_adjusted:
            feature_name_to_adjusted[key] = []
        feature_name_to_adjusted[key].append(adj_feat)
    
    feature_name_to_delete = {(f.name, f.canonical_id or "") for f in features_to_delete}
    
    updated_feature_instances = []
    for fi in analyzer.get("feature_instances", []):
        fname = fi.get("feature_name", "")
        fcanonical_id = fi.get("canonical_id", "")
        fstart = fi.get("start", 0)
        fend = fi.get("end", 0)
        
        # Skip deleted features
        if (fname, fcanonical_id) in feature_name_to_delete:
            continue
        
        # Find matching adjusted feature
        key = (fname, fcanonical_id)
        if key in feature_name_to_adjusted:
            # Find the adjusted feature that was derived from this original feature
            for adj_feat in feature_name_to_adjusted[key]:
                # Use first match (should only be one per feature)
                fi = dict(fi)  # Make a copy
                fi["start"] = adj_feat.start
                fi["end"] = adj_feat.end
                updated_feature_instances.append(fi)
                break
        else:
            # Keep original if not adjusted
            updated_feature_instances.append(fi)
    
    analyzer["feature_instances"] = updated_feature_instances
    
    gap_annotations = gap_annotations_for_modules(top_level, len(sequence))

    # Detect motifs (start/stop codons, 2A peptides, Kozak sequences)
    motif_annotations: List[Dict[str, Any]] = []
    try:
        # Convert Feature objects to dicts with role for motif detection
        feature_dicts = []
        feature_instances = analyzer.get("feature_instances", [])
        for fi in feature_instances:
            role = fi.get("role", "") if isinstance(fi, dict) else getattr(fi, "role", "")
            ftype = fi.get("feature_type", "") if isinstance(fi, dict) else getattr(fi, "feature_type", "")
            start = fi.get("start", 0) if isinstance(fi, dict) else getattr(fi, "start", 0)
            end = fi.get("end", 0) if isinstance(fi, dict) else getattr(fi, "end", 0)
            strand = fi.get("strand", 1) if isinstance(fi, dict) else getattr(fi, "strand", 1)
            fname = fi.get("feature_name", "") if isinstance(fi, dict) else getattr(fi, "feature_name", "")
            feature_dicts.append({
                "role": role,
                "feature_type": ftype,
                "start": start,
                "end": end,
                "strand": strand,
                "feature_name": fname,
            })
        motif_detector = MotifDetector(sequence, circular)
        detected_motifs = motif_detector.detect_all_motifs(feature_dicts)
        for motif in detected_motifs:
            motif_annotations.append({
                "name": motif.name,
                "start": motif.start,
                "end": motif.end,
                "direction": motif.strand,
                "color": _motif_color(motif.motif_type),
                "source": "motif_detector",
                "layer": "motif",
                "type": "motif",
                "motif_type": motif.motif_type,
                "description": motif.description,
                "sequence": motif.sequence,
            })
    except Exception as e:
        print(f"[WARN] Motif detection failed: {e}")

    # Debug: print motifs before filtering
    print(f"[DEBUG] Before filtering: {len(motif_annotations)} motifs")
    feature_instances = analyzer.get("feature_instances", [])
    # Debug: Show ALL feature roles to understand why CDS features might be missing
    print(f"[DEBUG] Total features: {len(feature_instances)}")
    all_roles = [(f.get('feature_name', '')[:30], f.get('role', ''), f.get('feature_type', '')) for f in feature_instances]
    print(f"[DEBUG] Feature roles (first 20): {all_roles[:20]}")
    cds_features_debug = [
        f for f in feature_instances
        if f.get("feature_type", "").lower() == "cds" or
           f.get("role", "") in ("expression_payload", "editing_payload", "reporter_payload", "selection_payload", "bacterial_marker")
    ]
    print(f"[DEBUG] CDS features for filtering: {[(f.get('feature_name'), f.get('start'), f.get('end'), f.get('role')) for f in cds_features_debug]}")
    
    # Filter motifs to only those at actual CDS boundaries
    motif_annotations = _filter_motifs_by_modules(motif_annotations, analyzer.get("module_instances", []), feature_instances, sequence)
    print(f"[DEBUG] After filtering: {len(motif_annotations)} motifs")

    return {
        "features_inferred": len(features),
        "motif_annotations": motif_annotations,
        "module_count": len(top_level) + len(cds_modules),
        "modules": [
            {
                "id": mod.id,
                "module_type": mod.module_type,
                "start": mod.start,
                "end": mod.end,
                "wraps": mod.wraps,
                "payload_id": mod.payload_id,
                "metadata": mod.metadata,
            }
            for mod in all_modules
        ],
        "cds_modules": [
            {
                "id": mod.id,
                "module_type": mod.module_type,
                "start": mod.start,
                "end": mod.end,
                "payload_id": mod.payload_id,
                "metadata": mod.metadata,
            }
            for mod in cds_modules
        ],
        "module_annotations": module_annotations,
        "feature_adjustments": feature_adjustments,
        "deleted_features": deleted_features,
        "gap_annotations": gap_annotations,
        "feature_instances": analyzer.get("feature_instances", []),
        "module_instances": analyzer.get("module_instances", []),
        "junctions": analyzer.get("junctions", []),
        "build_profiles": analyzer.get("build_profiles", []),
        "construct_graph": analyzer.get("construct_graph", {}),
        "capabilities": analyzer.get("capabilities", []),
        "intent_inference": analyzer.get("intent_inference", {}),
        "rule_findings": analyzer.get("rule_findings", []),
        "analyzer_summary": analyzer.get("summary", {}),
        "cds_filtered_features": all_filtered_features,
        "cds_boundary_corrections": all_boundary_corrections,
    }
