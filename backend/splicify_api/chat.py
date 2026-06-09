"""
Chat orchestration endpoint for AI Plasmid Design.
Mirrors the n8n workflow:
  1. extract sequences  →  fragments_in (complete sequences, multi-line joined)
  2. redact sequences   →  compact prompt for LLM
  3. parse intent       →  intent + gibson/pcr parameters (n8n system prompt)
  4. merge              →  fragments_in + routing/params
  5. build payload      →  proper request object for each workflow
  6. call endpoint      →  design_gibson_primers / design_primers / etc.
  7. explanation        →  optional AI summary
  8. return             →  {ok, reply, viz, viz_list, files, sessionId}
"""
from __future__ import annotations

import logging
import os
import base64
import uuid
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, File, Form, UploadFile

logger = logging.getLogger(__name__)

# Unified pre-design system imports
from .predesign import (
    DesignRequest,
    PartSpecification,
    TargetSpecification,
    InputSource,
    PartResolver,
    ResolvedPart,
    TargetPlasmidBuilder,
    CloningRouter,
    WorkflowMethod,
    get_knowledge_base,
)

from .explanation import generate_explanation
from .extractors import build_fragment_objects, extract_sequences, redact_sequences
from .gibson_primers import GibsonFragment, GibsonRequest, design_gibson_primers
from .intent import parse_intent
from .request_normalizer import normalize_request
from .feature_dna_resolver import get_feature_dna
from .plannotate_router import AnnotateRequest, annotate_genbank, _load_kb_records
from .restriction_cloning_designer import design_restriction_cloning
from .files import (
    build_gibson_files,
    build_sdm_files,
)
from .visualization import (
    build_gibson_viz,
    build_sdm_viz,
)
from .cloning.sdm_parser import SDMMutationParser
from .cloning.sdm_operator import SDMOperator
from .cloning.sgrna_oligo_designer import (
    assemble_sgrna_plasmid,
    design_sgrna_oligos,
    design_sgrna_oligos_from_vector,
    load_lenticrispr_v2,
    parse_genbank_features,
    parse_genbank_sequence,
)
from .cloning.golden_gate_primer_designer import (
    design_multi_fragment_assembly,
    design_single_fragment_replacement,
    design_scarless_deletion,
    build_design_response,
)
from .files import build_sdm_files, build_restriction_files
import re
from pathlib import Path
from . import _data

router = APIRouter()


# -------------------------
# Feature Lookup from Knowledge Base
# -------------------------

def lookup_genetic_features(feature_names: List[str], feature_type: str = 'CDS') -> List[Dict[str, Any]]:
    """
    Look up genetic features from the Module_Library_gb knowledge base.
    Now extracts specific features (CDS, promoter, etc.) from GenBank files.

    Args:
        feature_names: List of feature names to look up (e.g., ["GFP", "mCherry"])
        feature_type: Type of feature to extract ('CDS', 'promoter', 'terminator', 'misc_feature')

    Returns:
        List of dicts with {"name": str, "sequence": str, "type": str, "length": int, ...}
    """
    from Bio import SeqIO

    module_lib_path = _data.data_path("Module_Library_gb")
    if not module_lib_path.exists():
        logger.warning(f"Module library not found at {module_lib_path}")
        return []

    results = []

    # Common feature name mappings
    feature_mappings = {
        "ef1a": ["EF1a", "EF-1a", "elongation factor 1 alpha", "EF1alpha"],
        "mcherry": ["mCherry", "mcherry", "cherry"],
        "egfp": ["eGFP", "EGFP", "enhanced GFP", "TurboGFP", "Cycle 3 GFP"],
        "gfp": ["GFP", "green fluorescent protein", "TurboGFP", "eGFP", "Cycle 3 GFP"],
        "bgh": ["bGH", "BGH", "bovine growth hormone"],
        "polya": ["poly", "polyA", "polyadenylation"],
        "cmv": ["CMV", "cytomegalovirus"],
        "sv40": ["SV40", "simian virus 40"],
        "puro": ["puro", "puromycin"],
        "neo": ["neo", "neomycin"],
    }

    for feature_name in feature_names:
        feature_lower = feature_name.lower()

        # Determine search terms
        search_terms = [feature_name]
        for key, aliases in feature_mappings.items():
            if key in feature_lower:
                search_terms.extend(aliases)

        # Search in GenBank files
        found = False
        for gb_file in module_lib_path.rglob("*.gb"):
            if found:
                break
            if "Module_Library_JSON" in str(gb_file):
                continue

            try:
                record = SeqIO.read(str(gb_file), "genbank")

                # Extract features of requested type
                for feat in record.features:
                    if feat.type != feature_type:
                        continue

                    # Get feature name/label
                    feat_name = feat.qualifiers.get("label",
                                feat.qualifiers.get("gene",
                                feat.qualifiers.get("product", [""])))
                    if isinstance(feat_name, list):
                        feat_name = feat_name[0] if feat_name else ""

                    # Check if this matches our search
                    if any(term.lower() in feat_name.lower() for term in search_terms):
                        # Extract feature sequence
                        feat_seq = str(feat.extract(record.seq))

                        results.append({
                            "name": feat_name,
                            "sequence": feat_seq,
                            "type": feat.type,
                            "length": len(feat_seq),
                            "source": str(gb_file.relative_to(module_lib_path)),
                            "source_plasmid": record.name,
                            "location": f"{feat.location.start}-{feat.location.end}"
                        })

                        logger.info(f"Found {feat_name} ({feat.type}): {len(feat_seq)} bp from {gb_file.name}")
                        found = True
                        break

            except Exception as e:
                logger.warning(f"Error reading {gb_file}: {e}")

    return results


# ---------------------------------------------------------------------------
# Part-reference extraction and identification against the pLannotate feature KB.
# ---------------------------------------------------------------------------

# Tokens stripped from the *start* of each candidate fragment during extraction.
# Covers common prompt preambles like "design golden gate primers to assemble X + Y + Z".
_PART_PREAMBLE_STOPS = {
    "design", "create", "make", "build", "generate", "produce", "construct",
    "golden", "gate", "gibson", "assembly",
    "primers", "primer", "oligos", "oligo",
    "to", "assemble", "assembling", "for", "using", "with", "from", "of", "and",
    "the", "a", "an", "my", "our", "this", "these", "that",
    "please", "could", "would", "you", "can", "i", "we", "want",
}

# Role suffixes mapped to pLannotate feature_type. Longest suffixes first.
_ROLE_SUFFIXES = [
    ("polyadenylation signal", "polyA_signal"),
    ("poly(a) signal", "polyA_signal"),
    ("poly-a signal", "polyA_signal"),
    ("poly a signal", "polyA_signal"),
    ("polya signal", "polyA_signal"),
    ("poly(a)", "polyA_signal"),
    ("poly a", "polyA_signal"),
    ("polya", "polyA_signal"),
    ("terminator", "terminator"),
    ("promoter", "promoter"),
    ("enhancer", "enhancer"),
    ("intron", "intron"),
    ("origin", "rep_origin"),
    ("cds", "CDS"),
    ("gene", "CDS"),
    ("orf", "CDS"),
    ("tag", "CDS"),
]



# ---------------------------------------------------------------------------
# Canonical-request adapters for per-workflow handlers.
# ---------------------------------------------------------------------------

def _restriction_inputs_from_canonical(
    canonical_request,
    intent_result: Dict[str, Any],
    seq_data: Dict[str, Any],
) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[Tuple[str, str]]]:
    """
    Extract (insert_name, insert_sequence, vector_name, enzyme_override)
    for the restriction_cloning handler from the canonical request, with a
    defensive fallback to intent_result when canonical_request is None.

    Type-II enzymes are pulled preferentially over Type-IIs; if the user named
    exactly two Type-II cutters the tuple is returned for the designer's
    `enzyme_override` arg. >2 named enzymes → None so the designer auto-picks
    from the Type-II preferred list rather than silently truncating.
    """
    insert_name = None
    insert_sequence = None
    vector_name = None
    enzyme_override = None

    if canonical_request is not None:
        if canonical_request.vector and canonical_request.vector.name:
            vector_name = canonical_request.vector.name

        insert_parts = canonical_request.parts_by_role("insert")
        if not insert_parts:
            # Some prompts put the insert in a labeled fragment ("Frag1: ATGC…")
            # without a role keyword. Fall through to the first prompt-sequence
            # part we have.
            insert_parts = [p for p in canonical_request.parts
                             if p.source == "prompt_sequence"]
        if insert_parts:
            chosen = insert_parts[0]
            insert_name = chosen.name
            insert_sequence = chosen.sequence

        type_ii = [e for e in canonical_request.enzymes if e.kind == "type_ii"
                   and e.preference != "avoid"]
        if len(type_ii) == 2:
            enzyme_override = (type_ii[0].name, type_ii[1].name)

    # Defensive fallback — intent_result path mirrors the pre-migration logic.
    if insert_name is None and insert_sequence is None:
        rc_params = (intent_result or {}).get("restriction_cloning") or {}
        fallback_name = (rc_params.get("insert_name") or "").strip() or None
        fallback_seq = (rc_params.get("insert_sequence") or "").strip() or None
        if fallback_name or fallback_seq:
            insert_name = fallback_name
            insert_sequence = fallback_seq
    if vector_name is None:
        rc_params = (intent_result or {}).get("restriction_cloning") or {}
        vector_name = (rc_params.get("vector_name") or "").strip() or None
    if enzyme_override is None:
        rc_params = (intent_result or {}).get("restriction_cloning") or {}
        enz_list = rc_params.get("enzymes") or []
        if isinstance(enz_list, list) and len(enz_list) == 2:
            enzyme_override = (enz_list[0], enz_list[1])

    # Prompt-pasted literal sequences (seq_data) win when no canonical/intent
    # sequence is available — matches the prior handler's fall-through.
    if insert_sequence is None and seq_data and seq_data.get("count", 0) >= 1:
        insert_sequence = seq_data["sequences"][0].upper()

    return insert_name, insert_sequence, vector_name, enzyme_override


def _normalize_feature_text(s: str) -> str:
    """Normalize a feature name for comparison: lowercase, unicode → ascii, collapse punctuation to spaces."""
    if not s:
        return ""
    out = s.lower().replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    out = re.sub(r"poly\(a\)", "polya", out)
    out = re.sub(r"[\s\-_()/.,]+", " ", out)
    return re.sub(r"\s+", " ", out).strip()


def _features_from_annotate_response(ann_result) -> List[Dict[str, Any]]:
    """Convert an AnnotateResponse into the flat feature-dict shape the SDM
    parser expects. Pulls from both `annotations` (raw pLannotate hits, with
    "start..end" location strings) and `modules` (hierarchical module dicts)."""
    out: List[Dict[str, Any]] = []
    # 1. Raw annotations: location is "start..end"
    for ann in (getattr(ann_result, "annotations", None) or []):
        loc = getattr(ann, "location", "") or ""
        try:
            parts = loc.replace("..", ",").split(",")
            start = int(parts[0]) if parts and parts[0] else 0
            end = int(parts[1]) if len(parts) > 1 and parts[1] else start
        except (ValueError, IndexError):
            start, end = 0, 0
        name = (getattr(ann, "name", None)
                or (ann.qualifiers or {}).get("label")
                or (ann.qualifiers or {}).get("gene")
                or "")
        out.append({
            "type": getattr(ann, "type", "misc_feature") or "misc_feature",
            "start": start,
            "end": end,
            "strand": 1,
            "name": name,
        })
    # 2. Hierarchical modules: dict shape, may overlap with annotations
    for mod in (getattr(ann_result, "modules", None) or []):
        start = mod.get("start")
        end = mod.get("end")
        try:
            start_int = int(start) if start is not None else 0
            end_int = int(end) if end is not None else 0
        except (TypeError, ValueError):
            continue
        out.append({
            "type": mod.get("type", mod.get("module_type", "misc_feature")),
            "start": start_int,
            "end": end_int,
            "strand": 1,
            "name": (mod.get("name") or mod.get("label")
                     or mod.get("payload_id") or mod.get("module_type") or ""),
        })
    return out


def extract_part_candidates(message: str) -> List[Dict[str, Any]]:
    """Parse a free-text message into candidate part references.

    Handles list styles:
      - "X + Y + Z"
      - "X, Y, Z" / "X, Y, and Z"
      - "an X, a Y, and a Z"
      - "X-Y-Z" (compact hyphen-joined fallback)

    For each extracted chunk, strips preamble stopwords, then detects a trailing
    role keyword (promoter, polyA, terminator, etc.) and maps it to a pLannotate
    feature_type. Returns a list of {"name", "feature_type"} dicts.
    """
    msg = (message or "").strip().rstrip(".!?")
    if not msg:
        return []

    parts = re.split(r"\s*\+\s*|\s*,\s*|\s+and\s+|\s+then\s+", msg)
    parts = [p.strip() for p in parts if p and p.strip()]

    if len(parts) < 2:
        # Fallback: compact hyphen-joined form like "CMV-eGFP-bGHpolyA"
        for tok in msg.split():
            if tok.count("-") >= 2:
                parts = [p for p in tok.split("-") if p]
                break

    out: List[Dict[str, Any]] = []
    for raw in parts:
        words = re.split(r"\s+", raw)
        while words and words[0].lower().strip(".,;:") in _PART_PREAMBLE_STOPS:
            words.pop(0)
        if not words:
            continue
        name = " ".join(words).strip(".,;:!? ")
        if not name or len(name) < 2:
            continue

        feature_type: Optional[str] = None
        name_lower = name.lower()
        for suffix, ftype in _ROLE_SUFFIXES:
            if name_lower == suffix:
                name = ""
                break
            if name_lower.endswith(" " + suffix):
                feature_type = ftype
                name = name[: -(len(suffix) + 1)].strip()
                break

        if not name:
            continue
        out.append({"name": name, "feature_type": feature_type})

    return out


def identify_features_from_kb(candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Identify candidate part references in the pLannotate feature KB.

    Matching strategy per candidate (each filtered by feature_type when provided):
      1. Exact match (normalized) on feature_name / normalized_feature_name / sseqid / alternative_names.
      2. Substring match on the same fields; prefer shorter feature_name (more atomic).

    Returns a list of identified features with fields suitable for downstream
    workflows: name, feature_id, feature_type, sequence, length, source,
    descriptions, matched_via, query, requested_type.
    """
    records = _load_kb_records()
    if not records:
        logger.warning("pLannotate feature KB is empty or unavailable")
        return []

    results: List[Dict[str, Any]] = []
    for cand in candidates:
        raw_name = (cand.get("name") or "").strip()
        if not raw_name:
            continue
        qn = _normalize_feature_text(raw_name)
        if not qn:
            continue
        requested_type = cand.get("feature_type")

        # Primary pool: records of the requested feature_type. If the user used
        # imprecise terminology (e.g. "terminator" for a polyA signal) we
        # retry against a broader pool below.
        pool = [r for r in records if (not requested_type) or r.get("feature_type") == requested_type]
        if not pool:
            pool = records

        # Known type-confusion equivalences: a user saying "terminator" often
        # means a polyA signal in mammalian expression contexts, and vice versa.
        type_fallbacks = {
            "terminator": ["polyA_signal"],
            "polyA_signal": ["terminator"],
            "CDS": ["gene", "protein_generator"],
        }

        exacts: List[Dict[str, Any]] = []
        partials: List[tuple] = []  # (priority_len, record)
        for r in pool:
            names = [
                _normalize_feature_text(r.get("feature_name") or ""),
                _normalize_feature_text(r.get("normalized_feature_name") or ""),
                _normalize_feature_text(r.get("sseqid") or ""),
            ]
            names += [_normalize_feature_text(a) for a in (r.get("alternative_names") or [])]
            names = [n for n in names if n]
            if qn in names:
                exacts.append(r)
                continue
            qn_compact = qn.replace(" ", "")
            for n in names:
                n_compact = n.replace(" ", "")
                if qn in n.split() or qn_compact in n_compact:
                    partials.append((len(r.get("feature_name") or ""), r))
                    break

        # Fallback: if nothing matched the requested type, expand the pool to
        # include type equivalences (e.g. terminator <-> polyA_signal) and
        # rescan. This catches prompts where the user's role label maps to a
        # different KB feature_type.
        if requested_type and not exacts and not partials and requested_type in type_fallbacks:
            alt_types = type_fallbacks[requested_type]
            alt_pool = [r for r in records if r.get("feature_type") in alt_types]
            for r in alt_pool:
                names = [
                    _normalize_feature_text(r.get("feature_name") or ""),
                    _normalize_feature_text(r.get("normalized_feature_name") or ""),
                    _normalize_feature_text(r.get("sseqid") or ""),
                ]
                names += [_normalize_feature_text(a) for a in (r.get("alternative_names") or [])]
                names = [n for n in names if n]
                if qn in names:
                    exacts.append(r)
                    continue
                qn_compact = qn.replace(" ", "")
                for n in names:
                    n_compact = n.replace(" ", "")
                    if qn in n.split() or qn_compact in n_compact:
                        partials.append((len(r.get("feature_name") or ""), r))
                        break

        chosen = None
        if exacts:
            # Prefer records whose representative_sequence is non-empty.
            for r in exacts:
                if ((r.get("intrinsic_properties") or {}).get("sequence_derived", {}).get("representative_sequence")):
                    chosen = r
                    break
            if chosen is None:
                chosen = exacts[0]
        elif partials:
            partials.sort(key=lambda x: x[0])
            for _, r in partials:
                if ((r.get("intrinsic_properties") or {}).get("sequence_derived", {}).get("representative_sequence")):
                    chosen = r
                    break
            if chosen is None:
                chosen = partials[0][1]

        if chosen is None:
            continue

        raw_seq = (
            (chosen.get("intrinsic_properties") or {})
            .get("sequence_derived", {})
            .get("representative_sequence")
            or ""
        )
        # Ask the unified resolver for DNA — it preserves a direct-DNA hit,
        # reaches into feature_reference.fna, and back-translates from
        # feature_protein_kb.json / feature_protein.faa when neither of those
        # carries DNA. Without this step, post-GenoLIB-rebuild CDS entries
        # (eGFP, Cas9, selection markers) hand protein strings to cloning
        # workflows, which silently produce garbage primers.
        sseqid = chosen.get("sseqid") or ""
        resolved = get_feature_dna(
            sseqid,
            protein_sequence=(raw_seq if raw_seq and not all(c in "ACGTNacgtn" for c in raw_seq[:50]) else None),
            json_representative_sequence=raw_seq if raw_seq else None,
        )
        seq = (resolved.sequence or "").upper()
        descriptions = (chosen.get("source") or {}).get("descriptions", []) or []

        results.append({
            "name": chosen.get("feature_name"),
            "feature_id": chosen.get("feature_id"),
            "feature_type": chosen.get("feature_type"),
            "sequence": seq,
            "length": len(seq),
            "source": "plannotate_feature_kb",
            "descriptions": descriptions,
            "matched_via": "exact" if exacts else "substring",
            "query": raw_name,
            "requested_type": requested_type,
            "sequence_provenance": resolved.provenance,
            "sequence_provenance_notes": resolved.notes,
            "sequence_organism": resolved.organism,
        })

    return results


_HELP_TEXT = """I can help you with molecular biology design tasks:

**1. Gibson Assembly** — provide two or more DNA fragment sequences
> *"Design gibson assembly primers for: Frag1: ATGC..., Frag2: GCTA..."*

**2. PCR Primer Design** — provide a single template sequence
> *"Design PCR primers for this template: ATGCATGC..."*

**3. Batch PCR** — provide multiple template sequences
> *"Design primers for these templates: Template1: ATGC..., Template2: GCTA..."*

**4. Inventory-Based Gibson** — upload a target plasmid + inventory plasmid files

**5. Plasmid Annotation** — upload a GenBank file
> *"Annotate this plasmid"*

**6. Site-Directed Mutagenesis** — upload a GenBank file and describe the mutation
> *"Delete the His-tag"*
> *"Insert FLAG tag after position 500"*
> *"Change codon 45 from Arg to Ala"*
> *"Substitute AATTCC with GGCCAA"*

**7. sgRNA Golden Gate Cloning** — design oligos to clone a gRNA into lentiCRISPR v2 or similar vectors
> *"Design oligos to clone gRNA GAGTCCGAGCAGAAGAAGAA into lentiCRISPR v2"*
> *"Clone this guide RNA ATGCATGCATGCATGCATGC using Golden Gate"*

Please provide your sequences or upload files to get started!"""


# ---------------------------------------------------------------------------
# Unified Pre-Design Helper Functions
# ---------------------------------------------------------------------------

async def _build_design_request_from_chat(
    message: str,
    seq_data: Dict[str, Any],
    intent_result: Dict[str, Any],
    file: Optional[UploadFile],
    inventory_files: Optional[List[UploadFile]],
    session_id: str
) -> DesignRequest:
    """
    Build a DesignRequest from chat input.

    This converts the parsed chat message, extracted sequences, and files
    into a standardized DesignRequest for the pre-design pipeline.
    """
    parts = []

    # Parse parts from extracted sequences
    for i, seq in enumerate(seq_data.get("sequences", [])):
        label = seq_data.get("labels", [])[i] if i < len(seq_data.get("labels", [])) else f"Fragment{i+1}"
        parts.append(
            PartSpecification(
                name=label,
                source=InputSource.DIRECT_SEQUENCE,
                sequence=seq,
                specified_order=i + 1,
            )
        )

    # Check for feature names in message (common patterns)
    # Examples: "CMV promoter", "eGFP", "mCherry", "bGH polyA"
    feature_patterns = [
        r'\b(CMV|EF1a|SV40|T7|lac|ara)\s*(?:promoter)?\b',
        r'\b(eGFP|mCherry|mRuby|EGFP|GFP|RFP|YFP)\b',
        r'\b(bGH|SV40|polyA|poly\s*A)\s*(?:polyA|terminator|signal)?\b',
        r'\b(AmpR|KanR|NeoR|resistance)\b',
    ]

    feature_matches = []
    for pattern in feature_patterns:
        for match in re.finditer(pattern, message, re.IGNORECASE):
            feature_name = match.group(0)
            if feature_name not in feature_matches:
                feature_matches.append(feature_name)

    # Add feature name parts (if not already added as sequences)
    for feature_name in feature_matches:
        parts.append(
            PartSpecification(
                name=feature_name,
                source=InputSource.FEATURE_NAME,
                feature_name=feature_name.strip(),
            )
        )

    # Build target specification
    target = None
    if file:
        target = TargetSpecification(
            source="uploaded",
            uploaded_file_id=f"file_{session_id}",
            topology="circular",
        )
    elif parts:
        target = TargetSpecification(
            source="assembled",
            parts=parts,
            assembly_order="listed",
            topology="circular",
        )

    # Build inventory file IDs
    inventory_ids = []
    if inventory_files:
        for i, inv_file in enumerate(inventory_files):
            inventory_ids.append(f"inv_{session_id}_{i}")

    return DesignRequest(
        session_id=session_id,
        user_message=message,
        parts=parts,
        target=target,
        inventory_file_ids=inventory_ids,
        suggested_workflow=intent_result.get("intent"),
        metadata={
            "intent_result": intent_result,
            "seq_count": seq_data.get("count", 0),
        },
    )


async def _build_file_cache(
    file: Optional[UploadFile],
    inventory_files: Optional[List[UploadFile]],
    session_id: str
) -> Dict[str, Dict[str, str]]:
    """Build file cache for part resolution."""
    cache = {}

    if file:
        content = await file.read()
        try:
            text = content.decode("utf-8")
        except UnicodeDecodeError:
            text = content.decode("latin-1", errors="replace")

        cache[f"file_{session_id}"] = {
            "name": file.filename or "uploaded.gb",
            "content": text,
        }
        await file.seek(0)

    if inventory_files:
        for i, inv_file in enumerate(inventory_files):
            content = await inv_file.read()
            try:
                text = content.decode("utf-8")
            except UnicodeDecodeError:
                text = content.decode("latin-1", errors="replace")

            cache[f"inv_{session_id}_{i}"] = {
                "name": inv_file.filename or f"inventory_{i}.gb",
                "content": text,
            }
            await inv_file.seek(0)

    return cache


async def _execute_inventory_workflow(
    design_request: DesignRequest,
    file_cache: Dict[str, Dict[str, str]],
    message: str,
    session_id: str,
    include_explanation: bool,
    intent_result: Dict[str, Any],
) -> Optional[Dict[str, Any]]:
    """
    Execute inventory-based workflow (target + inventory → homology extraction).

    For now, this falls back to legacy inv_gib handling which has been enhanced
    with pre-design context. In the future, this will use HomologyExtractor to
    support multiple workflows (Gibson, Golden Gate, Restriction, Gateway).

    Returns:
        None to fall back to legacy inv_gib path (which now has proper context)
    """
    logger.info("[Pre-Design] Falling back to legacy inv_gib with enhanced context")
    return None


async def _execute_unified_predesign(
    message: str,
    seq_data: Dict[str, Any],
    intent_result: Dict[str, Any],
    file: Optional[UploadFile],
    inventory_files: Optional[List[UploadFile]],
    session_id: str,
    include_explanation: bool
) -> Optional[Dict[str, Any]]:
    """
    Execute the unified pre-design pipeline.

    Returns:
        Result dict with workflow output, or None to fall back to legacy code
    """
    try:
        # Build design request
        design_request = await _build_design_request_from_chat(
            message, seq_data, intent_result, file, inventory_files, session_id
        )

        logger.info(f"[Pre-Design] Created DesignRequest with {len(design_request.parts)} parts")

        if not design_request.parts and not design_request.target:
            return None

        # Build file cache
        file_cache = await _build_file_cache(file, inventory_files, session_id)
        context = {"file_cache": file_cache}

        # Check for inventory-based workflow (target + inventory files)
        has_target_file = design_request.target and design_request.target.source == "uploaded"
        has_inventory = len(design_request.inventory_file_ids) > 0

        if has_target_file and has_inventory:
            # Inventory-based workflow (homology extraction)
            logger.info("[Pre-Design] Detected inventory-based workflow (target + inventory)")
            return await _execute_inventory_workflow(
                design_request,
                file_cache,
                message,
                session_id,
                include_explanation,
                intent_result,
            )

        # Regular assembly workflow (direct parts)
        # Resolve parts
        resolver = PartResolver()

        try:
            kb = get_knowledge_base()
            if not kb._loaded:
                kb.load()
            resolver.knowledge_base = kb
            logger.info(f"[Pre-Design] KB loaded: {kb.get_total_sequences()} sequences")
        except Exception as e:
            logger.warning(f"[Pre-Design] KB not available: {e}")

        resolved_parts = await resolver.resolve_all(design_request.parts, context)
        logger.info(f"[Pre-Design] Resolved {len(resolved_parts)} parts")

        # Annotate every resolved part once (full pipeline, cached by sequence
        # hash). Downstream consumers — the workflow assessor, the response
        # builders, and the (later) LLM orchestrator — read modules + roles
        # off the annotation rather than re-running it.
        from .annotation_cache import annotate_cached
        part_annotations: List[Dict[str, Any]] = []
        for part in resolved_parts:
            try:
                ann = await annotate_cached(part.sequence, circular=False, depth="full")
                part_annotations.append({
                    "name": part.name,
                    "length": part.length,
                    "modules": ann.get("module_annotations") or ann.get("modules") or [],
                    "interactions": ann.get("interactions") or [],
                })
            except Exception as exc:
                logger.warning("[Pre-Design] part annotation failed for %s: %s", part.name, exc)
                part_annotations.append({"name": part.name, "length": part.length,
                                         "modules": [], "interactions": []})

        # Build target plasmid
        builder = TargetPlasmidBuilder()
        topology = "linear" if "linear" in message.lower() else "circular"
        target = builder.build_from_parts(resolved_parts, topology=topology)
        logger.info(f"[Pre-Design] Target: {target.length} bp, {topology}")

        # Annotate the assembled target at module+interaction depth only
        # (cloning-feature pass is skipped here — re-run at full depth on the
        # final response target if the workflow ships one).
        target_modules: List[Dict[str, Any]] = []
        target_interactions: List[Dict[str, Any]] = []
        try:
            # 2026-05-12: use the _llm-endpoint annotator so the verifier
            # sees full rule-based modules + interactions (the hierarchy
            # endpoint used by annotate_cached drops both, which blinded
            # the orientation-mismatch detector + auto-correct).
            from .annotation_cache import annotate_llm_cached as _alc_predesign
            target_ann = await _alc_predesign(
                target.sequence,
                circular=(target.topology == "circular"),
            )
            # Rule-based modules live under hierarchical_annotations with a
            # module_type tag in the _llm response shape.
            target_modules = [h for h in (target_ann.get("hierarchical_annotations") or [])
                              if h.get("module_type")]
            target_interactions = target_ann.get("interactions") or []
        except Exception as exc:
            logger.warning("[Pre-Design] target annotation failed: %s", exc)

        # ----- Design verification on the assembled target -------------
        # Runs after target annotation; results ride on predesign_context as
        # `target_verification` so per-intent handlers + the workflow trace
        # can render orientation / assembly / required-change suggestions.
        from .target_from_inventory_router import verify_target_design as _verify_target
        target_verification = None
        try:
            target_annotations_for_verify = (
                target_ann.get("annotations") if isinstance(target_ann, dict) else []
            ) or []
            target_verification = _verify_target(
                target_sequence=target.sequence,
                target_annotations=target_annotations_for_verify,
                target_modules=target_modules,
                target_interactions=target_interactions,
                target_name=getattr(target, "name", "target"),
            )
            logger.info(
                "[Pre-Design] target verification: passed=%s, warnings=%d (by kind: %s)",
                target_verification["passed"],
                target_verification["summary"]["total_warnings"],
                target_verification["summary"]["by_kind"],
            )
        except Exception as exc:
            logger.warning("[Pre-Design] target verification failed: %s", exc)

        # ----- v14: auto-correct KB-referenced parts on verification failure
        # Only act on KB-resolved features the user explicitly requested.
        target_corrections: List[Dict[str, Any]] = []
        try:
            from .target_from_inventory_router import auto_correct_kb_part_orientation
            kb_part_names: List[str] = []
            for p in resolved_parts:
                if getattr(p, "origin", None) == "knowledge_base":
                    if p.name:
                        kb_part_names.append(p.name)
                    if getattr(p, "canonical_id", None):
                        kb_part_names.append(p.canonical_id)
            for ident in (intent_result or {}).get("kb_resolved", {}).get("identified", []) or []:
                nm = ident.get("feature_name") or ident.get("name")
                if nm:
                    kb_part_names.append(nm)

            if kb_part_names and target_verification and not target_verification.get("passed"):
                async def _reannotate(seq):
                    # Use _llm annotator so auto-correct sees the same
                    # interaction graph + rule-based modules that the
                    # verifier uses.
                    from .annotation_cache import annotate_llm_cached as _alc_ac
                    ann_payload = await _alc_ac(
                        seq, circular=(target.topology == "circular"),
                    )
                    return (
                        ann_payload.get("annotations") or [],
                        [h for h in (ann_payload.get("hierarchical_annotations") or [])
                         if h.get("module_type")],
                        ann_payload.get("interactions") or [],
                    )

                target_ann_for_verify = (
                    target_ann.get("annotations") if isinstance(target_ann, dict) else []
                ) or []
                ac = await auto_correct_kb_part_orientation(
                    target_sequence=target.sequence,
                    target_annotations=target_ann_for_verify,
                    target_modules=target_modules,
                    target_interactions=target_interactions,
                    kb_resolved_part_names=kb_part_names,
                    reannotate_async=_reannotate,
                )
                target_corrections = ac.get("corrections") or []
                if ac.get("corrected_sequence"):
                    logger.info(
                        "[Pre-Design] applied %d auto-correction(s) for KB parts: %s",
                        len([c for c in target_corrections if c.get("result") == "applied"]),
                        ", ".join(c["feature"] for c in target_corrections
                                  if c.get("result") == "applied"),
                    )
                    # Adopt the corrected target — refresh sequence + length +
                    # re-annotation outputs the predesign_context expects.
                    target.sequence = ac["corrected_sequence"]
                    target.length = len(ac["corrected_sequence"])
                    from .annotation_cache import annotate_llm_cached as _alc_post
                    new_ann_payload = await _alc_post(
                        target.sequence, circular=(target.topology == "circular"),
                    )
                    target_modules = [h for h in (new_ann_payload.get("hierarchical_annotations") or [])
                                      if h.get("module_type")]
                    target_interactions = new_ann_payload.get("interactions") or []
                    # Re-run verification on the corrected target.
                    target_verification = _verify_target(
                        target_sequence=target.sequence,
                        target_annotations=(
                            new_ann_payload.get("annotations") or []
                        ),
                        target_modules=target_modules,
                        target_interactions=target_interactions,
                        target_name=getattr(target, "name", "target"),
                    )
        except Exception as exc:
            logger.warning("[Pre-Design] auto-correction failed: %s", exc)

        # Build the abstract PlasmidSpec from intent + uploaded modules and
        # diff it against the assembled target so handlers (and later the LLM
        # orchestrator) can see what is missing or unexpected.
        from .plasmid_spec import build_plasmid_spec, diff_spec_against_target
        uploaded_modules = [m for pa in part_annotations for m in pa["modules"]]
        plasmid_spec = build_plasmid_spec(
            message=message,
            intent_result=intent_result,
            uploaded_modules=uploaded_modules,
        )
        spec_diff = diff_spec_against_target(plasmid_spec, target_modules)
        logger.info(
            "[Pre-Design] spec: %d required / %d present; diff: %d satisfied / %d missing",
            len(plasmid_spec.modules_required), len(plasmid_spec.modules_present),
            len(spec_diff["satisfied"]), len(spec_diff["missing"]),
        )

        # Optional LLM orchestrator pass — no-op unless PLASMID_LLM_ORCHESTRATOR=1.
        # Surfaces any orchestrator suggestions on predesign_context for the
        # per-intent handler / response builder to render.
        from .llm_orchestrator import review_design as _orchestrator_review
        orchestrator_suggestions = await _orchestrator_review(
            spec=plasmid_spec, target_modules=target_modules, spec_diff=spec_diff,
        )

        # Route to optimal workflow
        router = CloningRouter()
        objective = "balanced"
        if "cheap" in message.lower() or "cost" in message.lower():
            objective = "cost"
        elif "fast" in message.lower() or "quick" in message.lower():
            objective = "time"

        # Map intent to preferred workflow method
        intent = intent_result.get("intent") if intent_result else None
        preferred_method = None
        if intent == "gibson_design":
            preferred_method = WorkflowMethod.GIBSON
        elif intent == "golden_gate_primer_design":
            preferred_method = WorkflowMethod.GOLDEN_GATE

        logger.info(f"[Pre-Design] Intent: {intent}, preferred_method: {preferred_method}")

        candidates = await router.route(
            target,
            resolved_parts,
            objective=objective,
            preferred_method=preferred_method
        )
        best = next((c for c in candidates if c.compatible), None)

        if not best:
            logger.warning("[Pre-Design] No compatible workflows, falling back")
            return None

        logger.info(f"[Pre-Design] Selected: {best.method.value}")

        # Build pre-design context for explanation
        predesign_context = {
            "selected_workflow": best.method.value,
            "target_length": target.length,
            "target_topology": target.topology,
            "num_parts": len(resolved_parts),
            "num_compatible_workflows": len([c for c in candidates if c.compatible]),
            "selected_cost": best.total_cost_usd,
            "selected_time_days": best.total_calendar_days,
            "selected_risk": best.overall_risk_score,
            "selection_reason": best.incompatibility_reasons if not best.compatible else [],
            "alternative_workflows": [
                {
                    "method": c.method.value,
                    "compatible": c.compatible,
                    "cost": c.total_cost_usd if c.compatible else None,
                    "time_days": c.total_calendar_days if c.compatible else None,
                    "incompatibility_reasons": c.incompatibility_reasons if not c.compatible else [],
                }
                for c in candidates[:3]  # Top 3 alternatives
            ],
            "parts_summary": [
                {
                    "name": part.name,
                    "length": part.length,
                    "source": part.source_detail or "user_provided",
                }
                for part in resolved_parts
            ],
            "plasmid_spec": plasmid_spec.to_dict(),
            "spec_diff": spec_diff,
            "part_annotations": part_annotations,
            "target_modules": target_modules,
            "target_interactions": target_interactions,
            "orchestrator_suggestions": [e.__dict__ for e in orchestrator_suggestions],
            "target_verification": target_verification,
            "target_corrections": target_corrections,
        }

        # Delegate to Gibson workflow (only one integrated for now)
        if best.method == WorkflowMethod.GIBSON:
            # Step 3b: predesign also flows through the canonical
            # WorkflowInput shape so legacy + predesign paths share
            # _build_gibson_request_from_workflow_input.
            from .workflow_input import build_from_predesign
            _gibson_workflow_input = build_from_predesign(
                canonical_request=None,
                intent="gibson_design",
                session_id=session_id,
                resolved_parts=resolved_parts,
                target_plasmid=target,
                part_annotations=part_annotations,
                target_modules=target_modules,
                target_interactions=target_interactions,
                plasmid_spec=predesign_context.get("plasmid_spec"),
                spec_diff=spec_diff,
                predesign_context=predesign_context,
            )
            _gibson_workflow_input.workflow_args["assembly"] = topology
            gibson_req = _build_gibson_request_from_workflow_input(
                _gibson_workflow_input, include_explanation,
            )
            result = design_gibson_primers(gibson_req)

            reply = f"[Pre-Design] Selected Gibson Assembly for {len(resolved_parts)} fragments."

            return {
                "result": result,
                "reply": reply,
                "viz": build_gibson_viz(result),
                "files": build_gibson_files(result),
                "workflow": "gibson",
                "predesign_context": predesign_context,
                # Step 2 of workflow_input migration: surface predesign
                # internals so callers can construct a WorkflowInput via
                # workflow_input.build_from_predesign(...). Pure additive —
                # no existing reader touches these keys.
                "resolved_parts": resolved_parts,
                "target_plasmid": target,
                "part_annotations": part_annotations,
                "target_modules": target_modules,
                "target_interactions": target_interactions,
                "plasmid_spec": plasmid_spec,
                "spec_diff": spec_diff,
            }

        # Other workflows fall back to legacy
        return None

    except Exception as e:
        logger.error(f"[Pre-Design] Error: {e}", exc_info=True)
        return None


# ---------------------------------------------------------------------------
# Payload builders  (mirrors n8n "Build X payload" code nodes)
# ---------------------------------------------------------------------------

def _build_gibson_request_from_workflow_input(
    workflow_input,
    include_explanation: bool,
) -> GibsonRequest:
    """Project a Gibson WorkflowInput into a GibsonRequest.

    Step 3b of the workflow_input migration. Reads `workflow_input.parts` for
    fragments and `workflow_input.workflow_args["primer_params"]` for tuning
    knobs (the same shape as the old `intent_result["gibson_design"]
    .primer_params`, lifted into `workflow_args` by the adapter).

    Behavior parity with the legacy `_build_gibson_request(seq_data,
    intent_result, ...)` is preserved — same kwargs, same defaults, same
    "only set when LLM provided" gating.
    """
    fragments = [
        GibsonFragment(name=p.name, sequence=p.sequence) for p in workflow_input.parts
    ]
    pp = (workflow_input.workflow_args.get("primer_params") or {})
    assembly = workflow_input.workflow_args.get("assembly") or "circular"

    kwargs: Dict[str, Any] = {}
    if assembly in ("linear", "circular"):
        kwargs["assembly"] = assembly
    if pp.get("overlap_len") is not None:
        kwargs["overlap_len"] = int(pp["overlap_len"])
    if pp.get("overlap_min_len") is not None:
        kwargs["overlap_min_len"] = int(pp["overlap_min_len"])
    if pp.get("overlap_max_len") is not None:
        kwargs["overlap_max_len"] = int(pp["overlap_max_len"])
    if pp.get("overlap_target_tm") is not None:
        t = float(pp["overlap_target_tm"])
        kwargs["overlap_tm_min"] = t - 5.0
        kwargs["overlap_tm_max"] = t + 5.0
    if pp.get("overlap_min_tm") is not None:
        kwargs["overlap_tm_min"] = float(pp["overlap_min_tm"])
    if pp.get("overlap_max_tm") is not None:
        kwargs["overlap_tm_max"] = float(pp["overlap_max_tm"])
    if pp.get("anneal_target_tm") is not None:
        kwargs["anneal_tm_target"] = float(pp["anneal_target_tm"])
    if pp.get("anneal_min_len") is not None:
        kwargs["anneal_min_len"] = int(pp["anneal_min_len"])
    if pp.get("anneal_max_len") is not None:
        kwargs["anneal_max_len"] = int(pp["anneal_max_len"])

    return GibsonRequest(
        fragments=fragments,
        session_id=workflow_input.session_id,
        include_ai_explanation=include_explanation,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Main endpoint
# ---------------------------------------------------------------------------

@router.post("/api/chat")
async def chat(
    message: str = Form(""),
    session_id: str = Form(""),
    include_ai_explanation: str = Form("true"),
    describe_plasmid_intent: str = Form("false"),
    file: Optional[UploadFile] = File(None),
    inventory_files: Optional[List[UploadFile]] = File(None),
):
    """
    Main chat orchestration endpoint.

    Accepts multipart/form-data:
      message                – user's text request
      session_id             – UUID (generated if absent)
      include_ai_explanation – "true" | "false"
      describe_plasmid_intent – "true" | "false"
      file                   – optional target GenBank file
      inventory_files        – optional inventory GenBank files
    """
    if not session_id:
        session_id = str(uuid.uuid4())

    include_explanation = include_ai_explanation.strip().lower() == "true"
    has_target = file is not None
    has_inventory = bool(inventory_files)

    # ── Step 1: Extract sequences (handles multi-line, joins all DNA per label) ──
    seq_data = extract_sequences(message)

    # ── Step 2: Redact DNA so the LLM receives a compact, parseable prompt ───────
    redacted_msg = redact_sequences(message)

    # ── Step 3: Parse intent + parameters (n8n system prompt via OpenAI/fallback) ─
    intent_result = await parse_intent(
        message,
        has_target=has_target,
        has_inventory=has_inventory,
        seq_count=seq_data["count"],
        redacted_message=redacted_msg,
    )
    intent = intent_result.get("intent", "unknown")
    logger.info(f"=== Chat Request ===")
    logger.info(f"  Message: {message[:100]}...")
    logger.info(f"  has_target={has_target}, has_inventory={has_inventory}, seq_count={seq_data['count']}")
    logger.info(f"  Intent: {intent}")
    logger.info(f"  SDM params: {intent_result.get('sdm', {})}")

    # ── Step 3.5: Build canonical CloningRequest (additive; handlers opt-in) ────
    try:
        canonical_request = await normalize_request(
            message=message,
            redacted_message=redacted_msg,
            intent_result=intent_result,
            seq_data=seq_data,
            has_target=has_target,
            has_inventory=has_inventory,
        )
        logger.info(
            f"  Canonical: parts={len(canonical_request.parts)} "
            f"vector={canonical_request.vector.name if canonical_request.vector else None} "
            f"enzymes={[e.name for e in canonical_request.enzymes]} "
            f"mutations={len(canonical_request.mutations)} "
            f"constraints={[(c.kind, c.target) for c in canonical_request.constraints]}"
        )
    except Exception as exc:
        logger.warning(f"  Canonical request normalization failed: {exc}")
        canonical_request = None


    # ── Target-from-Inventory pre-router (now scoring-only as of 2026-04-27) ──
    # The chat-level `target_from_inventory` intent was retired; the router
    # is still invoked here whenever target+inventory are both present so
    # downstream workflows can include the per-method scorecard in their reply.
    _routing_audit_md: str = ""
    _routing_audit_payload: Optional[Dict[str, Any]] = None
    _routing_context_for_llm: Optional[Dict[str, Any]] = None

    # Universal pre-router: when both a target and inventory are present,
    # evaluate every cloning workflow once up front and cache the reports.
    # Every downstream workflow branch then appends a "Method Selection
    # Context" section to its reply so the user always sees the scorecard
    # + rationale behind the chosen path, even when they named the method
    # explicitly. The chat-level `target_from_inventory` intent was removed
    # 2026-04-27; the underlying router (`target_from_inventory_router.py`)
    # is still reachable here for scoring context only.
    _WORKFLOWS_WITH_ROUTING_CONTEXT = {
        "gibson_design", "golden_gate_primer_design", "gateway_cloning",
        "sdm_design", "restriction_cloning",
    }

    async def _maybe_run_universal_router():
        """Run route_from_uploads once if target+inventory both present.
        Returns (routing_dict, routing_ctx_for_llm) or (None, None).
        """
        if not (has_target and has_inventory):
            return None, None
        try:
            from .target_from_inventory_router import route_from_uploads as _rfu
            _routing = await _rfu(file, inventory_files or [])
            _reports = _routing.get("reports") or []
            _chosen_wf = _routing.get("chosen_intent")
            _chosen_report = next(
                (r for r in _reports if r.get("workflow") == _chosen_wf), None
            )
            _ctx = {
                "chosen_workflow": _chosen_wf,
                "chosen_rationale": (_chosen_report or {}).get("rationale", ""),
                "chosen_score": (_chosen_report or {}).get("score"),
                "chosen_validation_mode": _routing.get("chosen_validation_mode"),
                "chosen_work_estimate": (_chosen_report or {}).get("work_estimate"),
                "chosen_success_estimate": (_chosen_report or {}).get("success_estimate"),
                "alternatives": [
                    {
                        "workflow": r.get("workflow"),
                        "feasible": r.get("feasible"),
                        "score": r.get("score"),
                        "rationale": (r.get("rationale") or "")[:240],
                    }
                    for r in _reports if r.get("workflow") != _chosen_wf
                ],
                "audit_markdown": _routing.get("audit_markdown", ""),
            }
            return _routing, _ctx
        except Exception as _re:
            logger.warning(f"[router] pre-evaluation failed: {_re}")
            return None, None

    def _method_selection_markdown(ctx: Optional[Dict[str, Any]], selected_intent: str) -> str:
        """Render the router-scored workflow table as a deterministic markdown
        block. Appended to each workflow reply so users see the alternatives.
        Called from every branch in _WORKFLOWS_WITH_ROUTING_CONTEXT.
        """
        if not ctx:
            return ""
        alts = ctx.get("alternatives") or []
        lines = ["", "---", "", "### Method Selection Context", ""]
        _cs = ctx.get("chosen_score")
        _ce = ctx.get("chosen_success_estimate")
        _cw = ctx.get("chosen_work_estimate")
        _cv = ctx.get("chosen_validation_mode") or ""
        _chosen_wf = ctx.get("chosen_workflow") or selected_intent
        lines.append(
            f"**Router pick:** `{_chosen_wf}` "
            + (f"(score {_cs:.3f}" if isinstance(_cs, (int, float)) else "(")
            + (f", success {_ce:.2f}" if isinstance(_ce, (int, float)) else "")
            + (f", work {_cw}" if _cw is not None else "")
            + (f", validation {_cv})" if _cv else ")")
        )
        if _chosen_wf != selected_intent:
            lines.append(
                f"_You invoked `{selected_intent}` directly; the router would have picked "
                f"`{_chosen_wf}` based on structural evaluation of your target + inventory._"
            )
        _rat = ctx.get("chosen_rationale") or ""
        if _rat:
            lines.append(f"**Rationale:** {_rat}")
        if alts:
            lines.append("")
            lines.append("**Alternative workflows evaluated:**")
            lines.append("")
            lines.append("| Workflow | Feasible | Score | Rationale |")
            lines.append("|---|---|---|---|")
            _ordered = sorted(alts, key=lambda a: (-int(bool(a.get("feasible"))),
                                                     -(a.get("score") or 0.0)))
            for a in _ordered:
                _s = a.get("score")
                _ss = f"{_s:.3f}" if isinstance(_s, (int, float)) else "—"
                _feas = "✅" if a.get("feasible") else "❌"
                _wf = a.get("workflow", "?")
                _ar = (a.get("rationale") or "").replace("|", "\\|")[:120]
                lines.append(f"| {_wf} | {_feas} | {_ss} | {_ar} |")
        return "\n".join(lines)

    # Run the universal router once up front (cached for all branches).
    if has_target and has_inventory:
        _routing_audit_payload, _routing_context_for_llm = await _maybe_run_universal_router()
        if _routing_audit_payload:
            _routing_audit_md = _routing_audit_payload.get("audit_markdown", "")

    # Intents whose existing /api/chat branches can execute directly from file uploads.
    _DISPATCHABLE_FROM_UPLOADS = {
        "gateway_cloning", "sdm_design", "restriction_cloning",
        "sgrna_golden_gate", "annotate_gb",
    }

    # Auto-route on implicit intent: if the user uploaded target+inventory but
    # did not name a method, the universal pre-router already scored every
    # workflow above. Promote the winning workflow to `intent` so the matching
    # branch dispatches without a separate `target_from_inventory` chat label.
    # Restored 2026-04-27 (the explicit chat intent stays removed).
    _IMPLICIT_INTENTS = {"unknown", "plasmid_design", ""}
    if (
        has_target and has_inventory
        and intent in _IMPLICIT_INTENTS
        and _routing_audit_payload
    ):
        _chosen = _routing_audit_payload.get("chosen_intent")
        _chosen_args = _routing_audit_payload.get("chosen_handler_args") or {}
        if _chosen in _DISPATCHABLE_FROM_UPLOADS:
            logger.info(
                f"[auto-route] implicit intent={intent!r} + target+inventory "
                f"-> promoting to chosen_intent={_chosen!r}"
            )
            intent = _chosen
            if _chosen == "sdm_design":
                sdm = intent_result.setdefault("sdm", {})
                if _chosen_args.get("mutation_type"):
                    sdm["mutation_type"] = _chosen_args["mutation_type"]
                if _chosen_args.get("target_position_start") is not None:
                    sdm["target_position_start"] = _chosen_args["target_position_start"]
                if _chosen_args.get("target_position_end") is not None:
                    sdm["target_position_end"] = _chosen_args["target_position_end"]
            elif _chosen == "sgrna_golden_gate":
                sgrna = intent_result.setdefault("sgrna", {})
                if _chosen_args.get("enzyme"):
                    sgrna["enzyme"] = _chosen_args["enzyme"]
            elif _chosen == "gateway_cloning":
                intent_result.setdefault("gateway", {}).update({
                    "reaction_type": _chosen_args.get("reaction_type"),
                    "gateway_variant": _chosen_args.get("gateway_variant"),
                })
            elif _chosen == "restriction_cloning":
                rc = intent_result.setdefault("restriction_cloning", {})
                if _chosen_args.get("enzymes"):
                    rc["enzymes"] = _chosen_args["enzymes"]
        else:
            logger.info(
                f"[auto-route] chosen_intent={_chosen!r} not dispatchable from uploads; "
                "leaving intent unchanged (audit will be appended to reply)."
            )
            if _chosen == "synthesis_fallback":
                _manifest = _chosen_args
                _routing_audit_md += (
                    "\n\n### Synthesis manifest\n\n"
                    f"- Total synth: **{_manifest.get('total_synth_bp', 0)} bp** in "
                    f"**{len(_manifest.get('synth_blocks', []))} block(s)**\n"
                    f"- Estimated cost: **${_manifest.get('est_cost_usd', 0):.0f}**\n"
                    f"- Inventory anchors: {len(_manifest.get('inventory_anchors', []))}\n"
                )

    try:
        result: Optional[Dict[str, Any]] = None
        viz: Optional[Dict[str, Any]] = None
        viz_list: Optional[List[Dict[str, Any]]] = None
        files: List[Dict[str, str]] = []
        reply = ""
        ident_prefix = ""  # Identification summary prepended to final reply

        # ── Unified Pre-Design System (Feature Flag) ──────────────────────────
        # Check if unified pre-design is enabled
        USE_UNIFIED_PREDESIGN = os.getenv("ENABLE_UNIFIED_PREDESIGN", "true").lower() == "true"

        # List of intents compatible with unified pre-design
        predesign_compatible_intents = [
            "gibson_design",
            "golden_gate_primer_design",
            # SDM has special handling, keep in legacy for now
        ]

        if USE_UNIFIED_PREDESIGN and intent in predesign_compatible_intents:
            logger.info("[Pre-Design] Attempting unified pre-design pipeline")
            predesign_result = await _execute_unified_predesign(
                message,
                seq_data,
                intent_result,
                file,
                inventory_files,
                session_id,
                include_explanation,
            )

            # If pre-design succeeded, use its result
            if predesign_result:
                logger.info("[Pre-Design] Using unified pre-design result")
                result = predesign_result.get("result")
                reply = predesign_result.get("reply", "")
                viz = predesign_result.get("viz")
                files = predesign_result.get("files", [])
                predesign_context = predesign_result.get("predesign_context")

                # Add explanation if requested
                if include_explanation and result:
                    explanation = await generate_explanation(
                        intent=intent_result.get("intent", "gibson_design"),
                        result=result,
                        user_message=message,
                        predesign_context=predesign_context,
                    )
                    if explanation:
                        reply += f"\n\n{explanation}"

                # Route the predesign result through normalize_response so the
                # universal {assembled.gb, parts_order.csv, protocol.csv,
                # workflow_trace.txt} files emit on this path too. Falls back
                # to the legacy raw-files dict on builder failure.
                try:
                    from .workflow_input import build_from_predesign
                    from .output_builders import build_gibson_output
                    from .output_normalizer import normalize_response
                    _wf_in_pre = build_from_predesign(
                        canonical_request=canonical_request,
                        intent=intent,
                        session_id=session_id,
                        resolved_parts=predesign_result.get("resolved_parts") or [],
                        target_plasmid=predesign_result.get("target_plasmid"),
                        part_annotations=predesign_result.get("part_annotations"),
                        target_modules=predesign_result.get("target_modules"),
                        target_interactions=predesign_result.get("target_interactions"),
                        plasmid_spec=(predesign_context or {}).get("plasmid_spec"),
                        spec_diff=predesign_result.get("spec_diff"),
                        predesign_context=predesign_context,
                    )
                    if intent == "gibson_design":
                        _wf_out_pre = build_gibson_output(
                            workflow_input=_wf_in_pre,
                            result=result,
                            llm_summary=reply,
                        )
                        _resp = normalize_response(_wf_out_pre, _wf_in_pre)
                        _resp.setdefault("intent", intent)
                        return _resp
                except Exception as _norm_exc:
                    logger.warning(
                        "[Pre-Design][output normalize] %s — falling back",
                        _norm_exc, exc_info=True,
                    )

                return {
                    "ok": True,
                    "reply": reply,
                    "viz": viz,
                    "viz_list": viz_list,
                    "files": files,
                    "sessionId": session_id,
                }
            else:
                logger.info("[Pre-Design] Falling back to legacy workflow")
                # Continue to legacy code below

        # ── Gateway Cloning ────────────────────────────────────────────────────
        if intent == "gateway_cloning":
            from .cloning.gateway_operator import GatewayOperator
            from .cloning.gateway_sites import scan_att_sites, GATEWAY_ATT_SITES

            logger.info("[Gateway] Processing Gateway cloning request")

            # Step 3h of workflow_input migration: file resolution + per-file
            # annotation (Step 2.75 already produced gateway_att hits) + KB
            # lookup of message-named inserts all live in
            # build_for_gateway_cloning. The handler reads the
            # GatewayOperator-shape modules from workflow_args and the
            # per-module att-site subtypes from each part's metadata, so
            # scan_att_sites no longer runs once per module.
            from .workflow_input_adapters import build_for_gateway_cloning
            workflow_input = await build_for_gateway_cloning(
                message=message,
                file=file,
                inventory_files=inventory_files,
                intent_result=intent_result,
                canonical_request=canonical_request,
                session_id=session_id,
            )
            modules = workflow_input.workflow_args["modules"]
            vector_found = workflow_input.workflow_args["vector_found"]
            insert_found = workflow_input.workflow_args["insert_found"]
            gateway_identified_inserts = (
                workflow_input.provenance.get("gateway_identified_inserts") or []
            )
            for ident in gateway_identified_inserts:
                logger.info(
                    "[Gateway] Added insert from pLannotate KB: %s -> %s (%d bp)",
                    ident["query"], ident["name"], ident["length"],
                )

            # Validate we have what we need
            if not vector_found:
                reply = "Gateway cloning requires a donor vector (like pDONR221) with attP sites. Please upload or reference a Gateway-compatible vector in your inventory."
                result = {"error": "no_gateway_vector"}
            elif not insert_found:
                reply = "Please specify the insert sequence. You can:\n- Upload a plasmid file\n- Reference a feature name (e.g., 'GFP', 'mCherry')\n- Provide a DNA sequence directly"
                result = {"error": "no_insert"}
            else:
                # Run Gateway operator
                logger.info(f"[Gateway] Running operator with {len(modules)} modules")
                operator = GatewayOperator()
                plan = operator.evaluate(modules, topology="circular")

                # 2026-05-12: orientation auto-correct via the unified
                # predesign helper (`auto_correct_kb_part_orientation`).
                # The same routine that runs inside `_execute_unified_predesign`
                # for Gibson/Golden Gate now runs against the Gateway product
                # — KB-resolved parts in the wrong orientation get
                # reverse-complemented before the assembled.gb is emitted.
                # The general `_check_cds_functional_context` orientation
                # detector flags any EXPRESSION_FEATURE_KINDS feature whose
                # natural partner sits on the opposite strand within window,
                # not just T7+CDS pairs.
                plan._cached_full = None
                if plan.feasible and plan.product_sequence:
                    try:
                        from .annotation_cache import annotate_llm_cached as _alc
                        from .target_from_inventory_router import (
                            auto_correct_kb_part_orientation as _autocorrect,
                        )

                        # Build kb-resolved part names — the auto-correct only
                        # touches features whose name matches one of these.
                        _kb_part_names = []
                        for _ident in (intent_result or {}).get("kb_resolved", {}).get("identified", []) or []:
                            _nm = _ident.get("feature_name") or _ident.get("name")
                            if _nm:
                                _kb_part_names.append(_nm)
                        for _m in modules:
                            _cid = _m.get("canonical_id") or _m.get("name")
                            if _cid and _m.get("role") == "insert":
                                _kb_part_names.append(_cid)

                        _full0 = await _alc(plan.product_sequence, circular=True)
                        plan._cached_full = _full0

                        async def _reanno(seq):
                            payload = await _alc(seq, circular=True)
                            return (
                                payload.get("annotations") or [],
                                [h for h in (payload.get("hierarchical_annotations") or [])
                                 if h.get("module_type")],
                                payload.get("interactions") or [],
                            )

                        if _kb_part_names:
                            _ac = await _autocorrect(
                                target_sequence=plan.product_sequence,
                                target_annotations=_full0.get("annotations") or [],
                                target_modules=[h for h in (_full0.get("hierarchical_annotations") or [])
                                                 if h.get("module_type")],
                                target_interactions=_full0.get("interactions") or [],
                                kb_resolved_part_names=_kb_part_names,
                                reannotate_async=_reanno,
                            )
                            _applied = [c for c in (_ac.get("corrections") or [])
                                        if c.get("result") == "applied"]
                            if _ac.get("corrected_sequence") and _applied:
                                logger.info(
                                    "[Gateway] auto-correct applied %d correction(s): %s",
                                    len(_applied),
                                    ", ".join(c["feature"] for c in _applied),
                                )
                                plan.product_sequence = _ac["corrected_sequence"]
                                # Re-annotate the corrected product for the viz layer.
                                plan._cached_full = await _alc(plan.product_sequence, circular=True)
                                for _c in _applied:
                                    plan.warnings.append(
                                        f"[orientation] Auto-reverse-complemented '{_c['feature']}' "
                                        f"({_c['action']}); see workflow_trace for details."
                                    )
                            else:
                                logger.info(
                                    "[Gateway] auto-correct found nothing to fix (KB names checked: %s)",
                                    _kb_part_names,
                                )
                    except Exception as _exc:
                        logger.warning("[Gateway] auto-correct skipped: %s", _exc, exc_info=True)

                if not plan.feasible:
                    reply = f"Gateway design failed: {'; '.join(plan.infeasibility_reasons)}"
                    result = {"error": "infeasible", "reasons": plan.infeasibility_reasons}
                else:
                    # Build response
                    reply_parts = []
                    reply_parts.append("## Gateway Cloning Design\n\n")
                    reply_parts.append(f"**Reaction Type:** {plan.reaction_type}\n")
                    reply_parts.append(f"**Modules:** {plan.fragment_count}\n\n")
                    
                    # Show input modules
                    reply_parts.append("### Input Modules\n\n")
                    for i, mod in enumerate(modules, 1):
                        reply_parts.append(f"{i}. **{mod['canonical_id']}** ({mod['role']})\n")
                        reply_parts.append(f"   - {mod['description']}\n")

                        # Read pre-scanned gateway_att subtypes from the
                        # adapter (Step 2.75 cloning_features). KB-resolved
                        # inserts have an empty att_subtypes list.
                        _part = workflow_input.part_by_name(mod["canonical_id"])
                        site_names = (_part.metadata.get("att_subtypes") or []) if _part else []
                        if site_names:
                            reply_parts.append(f"   - att sites: {', '.join(site_names)}\n")
                        reply_parts.append("\n")

                    # Junction information
                    reply_parts.append("### Recombination\n\n")
                    for jp in plan.junction_plans:
                        reply_parts.append(f"**Junction {jp.junction_index}:** {jp.left_module_name} → {jp.right_module_name}\n")
                        reply_parts.append(f"- Strategy: {jp.strategy}\n")
                        reply_parts.append(f"- att sites: {jp.left_att_site} + {jp.right_att_site} → {jp.product_left_site} + {jp.product_right_site}\n")
                        
                        if jp.warnings:
                            for warning in jp.warnings:
                                reply_parts.append(f"- ⚠️ {warning}\n")
                        reply_parts.append("\n")

                    # Primers if needed. If any identified insert is
                    # back-translated (no physical template), lead with
                    # a synthesis-first note and demote the primer header
                    # to a fallback path.
                    _gateway_synthesis_required = any(
                        (f.get("sequence_provenance") or "").startswith("backtranslated_")
                        for f in gateway_identified_inserts
                    )
                    if plan.primer_table:
                        if _gateway_synthesis_required:
                            reply_parts.append("### Synthesis-first path (recommended)\n\n")
                            reply_parts.append(
                                "At least one insert was back-translated from the protein KB "
                                "(no physical template exists). Order the attB-flanked fragment "
                                "shown in `gateway_synthesis_fragment.csv` as a gBlock / eBlock "
                                "and run BP clonase directly — **no PCR required**.\n\n"
                            )
                            reply_parts.append("### Alternative PCR path (fallback)\n\n")
                            reply_parts.append(
                                "If a physical template is acquired later, the following attB-tailed "
                                "primers would amplify it for BP clonase:\n\n"
                            )
                        else:
                            reply_parts.append("### Primers for attB Site Addition\n\n")
                            reply_parts.append("Design PCR primers to add Gateway attB sites to your insert:\n\n")
                        
                        for primer in plan.primer_table:
                            reply_parts.append(f"**{primer['primer_name']}**\n")
                            reply_parts.append(f"\n")
                            reply_parts.append(f"- Sequence: `{primer.get('sequence', 'N/A')}`\n")
                            if primer.get("att_site_tail") and primer.get("annealing_region"):
                                reply_parts.append(f"  - GGGG tail (4 bp) + attB site ({len(primer['att_site_tail'])} bp) + annealing region ({len(primer['annealing_region'])} bp)\n")
                            reply_parts.append(f"- Length: {primer['length']} bp\n")
                            reply_parts.append(f"- Tm (annealing): {primer.get('tm_anneal', 60):.1f}°C\n")
                            reply_parts.append(f"- Purpose: {primer['purpose']}\n\n")

                    # Product information
                    reply_parts.append("### Expected Product\n\n")
                    reply_parts.append(f"**{plan.reaction_type} Reaction Product:**\n\n")
                    reply_parts.append(f"- Size: {len(plan.product_sequence):,} bp\n")
                    
                    product_sites = scan_att_sites(plan.product_sequence, fuzzy_threshold=0)
                    if product_sites:
                        product_site_names = [s.site_type for s in product_sites]
                        reply_parts.append(f"- att sites: {', '.join(set(product_site_names))}\n")
                    
                    from .cloning.gateway_sites import scan_for_ccdb
                    product_ccdb = scan_for_ccdb(plan.product_sequence)
                    if len(product_ccdb) == 0:
                        reply_parts.append(f"- ✅ ccdB gene removed (Entry clone)\n")
                    else:
                        reply_parts.append(f"- ⚠️ ccdB gene present (requires DB3.1 strain)\n")

                    # Warnings
                    if plan.warnings:
                        reply_parts.append("\n### Warnings\n\n")
                        for warning in plan.warnings:
                            reply_parts.append(f"- ⚠️ {warning}\n")

                    reply = "".join(reply_parts)

                    # Prepend KB-identification audit line if we resolved any insert via the pLannotate KB.
                    if gateway_identified_inserts:
                        ident_lines = ["Identified insert(s) from pLannotate feature KB:"]
                        for f in gateway_identified_inserts:
                            ident_lines.append(
                                f"  - {f['query']} -> {f['name']} ({f['feature_type']}, {f['length']} bp) [{f['feature_id']}]"
                            )
                        reply = "\n".join(ident_lines) + "\n\n" + reply


                    # Generate PCR amplicon GenBank file
                    # This is the product of PCR with attB primers
                    pcr_amplicon_files = []

                    # Find the insert module (the one that gets primers).
                    # No att sites = insert; read the pre-scanned subtypes
                    # from each part's metadata.
                    insert_module = None
                    for mod in modules:
                        _p = workflow_input.part_by_name(mod["canonical_id"])
                        _subs = (_p.metadata.get("att_subtypes") or []) if _p else []
                        if not _subs:
                            insert_module = mod
                            break

                    if insert_module and plan.primer_table:
                        from Bio.Seq import Seq
                        from Bio.SeqRecord import SeqRecord
                        from Bio.SeqFeature import SeqFeature, FeatureLocation
                        import tempfile

                        insert_name = insert_module.get("canonical_id", "insert")
                        insert_seq = insert_module.get("sequence", "")

                        # Get primers
                        fwd_primer = None
                        rev_primer = None
                        for primer in plan.primer_table:
                            if "FWD" in primer.get("primer_name", ""):
                                fwd_primer = primer
                            elif "REV" in primer.get("primer_name", ""):
                                rev_primer = primer

                        if fwd_primer and rev_primer:
                            # Build PCR amplicon sequence
                            # Structure: GGGG + attB1 + insert + attB2_RC + GGGG
                            gggg_tail = "GGGG"
                            attB1 = fwd_primer.get("att_site_tail", "")
                            attB2_rc = rev_primer.get("att_site_tail", "")

                            # Extract RBS from forward primer if present
                            # Primer structure: GGGG + attB + CC + RBS + annealing
                            # RBS is between CC and annealing region
                            fwd_primer_seq = fwd_primer.get("sequence", "")
                            rbs_fwd = ""
                            
                            # Expected length without RBS: 4 (GGGG) + 25 (attB) + 2 (CC) + 20 (annealing) = 51bp
                            # With RBS: 51 + len(RBS)
                            if len(fwd_primer_seq) > 51:
                                # There is RBS - extract it
                                # RBS is after GGGG (4) + attB (25) + CC (2) = 31bp
                                rbs_start = 31
                                rbs_end = len(fwd_primer_seq) - 20  # Before annealing region
                                if rbs_end > rbs_start:
                                    rbs_fwd = fwd_primer_seq[rbs_start:rbs_end]
                            
                            # Frame-maintaining nucleotides
                            cc_spacer = "CC"  # After attB1
                            c_spacer = "C"    # Before attB2
                            
                            # Build PCR amplicon sequence with RBS and spacers
                            # Structure: GGGG + attB1 + CC + RBS + insert + C + attB2_RC + GGGG
                            pcr_amplicon_seq = gggg_tail + attB1 + cc_spacer + rbs_fwd + insert_seq + c_spacer + attB2_rc + gggg_tail

                            # Create SeqRecord
                            # Sanitize name for GenBank format (no spaces allowed)
                            insert_name_clean = insert_name.replace(" ", "_").replace("-", "_")
                            
                            amplicon_record = SeqRecord(
                                Seq(pcr_amplicon_seq),
                                id=f"{insert_name_clean}_PCR_amplicon",
                                name=f"{insert_name_clean}_PCR",
                                description=f"PCR amplicon of {insert_name} with Gateway attB sites"
                            )

                            # Set required GenBank annotations
                            amplicon_record.annotations["molecule_type"] = "DNA"
                            amplicon_record.annotations["topology"] = "linear"

                            # Add features
                            pos = 0

                            # 5' GGGG tail
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(pos, pos + 4),
                                type="misc_feature",
                                qualifiers={"label": "5' GGGG tail",
                                          "note": "Recommended by Gateway manual for PCR efficiency",
                                          "ApEinfo_fwdcolor": "#ff9999"}
                            ))
                            pos += 4

                            # attB1 site
                            attB1_start = pos
                            attB1_end = pos + len(attB1)
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(attB1_start, attB1_end),
                                type="protein_bind",
                                qualifiers={"label": "attB1",
                                          "note": "Gateway recombination site",
                                          "ApEinfo_fwdcolor": "#ff6b9d"}
                            ))
                            pos += len(attB1)

                            # CC frame spacer
                            spacer_cc_start = pos
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(spacer_cc_start, spacer_cc_start + 2),
                                type="misc_feature",
                                qualifiers={"label": "frame spacer (CC)",
                                          "note": "Frame-maintaining nucleotides (avoids stop codons)",
                                          "ApEinfo_fwdcolor": "#dddddd"}
                            ))
                            pos += 2

                            # RBS (Ribosome Binding Site) if present
                            rbs_start = pos
                            if rbs_fwd:
                                rbs_end = pos + len(rbs_fwd)
                                # Combined Shine-Dalgarno/Kozak for native expression
                                rbs_type = "RBS (Shine-Dalgarno/Kozak)"
                                rbs_note = "Native ribosome binding site for both bacterial and eukaryotic expression"

                                amplicon_record.features.append(SeqFeature(
                                    FeatureLocation(rbs_start, rbs_end),
                                    type="regulatory",
                                    qualifiers={"label": rbs_type,
                                              "note": rbs_note,
                                              "ApEinfo_fwdcolor": "#ffcc99"}
                                ))
                                pos += len(rbs_fwd)

                            # Insert (template)
                            insert_start = pos
                            insert_end = pos + len(insert_seq)
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(insert_start, insert_end),
                                type="CDS",
                                qualifiers={"label": insert_name,
                                          "note": "Insert from template",
                                          "ApEinfo_fwdcolor": "#85dae9"}
                            ))

                            # FWD primer annealing site (within insert)
                            anneal_len = len(fwd_primer.get("annealing_region", ""))
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(insert_start, insert_start + anneal_len),
                                type="primer_bind",
                                qualifiers={"label": f"{fwd_primer.get('primer_name', 'FWD')} annealing",
                                          "note": f"Forward primer binding site ({anneal_len} bp)",
                                          "ApEinfo_fwdcolor": "#a6d96a"}
                            ))

                            # REV primer annealing site (within insert)
                            rev_anneal_len = len(rev_primer.get("annealing_region", ""))
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(insert_end - rev_anneal_len, insert_end),
                                type="primer_bind",
                                qualifiers={"label": f"{rev_primer.get('primer_name', 'REV')} annealing",
                                          "note": f"Reverse primer binding site ({rev_anneal_len} bp)",
                                          "ApEinfo_fwdcolor": "#fdae61"}
                            ))

                            pos += len(insert_seq)

                            # C frame spacer
                            spacer_c_start = pos
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(spacer_c_start, spacer_c_start + 1),
                                type="misc_feature",
                                qualifiers={"label": "frame spacer (C)",
                                          "note": "Frame-maintaining nucleotide",
                                          "ApEinfo_fwdcolor": "#dddddd"}
                            ))
                            pos += 1

                            # attB2 site (reverse complement)
                            attB2_start = pos
                            attB2_end = pos + len(attB2_rc)
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(attB2_start, attB2_end),
                                type="protein_bind",
                                qualifiers={"label": "attB2",
                                          "note": "Gateway recombination site (reverse strand)",
                                          "ApEinfo_fwdcolor": "#ff6b9d"}
                            ))
                            pos += len(attB2_rc)

                            # 3' GGGG tail
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(pos, pos + 4),
                                type="misc_feature",
                                qualifiers={"label": "3' GGGG tail",
                                          "note": "Recommended by Gateway manual for PCR efficiency",
                                          "ApEinfo_fwdcolor": "#ff9999"}
                            ))

                            # FWD primer extension (GGGG + attB1 + CC + RBS)
                            fwd_ext_end = rbs_start + len(rbs_fwd) if rbs_fwd else spacer_cc_start + 2
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(0, fwd_ext_end),
                                type="misc_feature",
                                qualifiers={"label": "FWD primer extension",
                                          "note": f"Added by {fwd_primer.get('primer_name', 'FWD primer')} (GGGG + attB1 + CC + RBS)",
                                          "ApEinfo_fwdcolor": "#b3e2cd"}
                            ))

                            # REV primer extension (C + attB2 + GGGG)
                            rev_ext_start = spacer_c_start
                            amplicon_record.features.append(SeqFeature(
                                FeatureLocation(rev_ext_start, pos + 4),
                                type="misc_feature",
                                qualifiers={"label": "REV primer extension",
                                          "note": f"Added by {rev_primer.get('primer_name', 'REV primer')} (C + attB2 + GGGG)",
                                          "ApEinfo_fwdcolor": "#fdcdac"}
                            ))


                            # Write to temp file
                            with tempfile.NamedTemporaryFile(mode='w', suffix='.gb', delete=False) as tmp:
                                from Bio import SeqIO as BioSeqIO
                                BioSeqIO.write(amplicon_record, tmp, "genbank")
                                tmp_path = tmp.name

                            # Read back as string
                            with open(tmp_path, 'r') as f:
                                gb_content = f.read()

                            os.unlink(tmp_path)

                            # Add to files list
                            pcr_amplicon_files.append({
                                "fileName": f"{insert_name_clean}_PCR_amplicon.gb",
                                "mimeType": "application/octet-stream",
                                "dataBase64": base64.b64encode(gb_content.encode("utf-8")).decode("ascii")
                            })

                            # Synthesis-ready CSV — the amplicon sequence is
                            # already the gBlock the user would order. List it
                            # as row 0, then the PCR primers beneath so both
                            # paths are visible in one place.
                            synth_csv_lines = [
                                "Row,Name,Length,Sequence,Notes",
                                (
                                    f"1,{insert_name_clean}_synthesis_fragment,"
                                    f"{len(pcr_amplicon_seq)},{pcr_amplicon_seq},"
                                    "\"Synthesis-ready linear fragment: GGGG + attB1 + insert + attB2 + GGGG. "
                                    "Order as gBlock/eBlock and run BP clonase directly — no PCR required.\""
                                ),
                            ]
                            if fwd_primer:
                                fwd_seq = fwd_primer.get("sequence", "")
                                synth_csv_lines.append(
                                    f"2,{fwd_primer.get('primer_name','FWD')},"
                                    f"{len(fwd_seq)},{fwd_seq},"
                                    "\"attB1-tailed forward primer (PCR path)\""
                                )
                            if rev_primer:
                                rev_seq = rev_primer.get("sequence", "")
                                synth_csv_lines.append(
                                    f"3,{rev_primer.get('primer_name','REV')},"
                                    f"{len(rev_seq)},{rev_seq},"
                                    "\"attB2-tailed reverse primer (PCR path)\""
                                )
                            synth_csv = "\n".join(synth_csv_lines)
                            pcr_amplicon_files.append({
                                "fileName": "gateway_synthesis_fragment.csv",
                                "mimeType": "text/csv",
                                "dataBase64": base64.b64encode(synth_csv.encode("utf-8")).decode("ascii"),
                            })

                            logger.info(f"[Gateway] Generated PCR amplicon file: {len(pcr_amplicon_seq)} bp")

                    # Merge PCR amplicon files with any other files
                    files.extend(pcr_amplicon_files)

                    # Always emit a Gateway product GenBank for download so
                    # the user can inspect the assembled plasmid even when
                    # the response viz fails to render.
                    if plan.product_sequence:
                        try:
                            from .files import _make_genbank as _gw_make_gb
                            # Merge re-annotated bare features (post-auto-correct)
                            # with the fragment / att overlay annotations so the
                            # downloaded .gb shows the corrected GFP orientation.
                            _gb_anns = []
                            _full_for_gb = getattr(plan, "_cached_full", None)
                            if _full_for_gb:
                                for a in (_full_for_gb.get("annotations") or []):
                                    if (a.get("layer") or "") == "cloning_feature":
                                        continue
                                    _gb_anns.append(a)
                            for a in (plan.product_annotations or []):
                                _gb_anns.append(a)
                            _gw_gb = _gw_make_gb(
                                seq=plan.product_sequence,
                                name=f"{modules[0]['canonical_id']}_Gateway_Product"[:16],
                                description=f"Gateway {plan.reaction_type} product",
                                annotations=_gb_anns,
                                topology="circular",
                            )
                            files.append({
                                "fileName": f"{modules[0]['canonical_id']}_gateway_product.gb",
                                "mimeType": "application/octet-stream",
                                "dataBase64": base64.b64encode(_gw_gb).decode("ascii"),
                            })
                        except Exception as _exc:
                            logger.warning("[Gateway] product .gb emit failed: %s", _exc)

                    # Build result for visualization
                    result = {
                        "method": "gateway_cloning",
                        "reaction_type": plan.reaction_type,
                        "feasible": True,
                        "product_sequence": plan.product_sequence,
                        "product_size": len(plan.product_sequence),
                        "primers": plan.primer_table,
                        "junctions": [
                            {
                                "junction_index": jp.junction_index,
                                "left_module": jp.left_module_name,
                                "right_module": jp.right_module_name,
                                "left_att_site": jp.left_att_site,
                                "right_att_site": jp.right_att_site,
                                "product_left_site": jp.product_left_site,
                                "product_right_site": jp.product_right_site,
                                "strategy": jp.strategy
                            }
                            for jp in plan.junction_plans
                        ],
                        "annotations": plan.product_annotations if plan.product_annotations else []
                    }

                    # 2026-05-06: redesigned Gateway viz layering.
                    # Base layer (GRAY): full re-annotation of the product —
                    #   bare features + modules, excluding cloning_features
                    #   (restriction sites / PCR warnings). Each entry is
                    #   tinted gray so the structural map is legible without
                    #   competing with the colored fragment overlays.
                    # Overlay layer (COLORED): plan.product_annotations,
                    #   which carries fragment_annotations (insert /
                    #   donor_backbone) + att-site overlays + ccdB warnings
                    #   in their distinct colors.
                    _GRAY_BASE_FEATURE = "#9CA3AF"   # mid-gray (Tailwind gray-400)
                    _GRAY_BASE_MODULE = "#D1D5DB"    # lighter gray for modules
                    gateway_product_annotations = []
                    _full = None
                    if plan.product_sequence:
                        try:
                            # 2026-05-06: reuse the post-auto-correct _full
                            # cached on plan (avoids a second annotate call
                            # of the corrected sequence). Falls back to a
                            # fresh annotate_llm_cached if the cache wasn't
                            # populated.
                            from .annotation_cache import annotate_llm_cached, annotate_cached
                            _full = getattr(plan, "_cached_full", None)
                            if _full is None:
                                _full = await annotate_llm_cached(
                                    plan.product_sequence, circular=True,
                                )
                            # Bare features (annotations field) — gray base
                            for a in (_full.get("annotations") or []):
                                if (a.get("layer") or "") == "cloning_feature":
                                    continue
                                gateway_product_annotations.append({
                                    "name": a.get("name") or a.get("type") or "feature",
                                    "start": int(a.get("start", 0)),
                                    "end": int(a.get("end", 0)),
                                    "direction": -1 if (a.get("direction") == -1 or a.get("strand") == -1) else 1,
                                    "type": a.get("type"),
                                    "color": _GRAY_BASE_FEATURE,
                                    "layer": "feature",
                                    "metadata": {"viz_layer": "base_grayed"},
                                })
                            # Module-level annotations — slightly lighter gray.
                            # The _llm endpoint surfaces rule-based modules in
                            # hierarchical_annotations entries that carry a
                            # module_type tag.
                            _mods = [h for h in (_full.get("hierarchical_annotations") or [])
                                     if h.get("module_type")]
                            for m in _mods:
                                if (m.get("layer") or "") == "cloning_feature":
                                    continue
                                gateway_product_annotations.append({
                                    "name": m.get("name") or m.get("module_type") or "feature",
                                    "start": int(m.get("start", 0)),
                                    "end": int(m.get("end", 0)),
                                    "direction": -1 if m.get("strand") == -1 else 1,
                                    "module_type": m.get("module_type"),
                                    "color": _GRAY_BASE_MODULE,
                                    "layer": "module",
                                    "metadata": {"viz_layer": "base_grayed"},
                                })
                        except Exception as _exc:
                            logger.warning("[Gateway] full re-annotation failed: %s", _exc)
                    # Overlay the colored fragment / att-site annotations from
                    # plan.product_annotations — these win over the gray base
                    # when the frontend de-duplicates by span.
                    for a in (plan.product_annotations or []):
                        ac = dict(a)
                        ac.setdefault("metadata", {})
                        if isinstance(ac.get("metadata"), dict):
                            ac["metadata"]["viz_layer"] = "overlay_colored"
                        gateway_product_annotations.append(ac)


                    # Build visualization
                    viz = {
                        "type": "design",
                        "title": f"{modules[0]['canonical_id']}_Gateway_Product",
                        "sequence": plan.product_sequence,
                        "topology": "circular",
                        "total_length": len(plan.product_sequence) if plan.product_sequence else 0,
                        "annotations": gateway_product_annotations,
                        "method": "gateway_cloning"
                    }
                    logger.info(f"[Gateway viz] annotations count: {len(viz.get('annotations', []))}")
                    if not viz["annotations"]:
                        logger.info("[Gateway viz] WARNING: No annotations in viz!")



        elif intent == "gibson_design" and seq_data["count"] >= 2:
            # Step 3b of workflow_input migration: extracted fragments +
            # intent params are projected into a WorkflowInput first; the
            # GibsonRequest is then built from the canonical input shape.
            from .workflow_input_adapters import build_for_gibson_design
            workflow_input = build_for_gibson_design(
                seq_data=seq_data,
                intent_result=intent_result,
                canonical_request=canonical_request,
                session_id=session_id,
            )
            request = _build_gibson_request_from_workflow_input(
                workflow_input, include_explanation,
            )
            result = design_gibson_primers(request)

            n_frags = len(request.fragments)
            n_primers = sum(
                1 for p in result.get("primers_by_fragment", []) if p.get("needs_primers")
            )
            reply = (
                f"Gibson assembly designed for {n_frags} fragment(s) — "
                f"{n_primers} fragment(s) require new primers."
            )
            viz = build_gibson_viz(result)
            files = build_gibson_files(result)

        # ── Plasmid design from natural language description ──────────────────
        elif intent == "plasmid_design":
            # Step 3f of workflow_input migration: PlasmidSpec is built by
            # the adapter and stashed on workflow_input.provenance.
            from .describe_plasmid_handler import describe_plasmid
            from .workflow_input_adapters import build_for_plasmid_design
            workflow_input = build_for_plasmid_design(
                message=message,
                intent_result=intent_result,
                canonical_request=canonical_request,
                session_id=session_id,
            )
            design_result = await describe_plasmid(workflow_input)
            reply = design_result["reply"]
            viz = design_result.get("viz")
            files = design_result.get("files") or []

        # ── Site-directed mutagenesis ─────────────────────────────────────────
        # SDM can use either target file OR inventory file (first one) as the plasmid to modify
        elif intent == "sdm_design" and (has_target or has_inventory):
            # Step 3g of workflow_input migration: file resolution + GB parse
            # + features extraction live in build_for_sdm_design. The handler
            # keeps the SDMMutationParser logic, conditional annotation
            # augmentation, and SDMOperator call unchanged — they read
            # plasmid_seq, features, file_text, record (locus name) from
            # workflow_input.target.
            from .workflow_input_adapters import build_for_sdm_design
            workflow_input = await build_for_sdm_design(
                file=file,
                inventory_files=inventory_files,
                has_target=has_target,
                has_inventory=has_inventory,
                intent_result=intent_result,
                canonical_request=canonical_request,
                session_id=session_id,
            )

            if workflow_input.target is None:
                return {
                    "ok": False,
                    "reply": "Please upload a plasmid file (GenBank format) to modify.",
                    "sessionId": session_id,
                }
            _parse_err = workflow_input.target.metadata.get("parse_error")
            if _parse_err:
                return {
                    "ok": False,
                    "reply": f"Failed to parse GenBank file: {_parse_err}",
                    "sessionId": session_id,
                }

            plasmid_seq = workflow_input.target.sequence
            features = list(workflow_input.target.annotations or [])
            file_text = workflow_input.target.gb_text or ""
            _record_name = workflow_input.target.metadata.get("record_name") or workflow_input.target.name or "plasmid"
            logger.info("SDM: Using %s (%d bp, %d features)",
                        workflow_input.target.source_file or _record_name,
                        len(plasmid_seq), len(features))

            # Parse mutation request
            sdm_parser = SDMMutationParser()
            sdm_params = workflow_input.workflow_args.get("sdm_params") or {}

            # Always run the full annotation pipeline when the user is
            # targeting a named feature. Many SDM targets (His-tag, FLAG, HA,
            # V5, NLS, P2A/T2A, kozak, ...) are motifs or CDS submodules that
            # aren't in the plasmid's GenBank features list — they only show
            # up after the annotation pipeline runs. Cached by sequence hash
            # so this is cheap on repeat calls.
            needs_feature_lookup = bool(
                sdm_params.get("target_feature_name")
                or sdm_params.get("terminus")
            )
            if needs_feature_lookup:
                logger.info(
                    "[SDM] Augmenting features via annotation pipeline (existing %d, "
                    "targeting %r)", len(features), sdm_params.get("target_feature_name"),
                )
                try:
                    from .annotation_cache import annotate_cached
                    ann_full = await annotate_cached(
                        plasmid_seq, circular=True, depth="full",
                    )
                    # Merge module annotations (Pol II cassettes, lentiviral
                    # modules, NLS / tag / linker / 2A submodules, etc.)
                    for m in (ann_full.get("module_annotations")
                              or ann_full.get("modules") or []):
                        features.append({
                            "type": m.get("type") or m.get("module_type") or "misc_feature",
                            "start": int(m.get("start", 0)),
                            "end": int(m.get("end", 0)),
                            "strand": m.get("strand", 1),
                            "name": (m.get("name") or m.get("module_type")
                                     or m.get("label") or ""),
                        })
                    # Merge motifs (Kozak, His-tag, FLAG, HA, V5, NLS, ...)
                    for mo in (ann_full.get("motifs") or []):
                        features.append({
                            "type": "motif",
                            "start": int(mo.get("start", 0)),
                            "end": int(mo.get("end", 0)),
                            "strand": mo.get("direction", 1),
                            "name": mo.get("name") or mo.get("motif_type") or "",
                        })
                    # Merge flat features (CDS, promoters, etc.) so single-name
                    # KB hits ("eGFP", "Cas9") still resolve.
                    for f in (ann_full.get("annotations") or []):
                        features.append({
                            "type": f.get("feat_type") or f.get("type") or "misc_feature",
                            "start": int(f.get("start", 0)),
                            "end": int(f.get("end", 0)),
                            "strand": f.get("direction", 1),
                            "name": f.get("name") or f.get("feature_name") or "",
                        })
                    logger.info("[SDM] Feature count after annotation: %d", len(features))
                except Exception as e:
                    logger.warning(f"[SDM] Pre-emptive annotation failed: {e}")

            mutation_spec = await sdm_parser.parse_mutation_request(
                message=message,
                plasmid_sequence=plasmid_seq,
                plasmid_features=features,
                sdm_params=sdm_params,
            )

            # Fallback annotation if first parse couldn't find the target
            if "NEEDS_ANNOTATION_PIPELINE" in (mutation_spec.warnings or []):
                logger.info("[SDM] Retry annotation (parser requested it)")
                try:
                    ann_req = AnnotateRequest(gb_text=file_text, session_id=session_id)
                    ann_result = await annotate_genbank(ann_req)
                    features.extend(_features_from_annotate_response(ann_result))
                except Exception as e:
                    logger.warning(f"[SDM] Fallback annotation failed: {e}")

                mutation_spec = await sdm_parser.parse_mutation_request(
                    message=message,
                    plasmid_sequence=plasmid_seq,
                    plasmid_features=features,
                    sdm_params=sdm_params,
                )
            
            if mutation_spec.confidence < 0.5:
                reply = (
                    f"I couldn't confidently identify the mutation target. "
                    f"Reason: {mutation_spec.reasoning}\n\n"
                    f"Please specify more precisely, for example:\n"
                    f"- 'Delete the His-tag'\n"
                    f"- 'Insert FLAG tag at position 500'\n"
                    f"- 'Change codon 45 from Arg to Ala'\n"
                    f"- 'Substitute AATTCC with GGCCAA'"
                )
            else:
                # Design SDM primers
                sdm_op = SDMOperator()
                plan = sdm_op.evaluate(
                    template_seq=plasmid_seq,
                    old_seq=mutation_spec.old_sequence,
                    new_seq=mutation_spec.new_sequence,
                    template_name=_record_name,
                    insertion_position=mutation_spec.target_start if mutation_spec.mutation_type == "insertion" else None,
                )
                
                if not plan.feasible:
                    reasons = "; ".join(plan.infeasibility_reasons)
                    reply = f"SDM is not feasible for this edit: {reasons}"
                else:
                    result = plan.to_dict()
                    viz = build_sdm_viz(result)
                    files = build_sdm_files(result)
                    reply = plan.summary
                    
                    if plan.warnings:
                        reply += "\n\nWarnings:\n" + "\n".join(f"- {w}" for w in plan.warnings)
                    
                    if mutation_spec.feature_context:
                        reply += f"\n\nTargeting feature: {mutation_spec.feature_context}"

        # ── sgRNA Golden Gate oligo design ─────────────────────────────────────
        elif intent == "sgrna_golden_gate":
            # Step 3e of workflow_input migration: gRNA + name + vector +
            # enzyme resolution lives in build_for_sgrna_golden_gate. Vector
            # features (parsed by parse_genbank_features) ride on
            # vector.metadata["features"]. Validation of the gRNA (length /
            # alphabet) stays here so the user-facing error messages remain
            # exactly the same.
            from .workflow_input_adapters import build_for_sgrna_golden_gate
            workflow_input = await build_for_sgrna_golden_gate(
                message=message,
                file=file,
                inventory_files=inventory_files,
                has_target=has_target,
                has_inventory=has_inventory,
                seq_data=seq_data,
                intent_result=intent_result,
                canonical_request=canonical_request,
                session_id=session_id,
            )
            grna_seq = workflow_input.workflow_args.get("grna_sequence")
            grna_name = workflow_input.workflow_args.get("grna_name") or "sgRNA"

            if not grna_seq:
                reply = (
                    "I can help you design oligos for sgRNA Golden Gate cloning!\n\n"
                    "Please provide a gRNA target sequence (17-25 bp, without PAM). For example:\n"
                    "> *Design oligos to clone gRNA GAGTCCGAGCAGAAGAAGAA into lentiCRISPR v2*"
                )
            else:
                # Validate gRNA
                grna_seq = grna_seq.strip().upper()
                invalid_chars = set(grna_seq) - set("ACGTN")
                if invalid_chars:
                    reply = f"The gRNA sequence contains invalid characters: {invalid_chars}. Only A, C, G, T are allowed."
                elif len(grna_seq) < 17 or len(grna_seq) > 30:
                    reply = f"The gRNA sequence should be 17-30 bp, but got {len(grna_seq)} bp."
                else:
                    # None == "predesign tries all supported Type IIs enzymes".
                    # Only locked to a specific enzyme if the user/canonical
                    # request explicitly named one.
                    enzyme = workflow_input.workflow_args.get("enzyme") or None
                    if workflow_input.vector is None:
                        logger.error("Failed to load vector for assembly (no upload + default fallback failed).")
                        reply = "Failed to load vector for assembly. Please try again."
                        vector_seq = None
                        vector_features = None
                        vector_name = "lentiCRISPR v2"
                    else:
                        vector_seq = workflow_input.vector.sequence
                        vector_features = workflow_input.vector.metadata.get("features") or []
                        vector_name = workflow_input.vector.name
                        logger.info(
                            f"sgRNA: Using vector: {vector_name} ({len(vector_seq)} bp) "
                            f"[source={workflow_input.vector.source}]"
                        )

                    if vector_seq:
                        # ---- Predesign: cassette recognition + Type IIs pair
                        # + sticky-end derivation + vector validation -------
                        from .annotation_cache import annotate_cached
                        from .cloning.sgrna_predesign import (
                            predesign_sgrna_vector,
                            validate_assembled_cassette,
                        )

                        try:
                            _vec_ann = await annotate_cached(
                                vector_seq, circular=True, depth="full",
                            )
                        except Exception as _exc:
                            logger.warning("sgRNA vector annotation failed: %s", _exc)
                            _vec_ann = {}

                        vec_annotations = (
                            _vec_ann.get("annotations")
                            or _vec_ann.get("features")
                            or []
                        )
                        vec_modules = (
                            _vec_ann.get("module_annotations")
                            or _vec_ann.get("modules")
                            or []
                        )

                        predesign = predesign_sgrna_vector(
                            vector_sequence=vector_seq,
                            annotations=vec_annotations,
                            modules=vec_modules,
                            requested_enzyme=enzyme,
                            grna_starts_with_g=grna_seq.startswith("G"),
                        )

                        # Surface predesign on workflow_input so downstream
                        # output builders / trace can read it.
                        workflow_input.workflow_args["predesign"] = predesign.to_dict()

                        if not predesign.validation.get("passed"):
                            warn_lines = "\n".join(
                                f"- {w}" for w in predesign.validation.get("warnings", [])
                            )
                            reply = (
                                "I couldn't validate the vector for sgRNA Golden Gate cloning.\n\n"
                                f"**Cassette resolution**: {predesign.cassette_kind or 'none'}\n"
                                f"**Enzyme tried**: {predesign.enzyme}\n\n"
                                "Pre-design checks:\n"
                                f"{warn_lines}\n"
                            )
                        else:
                            # Use enzyme + sticky ends from predesign — never
                            # the legacy GGTG/GTTT defaults.
                            enzyme = predesign.enzyme
                            oligo_design = design_sgrna_oligos(
                                grna_sequence=grna_seq,
                                five_prime_overhang=predesign.five_prime_overhang,
                                three_prime_overhang=predesign.three_prime_overhang,
                                enzyme=enzyme,
                                prepend_g_for_pol3=predesign.prepend_g,
                            )

                            # ---- Workflow: in-silico assembly using the
                            # predesign-resolved cut sites + sticky ends ----
                            try:
                                assembled = assemble_sgrna_plasmid(
                                    vector_sequence=vector_seq,
                                    grna_sequence=grna_seq,
                                    grna_name=grna_name,
                                    enzyme=enzyme,
                                    original_features=vector_features,
                                    upstream_site_pos=predesign.upstream_site_pos,
                                    downstream_site_pos=predesign.downstream_site_pos,
                                    five_prime_overhang=predesign.five_prime_overhang,
                                    three_prime_overhang=predesign.three_prime_overhang,
                                    cassette_kind=predesign.cassette_kind,
                                    promoter_name=predesign.promoter_name,
                                    promoter_end=predesign.promoter_end,
                                    prepend_g=predesign.prepend_g,
                                    vector_name=vector_name,
                                )

                                # Re-annotate the assembled plasmid (full
                                # pipeline) and validate that the inserted
                                # gRNA sits inside a recognised expression
                                # cassette.
                                full_annotations = list(assembled.annotations or [])
                                assembled_modules: List[Dict[str, Any]] = []
                                try:
                                    _full = await annotate_cached(
                                        assembled.sequence, circular=True, depth="full",
                                    )
                                    assembled_modules = (
                                        _full.get("module_annotations")
                                        or _full.get("modules") or []
                                    )
                                    _seen = {(a.get("name"), a.get("start"), a.get("end"))
                                             for a in full_annotations}
                                    for m in assembled_modules:
                                        key = (m.get("name") or m.get("module_type"),
                                               m.get("start"), m.get("end"))
                                        if key in _seen:
                                            continue
                                        full_annotations.append({
                                            "name": m.get("name") or m.get("module_type") or "feature",
                                            "start": int(m.get("start", 0)),
                                            "end": int(m.get("end", 0)),
                                            "direction": -1 if m.get("strand") == -1 else 1,
                                            "color": m.get("color") or "#6B7280",
                                        })
                                except Exception as _exc:
                                    logger.warning("sgRNA full re-annotation failed: %s", _exc)

                                cassette_check = validate_assembled_cassette(
                                    assembled_modules=assembled_modules,
                                    insert_start=assembled.insert_start,
                                    insert_end=assembled.insert_end,
                                    cassette_kind=predesign.cassette_kind,
                                )
                                assembled.cassette_validation = cassette_check

                                viz = {
                                    "type": "design",
                                    "title": f"{vector_name} + {grna_name} sgRNA (Assembled)",
                                    "sequence": assembled.sequence,
                                    "topology": "circular",
                                    "total_length": assembled.total_length,
                                    "annotations": full_annotations,
                                    "restriction_sites": assembled.restriction_sites,
                                    "ligation_junctions": assembled.ligation_junctions,
                                }

                                # Compute the transcribed spacer length: the
                                # leading Gs at the END of the upstream
                                # top-strand sticky-end window get carried
                                # into the U6 transcript before the inserted
                                # spacer (e.g. CCGG → 2 leading Gs). Plus the
                                # insert (gRNA, with G prepended only when
                                # required by adjust_grna_for_pol3).
                                trailing_g = 0
                                for _ch in reversed(predesign.five_prime_overhang or ""):
                                    if _ch == "G":
                                        trailing_g += 1
                                    else:
                                        break
                                transcribed_spacer = trailing_g + len(oligo_design.effective_grna)
                                if oligo_design.grna_was_modified:
                                    spacer_breakdown = (
                                        f"{trailing_g} G(s) from vector + 1 prepended G + "
                                        f"{len(oligo_design.grna_sequence)} nt gRNA"
                                    )
                                elif trailing_g:
                                    spacer_breakdown = (
                                        f"{trailing_g} G(s) from vector + "
                                        f"{len(oligo_design.grna_sequence)} nt gRNA — "
                                        f"vector supplies the initiating G"
                                    )
                                else:
                                    spacer_breakdown = (
                                        f"{len(oligo_design.grna_sequence)} nt gRNA used verbatim "
                                        f"(starts with G)"
                                    )

                                grna_mod_line = ""
                                if oligo_design.grna_was_modified:
                                    grna_mod_line = (
                                        f"\n**gRNA modified for Pol III**: "
                                        f"{oligo_design.grna_modification_note}"
                                    )
                                # Mirror modification on workflow_args so the
                                # workflow trace renders it.
                                workflow_input.workflow_args["grna_modification"] = {
                                    "modified": oligo_design.grna_was_modified,
                                    "note": oligo_design.grna_modification_note,
                                    "input_grna": oligo_design.grna_sequence,
                                    "effective_grna": oligo_design.effective_grna,
                                    "transcribed_spacer_length": transcribed_spacer,
                                    "vector_leading_g_count": trailing_g,
                                }
                                assembly_block_ok = True
                            except Exception as e:
                                logger.warning(f"Assembly failed, showing oligos only: {e}")
                                assembly_block_ok = False
                                grna_mod_line = ""
                                transcribed_spacer = 0
                                spacer_breakdown = ""
                                # viz stays None

                            if assembly_block_ok:
                                # Header — strict structured output the user
                                # specified: target gRNA, type IIs sites,
                                # vector, filler bp, insert bp, final size,
                                # cassette validation, transcribed spacer
                                # length. No other sections in the reply
                                # except annealed-insert; the ### Files block
                                # is appended downstream by
                                # output_builders.build_sgrna_output. Protocol
                                # details live in protocol.csv.
                                use_neb_10beta = assembled.total_length > 10000
                                upstream_pos = predesign.upstream_site_pos
                                downstream_pos = predesign.downstream_site_pos
                                up_orient = predesign.upstream_orientation
                                down_orient = predesign.downstream_orientation
                                cassette_msg = cassette_check["message"]
                                cassette_pass = cassette_check["passed"]

                                advisories = (predesign.validation or {}).get("advisories", []) or []
                                advisory_block = ""
                                if advisories:
                                    advisory_block = (
                                        "\n**⚠ Off-target risk advisory**:\n"
                                        + "\n".join(f"- {a}" for a in advisories)
                                        + "\n"
                                    )

                                reply = f"""## sgRNA Golden Gate Oligo Design

**Target gRNA**: {oligo_design.grna_sequence} ({grna_name})
**Type IIs sites**: {enzyme} at {upstream_pos} ({up_orient}) / {downstream_pos} ({down_orient})
**Vector**: {vector_name}
**Filler removed**: {assembled.filler_removed_bp} bp
**gRNA insert**: {assembled.insert_length} bp{grna_mod_line}
**Transcribed spacer**: {transcribed_spacer} nt ({spacer_breakdown})
**Final plasmid size**: {assembled.total_length:,} bp
**Cassette validation**: {'PASSED — ' if cassette_pass else 'FAILED — '}{cassette_msg}
{advisory_block}
### Annealed Insert Structure

```
{oligo_design.annealed_product_display}
```
"""

                                if oligo_design.warnings:
                                    reply += "\n**Warnings**: " + "; ".join(oligo_design.warnings)

                                # Build downloadable files. The oligo CSV stays
                                # as a legacy file for back-compat; its rows
                                # are merged into parts_order.csv by the
                                # output_builders pipeline. Protocol details
                                # are rendered into protocol.csv inside
                                # build_sgrna_output (no protocol.txt needed).
                                oligo_csv = f"""Name,Sequence,Length,Tm
{grna_name}_sgRNA_Fwd,{oligo_design.forward_oligo},{len(oligo_design.forward_oligo)},{oligo_design.forward_tm}
{grna_name}_sgRNA_Rev,{oligo_design.reverse_oligo},{len(oligo_design.reverse_oligo)},{oligo_design.reverse_tm}
"""
                                files = [
                                    {
                                        "fileName": f"{grna_name}_sgrna_oligos.csv",
                                        "mimeType": "text/csv",
                                        "dataBase64": base64.b64encode(oligo_csv.encode("utf-8")).decode("ascii"),
                                    },
                                ]

                                # Add assembled GenBank file if assembly succeeded
                                if viz:
                                    assembled_gb = _build_sgrna_genbank(
                                        sequence=assembled.sequence,
                                        annotations=assembled.annotations,
                                        grna_name=grna_name,
                                        grna_sequence=grna_seq,
                                        vector_name=vector_name,
                                    )
                                    files.append({
                                        "fileName": f"{grna_name}_assembled.gb",
                                        "mimeType": "application/octet-stream",
                                        "dataBase64": base64.b64encode(assembled_gb.encode("utf-8")).decode("ascii"),
                                    })

        # ── Golden Gate Primer Design ─────────────────────────────────────────
        elif intent == "golden_gate_primer_design":
            # Step 3c of workflow_input migration: fragment resolution + KB
            # lookup are encapsulated in build_for_golden_gate_primer_design.
            # The handler reads workflow_input.parts for fragments, workflow_args
            # for tuning knobs, and provenance for the KB audit report.
            from .workflow_input_adapters import build_for_golden_gate_primer_design
            workflow_input = build_for_golden_gate_primer_design(
                message=message,
                seq_data=seq_data,
                intent_result=intent_result,
                canonical_request=canonical_request,
                session_id=session_id,
            )
            fragments = [
                {"name": p.name, "sequence": p.sequence}
                for p in workflow_input.parts
            ]
            workflow_type = workflow_input.workflow_args["workflow_type"]
            enzyme = workflow_input.workflow_args["enzyme"]
            identification_report = workflow_input.provenance.get("identification_report")
            unresolved_candidates: List[Dict[str, Any]] = (
                (identification_report or {}).get("unresolved") or []
            )
            if identification_report:
                logger.info(
                    "[golden_gate] identified %d/%d from pLannotate KB: %s",
                    len(identification_report.get("identified") or []),
                    len(identification_report.get("candidates") or []),
                    [
                        (f["query"], f["name"], f["feature_type"], f["length"])
                        for f in (identification_report.get("identified") or [])
                    ],
                )

            if workflow_type == "multi_fragment" and len(fragments) >= 2:
                try:
                    design_result = design_multi_fragment_assembly(
                        fragments=fragments,
                        topology="circular",
                        enzyme=enzyme,
                    )
                    response = build_design_response(design_result)

                    reply = response["reply"]
                    # Build identification summary — applied after optional LLM
                    # explanation so the KB audit trail is always visible.
                    if identification_report and identification_report.get("identified"):
                        ident_lines = ["Identified from pLannotate feature KB:"]
                        for f in identification_report["identified"]:
                            ident_lines.append(
                                f"  - {f['query']} -> {f['name']} ({f['feature_type']}, {f['length']} bp) [{f['feature_id']}]"
                            )
                        ident_prefix = "\n".join(ident_lines) + "\n\n"
                    viz = response.get("viz")
                    files = response.get("files", [])
                    result = response

                except Exception as e:
                    logger.error(f"Golden Gate design failed: {e}")
                    import traceback
                    traceback.print_exc()
                    reply = f"Golden Gate primer design failed: {str(e)}"

            elif workflow_type == "scarless_deletion" and len(fragments) >= 1:
                # For deletion: need template sequence and deletion coordinates
                reply = (
                    "To design a scarless deletion, please provide:\n"
                    "- Template sequence\n"
                    "- Deletion start and end positions\n\n"
                    "Example: 'Design Golden Gate primers to delete bp 100-200 from this template: ATGC...'"
                )

            elif workflow_type == "single_fragment" and seq_data["count"] >= 2:
                # For replacement: template + insert
                reply = (
                    "To design a single fragment replacement, please provide:\n"
                    "- Template sequence\n"
                    "- Insert sequence\n"
                    "- Insertion position\n\n"
                    "Example: 'Design Golden Gate primers to insert this sequence at position 500: ...'"
                )

            else:
                ident_summary = ""
                if identification_report:
                    ident_lines = []
                    ident = identification_report.get("identified", [])
                    unresolved = identification_report.get("unresolved", [])
                    if ident:
                        ident_lines.append("Identified from pLannotate feature KB:")
                        for f in ident:
                            ident_lines.append(
                                f"  - {f['query']} -> {f['name']} ({f['feature_type']}, {f['length']} bp) [{f['feature_id']}]"
                            )
                    if unresolved:
                        ident_lines.append("")
                        ident_lines.append("Could not identify:")
                        for c in unresolved:
                            t = c.get("feature_type") or "any"
                            ident_lines.append(f"  - {c['name']} (looked up as type={t})")
                    if ident_lines:
                        ident_summary = "\n".join(ident_lines) + "\n\n"

                if len(fragments) == 0:
                    reply = (
                        ident_summary +
                        "I couldn't resolve any of the parts you mentioned into KB features with sequences.\n\n"
                        "You can either:\n"
                        "1. Provide DNA sequences directly:\n"
                        "   'Design Golden Gate primers for: Fragment1: ATGC..., Fragment2: GCTA...'\n\n"
                        "2. Use feature names with role suffixes so I can type-filter the lookup:\n"
                        "   - Promoters: CMV promoter, EF1a promoter, SV40 promoter\n"
                        "   - CDS: eGFP, mCherry, Cas9\n"
                        "   - PolyA: bGH polyA, SV40 polyA\n\n"
                        "Example: 'Design Golden Gate primers to assemble CMV promoter + eGFP + bGH polyA'"
                    )
                elif len(fragments) == 1:
                    reply = (
                        ident_summary +
                        f"Only 1 fragment resolved ({fragments[0]['name']}), but Golden Gate assembly "
                        "requires at least 2 fragments.\n\nPlease specify additional fragments or sequences."
                    )
                else:
                    reply = (
                        "For Golden Gate multi-fragment assembly, please provide at least 2 DNA sequences or feature names."
                    )

        # ── Plasmid annotation ────────────────────────────────────────────────
        elif intent == "annotate_gb" and has_target:
            # Route the chat-integrated annotate_gb demo through the full hybrid
            # pipeline (Steps 1 + 2 + 2.5 + 2.6 + 2.75) so the output matches the
            # standalone /plannotate/annotate_sequence_llm endpoint.
            #
            # Step 3a of workflow_input migration: inputs come from a
            # WorkflowInput built by `build_for_annotate_gb`; the handler reads
            # win.target instead of re-parsing the raw upload. Behavior parity
            # with the legacy block is preserved (same defaults, same error
            # message format).
            from .workflow_input_adapters import build_for_annotate_gb
            workflow_input = await build_for_annotate_gb(
                file=file,
                canonical_request=canonical_request,
                session_id=session_id,
            )
            sequence = workflow_input.target.sequence
            is_circular = workflow_input.target.topology == "circular"
            plasmid_name = workflow_input.target.name
            gb_text = workflow_input.target.gb_text or ""
            _parse_err = workflow_input.target.metadata.get("parse_error")
            if _parse_err:
                reply = f"Failed to parse GenBank file: {_parse_err}"

            if sequence:
                # Single-pass annotation via the LLM-pipeline cache. Subsequent
                # requests against the same sequence (e.g. demo re-runs)
                # short-circuit on the in-process cache instead of re-running
                # the ~4 s pipeline.
                from .annotation_cache import annotate_llm_cached
                logger.info(
                    "[annotate_gb] running annotation pipeline (cache-keyed by "
                    "sequence hash, %d bp, circular=%s)",
                    len(sequence), is_circular,
                )
                ann_result = await annotate_llm_cached(
                    sequence, circular=is_circular,
                )

                if not ann_result.get("ok"):
                    reply = "Annotation failed: " + str(ann_result.get("error", "unknown error"))
                else:
                    flat = ann_result.get("annotations") or []
                    hier = ann_result.get("hierarchical_annotations") or []
                    # Merge both: flat pLannotate features (layer=feature) + hierarchical
                    # modules/CDS submodules/cloning features. Without both, the demo
                    # was missing all pLannotate feature rows (bla, CMV, promoters, etc.).
                    raw_anns = list(flat) + list(hier)
                    normalized = []
                    for a in raw_anns:
                        start = a.get("start")
                        end = a.get("end")
                        if start is None or end is None:
                            continue
                        entry = {
                            "name": a.get("name") or a.get("label") or "feature",
                            "start": int(start),
                            "end": int(end),
                            "direction": a.get("direction", a.get("strand", 1)),
                            "color": a.get("color", "#7C3AED"),
                            "layer": a.get("layer", "feature"),
                        }
                        if a.get("module_type"):
                            entry["module_type"] = a["module_type"]
                        normalized.append(entry)

                    summary = ann_result.get("summary", {})
                    feature_count = summary.get("plannotate_feature_count", len(flat))
                    # Aggregate "module" count for the user-facing summary:
                    # rule-based modules + CDS submodules (NLS, P2A tags,
                    # protein domains) + mammalian Pol II cassettes. The
                    # standalone /plannotate/annotate_sequence_llm summary
                    # reports these separately; aggregating gives a more
                    # accurate view of "module-level structure detected".
                    module_count = (
                        summary.get("module_count", 0)
                        + summary.get("cds_submodule_count", 0)
                        + summary.get("mammalian_pol2_count", 0)
                    )
                    interaction_count = summary.get(
                        "interaction_count",
                        len(ann_result.get("interactions") or []),
                    )
                    cloning_count = summary.get("cloning_feature_count", 0)
                    reply = (
                        f"Plasmid annotated — {feature_count} feature(s), "
                        f"{module_count} module(s), "
                        f"{interaction_count} interaction(s), "
                        f"{cloning_count} cloning feature(s) "
                        f"({len(normalized)} total annotations)."
                    )
                    viz = {
                        "type": "annotation",
                        "sequence": sequence,
                        "annotations": normalized,
                        "circular": is_circular,
                        "title": plasmid_name,
                        # Surface interactions so the frontend renders
                        # module-module edges (matches standalone
                        # /plannotate/annotate_sequence_llm output).
                        "interactions": ann_result.get("interactions") or [],
                        "cloning_features": ann_result.get("cloning_features") or {},
                        "hierarchical_annotations": hier,
                    }

        # ── Restriction Cloning ───────────────────────────────────────────────
        elif intent == "restriction_cloning":
            # Step 3d of workflow_input migration: vector file read, GB
            # parse, insert resolution + KB lookup, and enzyme-override
            # extraction are all encapsulated in
            # build_for_restriction_cloning. The handler reads
            # workflow_input.vector for the vector + gb_text,
            # parts_by_role("insert") for the insert sequence, and
            # provenance for the KB record (used by the synthesis-first /
            # audit blocks below).
            from .workflow_input_adapters import build_for_restriction_cloning
            workflow_input = await build_for_restriction_cloning(
                file=file,
                inventory_files=inventory_files,
                has_target=has_target,
                has_inventory=has_inventory,
                seq_data=seq_data,
                intent_result=intent_result,
                canonical_request=canonical_request,
                session_id=session_id,
            )

            insert_name = workflow_input.workflow_args.get("insert_name")
            vector_name = workflow_input.workflow_args.get("vector_name") or ""
            enzyme_override = workflow_input.workflow_args.get("enzyme_override")
            identified_record = workflow_input.provenance.get("identified_record")
            insert_parts = workflow_input.parts_by_role("insert")
            insert_seq = insert_parts[0].sequence if insert_parts else None
            resolved_label = insert_parts[0].name if insert_parts else None

            if workflow_input.vector is None:
                reply = (
                    "Restriction cloning needs a vector GenBank file (upload it as the target "
                    "or inventory). Then describe the insert, e.g. \"clone eGFP into pUC19 with HindIII/EcoRI\"."
                )
            else:
                vector_gb_text = workflow_input.vector.gb_text or ""

                if not insert_seq:
                    reply = (
                        (f"I couldn't find a sequence for insert \"{insert_name}\" in the pLannotate KB.\n\n"
                         if insert_name else
                         "I couldn't determine the insert sequence from your prompt.\n\n")
                        + "You can:\n"
                        "- Name a KB feature, e.g. \"clone eGFP into pUC19\"\n"
                        "- Paste a literal DNA sequence\n"
                        "- Upload a source plasmid and specify the feature to clone (future: auto-extract)"
                    )
                else:
                    design = design_restriction_cloning(
                        vector_gb_text=vector_gb_text,
                        insert_seq=insert_seq,
                        insert_name=resolved_label or (insert_name or "insert"),
                        vector_name=vector_name,
                        enzyme_override=enzyme_override,
                    )
                    if design.get("error"):
                        reply = f"Restriction cloning not feasible: {design['error']}"
                    else:
                        # [CANONICAL_PROVENANCE_PLUMB] — propagate the KB
                        # provenance tag into design.metadata so the file
                        # builder can label the synthesis fragment row.
                        if identified_record and isinstance(design.get("metadata"), dict):
                            design["metadata"]["sequence_provenance"] = identified_record.get("sequence_provenance")
                            design["metadata"]["sequence_organism"] = identified_record.get("sequence_organism")
                        viz = design.get("viz")
                        result = design
                        files = build_restriction_files(design)
                        reply = design.get("reply", "")
                        if identified_record:
                            # [CANONICAL_PROVENANCE_NOTE]
                            prov = identified_record.get("sequence_provenance") or ""
                            prov_line = ""
                            if prov.startswith("backtranslated_"):
                                org = prov.replace("backtranslated_", "")
                                prov_line = (
                                    f"\nNote: DNA for {identified_record['name']} "
                                    f"was back-translated from the protein KB using "
                                    f"{org} codon preferences. This is synthetic DNA, "
                                    f"not identical to any physical reference strain.\n"
                                )
                                # Synthesis-first: the insert has no physical template,
                                # so the PCR primer path is moot. Build the gBlock ourselves
                                # and prepend a recommended-path block to the reply.
                                try:
                                    from .synthetic_fragment_builder import build_restriction_synthesis_fragment
                                    _meta = (design or {}).get("metadata") or {}
                                    _frag = build_restriction_synthesis_fragment(
                                        insert_seq=_meta.get("insert_sequence", "") or insert_seq or "",
                                        left_enzyme=_meta.get("left_enzyme", "") or "",
                                        right_enzyme=_meta.get("right_enzyme", "") or "",
                                    )
                                    if _frag.get("sequence"):
                                        _syn_block = (
                                            "**Synthesis-first path (recommended — no physical template exists):**\n"
                                            f"Order this linear fragment ({_frag['length']} bp) as a gBlock / eBlock, "
                                            f"then digest with {_meta.get('left_enzyme','')} + {_meta.get('right_enzyme','')} "
                                            f"alongside the vector and ligate — **no PCR required**.\n\n"
                                            f"```\n{_frag['sequence']}\n```\n\n"
                                            "_The PCR primer design below is kept as a fallback in case a "
                                            "physical template for the insert is acquired later._\n\n"
                                        )
                                        design["reply"] = _syn_block + (design.get("reply") or "")
                                        reply = design["reply"]
                                except Exception as _synth_exc:
                                    logger.warning(f"[restriction] synthesis-first note skipped: {_synth_exc}")
                            elif prov == "curated":
                                prov_line = (
                                    f"\nDNA source: curated reference for "
                                    f"{identified_record['name']}.\n"
                                )
                            reply = (
                                f"Identified insert from pLannotate feature KB: "
                                f"{identified_record['query']} -> {identified_record['name']} "
                                f"({identified_record['feature_type']}, {identified_record['length']} bp) "
                                f"[{identified_record['feature_id']}]"
                                f"{prov_line}\n"
                            ) + reply

        # ── Unknown / fallback ────────────────────────────────────────────────
        else:
            reply = _HELP_TEXT

        # ── Step 7: Optional AI explanation ──────────────────────────────────
        # Skip explanation for workflows that build their own structured responses.
        if include_explanation and result and intent not in (
            "unknown", "", "gateway_cloning", "restriction_cloning",
        ):
            # Build a predesign_context with the router audit when available so
            # the LLM can cite the alternative workflows evaluated + why each
            # was rejected. Supported by gibson_design and inv_gib prompt
            # builders today; others will just ignore the extra kwarg.
            _step7_predesign = None
            if _routing_context_for_llm:
                _step7_predesign = {
                    "routing_audit": _routing_context_for_llm,
                    "num_compatible_workflows": 1 + sum(
                        1 for a in _routing_context_for_llm.get("alternatives", [])
                        if a.get("feasible")
                    ),
                }
            try:
                explanation = await generate_explanation(
                    intent, result, message, predesign_context=_step7_predesign,
                )
            except TypeError:
                # Older signature fallback (positional only).
                explanation = await generate_explanation(intent, result, message)
            if explanation:
                reply = explanation

        if ident_prefix:
            reply = ident_prefix + reply

        if _routing_audit_md:
            reply = _routing_audit_md + "\n\n---\n\n" + (reply or "")

        # ────────────────────────────────────────────────────────────────
        # Output normalization (Step 4 of the cloning_workflows refactor):
        # for the 6 cloning workflows, project handler intermediate state
        # into the universal {assembled.gb, parts_order.csv, protocol.csv,
        # workflow_trace.txt} shape via output_builders + normalize_response.
        # annotate_gb and plasmid_design intentionally pass through unchanged.
        # ────────────────────────────────────────────────────────────────
        _normalized_response = None
        _NORMALIZED_INTENTS = {
            "gibson_design", "golden_gate_primer_design", "restriction_cloning",
            "sgrna_golden_gate", "gateway_cloning", "sdm_design",
        }
        _local = locals()
        _wf_in = _local.get("workflow_input")
        if intent in _NORMALIZED_INTENTS and _wf_in is not None:
            # Predesign gate: evaluate the named workflow's compatibility
            # and stash on the workflow_input provenance so the trace
            # records why this workflow was considered feasible (or not).
            try:
                _gate = await _predesign_gate(
                    intent=intent, workflow_input=_wf_in, message=message,
                )
                if _gate:
                    _wf_in.provenance["predesign_evaluation"] = _gate
                    if not _gate["compatible"]:
                        # Hard error: user explicitly named a workflow that
                        # the predesign deems infeasible.
                        return {
                            "ok": False,
                            "intent": intent,
                            "sessionId": session_id,
                            "reply": (
                                f"Predesign evaluation rejected `{intent}` for the inputs you provided.\n\n"
                                + "Reasons:\n" + "\n".join(f"- {r}" for r in _gate["reasons"])
                                + "\n\nPick a different workflow or adjust the inputs."
                            ),
                            "viz": None, "viz_list": None, "files": None,
                        }
            except Exception as _gex:
                logger.warning("[predesign gate] outer failure: %s", _gex)
            try:
                from .output_builders import (
                    build_gibson_output, build_golden_gate_output,
                    build_restriction_output, build_sgrna_output,
                    build_gateway_output, build_sdm_output,
                )
                from .output_normalizer import normalize_response
                _wf_out = None
                if intent == "gibson_design" and result is not None:
                    _wf_out = build_gibson_output(
                        workflow_input=_wf_in, result=result, llm_summary=reply,
                    )
                elif intent == "golden_gate_primer_design":
                    _design = _local.get("design_result")
                    _resp = _local.get("response")
                    if _design is not None and _resp is not None:
                        _wf_out = build_golden_gate_output(
                            workflow_input=_wf_in,
                            response=_resp, design=_design,
                            identification_report=_local.get("identification_report"),
                            llm_summary=reply,
                        )
                elif intent == "restriction_cloning":
                    _design = _local.get("design")
                    if _design is not None and isinstance(_design, dict):
                        _wf_out = build_restriction_output(
                            workflow_input=_wf_in, design=_design,
                            identified_record=_local.get("identified_record"),
                            llm_summary=reply,
                        )
                elif intent == "sgrna_golden_gate":
                    _od = _local.get("oligo_design")
                    if _od is not None:
                        _wf_out = build_sgrna_output(
                            workflow_input=_wf_in,
                            oligo_design=_od,
                            assembled=_local.get("assembled"),
                            grna_name=_local.get("grna_name", "sgRNA"),
                            grna_seq=_local.get("grna_seq", "") or "",
                            enzyme=_local.get("enzyme", "BsmBI"),
                            vector_name=_local.get("vector_name", "lentiCRISPR v2"),
                            legacy_files=files or [],
                            viz=viz, llm_summary=reply,
                        )
                elif intent == "gateway_cloning":
                    _plan = _local.get("plan")
                    _modules = _local.get("modules")
                    if _plan is not None and _modules is not None:
                        _wf_out = build_gateway_output(
                            workflow_input=_wf_in,
                            plan=_plan, modules=_modules,
                            legacy_files=files or [], viz=viz,
                            gateway_identified_inserts=_local.get("gateway_identified_inserts"),
                            llm_summary=reply,
                        )
                elif intent == "sdm_design":
                    _plan = _local.get("plan")
                    _mut = _local.get("mutation_spec")
                    if _plan is not None and _mut is not None and result is not None:
                        _wf_out = build_sdm_output(
                            workflow_input=_wf_in,
                            result=result, plan=_plan, mutation_spec=_mut,
                            viz=viz,
                            template_seq=_local.get("plasmid_seq", "") or "",
                            template_name=_local.get("_record_name", "plasmid"),
                            llm_summary=reply,
                        )
                if _wf_out is not None:
                    _normalized_response = normalize_response(_wf_out, _wf_in)
            except Exception as _norm_exc:
                logger.warning(
                    "[output normalize] %s — falling back to legacy response shape",
                    _norm_exc, exc_info=True,
                )

        if _normalized_response is not None:
            _normalized_response.setdefault("intent", intent)
            if _routing_audit_payload is not None:
                _normalized_response["routing_audit"] = _routing_audit_payload
            return _normalized_response

        response = {
            "ok": True,
            "reply": reply,
            "sessionId": session_id,
            "intent": intent,
            "viz": viz,
            "viz_list": viz_list,
            "files": files or None,
        }
        if _routing_audit_payload is not None:
            response["routing_audit"] = _routing_audit_payload
        return response

    except Exception as exc:
        import traceback
        traceback.print_exc()
        return {
            "ok": False,
            "reply": f"Error processing request: {exc}",
            "sessionId": session_id,
            "intent": intent,
            "viz": None,
            "viz_list": None,
            "files": None,
        }


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Predesign gate — evaluates the user-named workflow's compatibility
# ---------------------------------------------------------------------------
async def _predesign_gate(
    *, intent: str,
    workflow_input,
    message: str,
) -> Optional[Dict[str, Any]]:
    """Run a feasibility check for the named workflow only. Returns a
    `predesign_evaluation` dict the trace builder can render. Returns
    None on best-effort failure (predesign isn't applicable / blew up).

    The dict shape:
        {
            "method": str,
            "compatible": bool,
            "reasons": List[str],
            "evaluated": bool,
        }
    """
    try:
        from .predesign.cloning_router import CloningRouter, WorkflowMethod
        from .predesign.part_resolver import ResolvedPart
        from .predesign.target_builder import TargetPlasmid
        from .predesign.design_request import InputSource

        # Project workflow_input.parts → ResolvedPart list and the target →
        # TargetPlasmid for the compatibility check. We only need sequence
        # + length here.
        parts: List[ResolvedPart] = []
        for p in (workflow_input.parts or []):
            try:
                parts.append(ResolvedPart(
                    name=p.name, sequence=p.sequence,
                    length=p.length or len(p.sequence or ""),
                    source=InputSource.DIRECT_SEQUENCE,
                    role=p.role, origin=p.source,
                ))
            except Exception:
                continue
        target_seq = ""
        topology = "circular"
        if workflow_input.target is not None:
            target_seq = workflow_input.target.sequence
            topology = workflow_input.target.topology or "circular"
        elif workflow_input.vector is not None:
            target_seq = workflow_input.vector.sequence
            topology = workflow_input.vector.topology or "circular"
        elif parts:
            target_seq = "".join(p.sequence for p in parts)

        if not target_seq and not parts:
            return None
        target = TargetPlasmid(
            sequence=target_seq or "N" * 100,
            length=len(target_seq) or 100,
            topology=topology,
            parts=parts,
        )

        # Inline per-intent compatibility check — keeps the gate independent
        # of the router's annotations-dict shape (which was designed for the
        # auto-router, not single-intent feasibility).
        reasons: List[str] = []
        n_parts = len(parts)
        part_lens = [p.length or len(p.sequence or "") for p in parts]
        target_seq_local = (target.sequence or "") if target else ""

        if intent == "gibson_design":
            if n_parts < 2:
                reasons.append(f"Gibson needs >=2 fragments (got {n_parts}).")
            if n_parts > 6:
                reasons.append(f"Gibson typically handles <=6 fragments (got {n_parts}).")
            if part_lens and min(part_lens) < 100:
                reasons.append(f"Shortest fragment {min(part_lens)} bp is below the 100 bp Gibson minimum.")
        elif intent == "golden_gate_primer_design":
            if n_parts < 2:
                reasons.append(f"Golden Gate needs >=2 fragments (got {n_parts}).")
            if n_parts > 24:
                reasons.append(f"Golden Gate >24 fragments is impractical (got {n_parts}).")
        elif intent == "restriction_cloning":
            if not workflow_input.vector and not workflow_input.target:
                reasons.append("Restriction cloning needs a vector GenBank.")
        elif intent == "sgrna_golden_gate":
            grna = (workflow_input.workflow_args or {}).get("grna_sequence") or ""
            if not grna:
                reasons.append("No gRNA sequence provided.")
            elif len(grna) < 17 or len(grna) > 30:
                reasons.append(f"gRNA length {len(grna)} outside 17-30 bp.")
        elif intent == "sdm_design":
            if workflow_input.target is None:
                reasons.append("SDM needs an uploaded target plasmid.")
        elif intent == "gateway_cloning":
            if n_parts < 1 and workflow_input.vector is None and workflow_input.target is None:
                reasons.append("Gateway needs a donor/destination + insert.")
            try:
                from .cloning.gateway_sites import scan_att_sites
                # Need at least one att-bearing module (donor / destination /
                # entry) anywhere in the inputs.
                any_atts = bool(scan_att_sites(target_seq_local, fuzzy_threshold=0)) if target_seq_local else False
                if not any_atts:
                    for pv in (workflow_input.parts or []):
                        if scan_att_sites(pv.sequence or "", fuzzy_threshold=0):
                            any_atts = True
                            break
                if not any_atts and workflow_input.vector is not None:
                    if scan_att_sites(workflow_input.vector.sequence or "", fuzzy_threshold=0):
                        any_atts = True
                if not any_atts:
                    reasons.append("No att sites detected in target / vector / parts. Gateway needs an attB/attP/attL/attR-bearing donor or destination.")
            except Exception as _exc:
                reasons.append(f"att-site scan failed: {_exc}")

        return {
            "method": intent,
            "compatible": len(reasons) == 0,
            "reasons": reasons,
            "evaluated": True,
        }
    except Exception as exc:
        logger.warning("[predesign gate] %s — %s", intent, exc)
        return None

def _build_sgrna_genbank(
    sequence: str,
    annotations: List[Dict],
    grna_name: str,
    grna_sequence: str,
    vector_name: str,
) -> str:
    """Build a GenBank format file for the assembled sgRNA plasmid."""
    title = f"{vector_name}_{grna_name}".replace(" ", "_")[:16]

    features = [
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(sequence)}",
        '                     /organism="synthetic construct"',
        '                     /mol_type="other DNA"',
    ]

    for ann in annotations:
        start = int(ann.get("start", 0)) + 1  # Convert to 1-based
        end = int(ann.get("end", 0))
        if end <= start:
            continue

        location = f"{start}..{end}"
        if int(ann.get("direction", 1)) < 0:
            location = f"complement({location})"

        feat_type = "misc_feature"
        name = str(ann.get("name", "feature"))

        # Use appropriate feature types
        if "gRNA" in name or "sgRNA" in name:
            feat_type = "misc_RNA"
        elif "junction" in name.lower():
            feat_type = "misc_feature"
        elif "promoter" in name.lower():
            feat_type = "promoter"
        elif "CDS" in name or "Cas9" in name or "PuroR" in name or "AmpR" in name or "BleoR" in name:
            feat_type = "CDS"

        features.append(f"     {feat_type.ljust(15)} {location}")
        features.append(f'                     /label="{name.replace(chr(34), chr(39))}"')

        if "junction" in name.lower():
            features.append('                     /note="Ligation junction from Golden Gate cloning"')

    # Build origin section
    origin = ["ORIGIN"]
    seq_lower = sequence.lower()
    for i in range(0, len(seq_lower), 60):
        chunk = seq_lower[i:i + 60]
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        origin.append(f"{str(i + 1).rjust(9)} {groups}")

    return "\n".join([
        f"LOCUS       {title.ljust(16)} {str(len(sequence)).rjust(6)} bp    DNA     circular SYN",
        f"DEFINITION  {vector_name} with {grna_name} gRNA cloned via BsmBI Golden Gate.",
        "ACCESSION   .",
        "VERSION     .",
        "KEYWORDS    CRISPR; sgRNA; Golden Gate.",
        "SOURCE      synthetic DNA construct",
        "  ORGANISM  synthetic DNA construct",
        f"COMMENT     gRNA target sequence: {grna_sequence}",
        "            Cloned using BsmBI Golden Gate assembly.",
        *features,
        *origin,
        "//",
    ])
