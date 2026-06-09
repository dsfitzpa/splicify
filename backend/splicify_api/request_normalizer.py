"""
Canonical request normalizer.

Takes the raw user inputs plus the existing intent-parser output and produces
a `CloningRequest`. Works in two layers:

1. **Deterministic merge** (always runs). Pulls parts / sequences / enzymes /
   mutations out of the data we already extracted: `seq_data` from
   `extractors.extract_sequences`, and `intent_result` from
   `intent.parse_intent`. Also infers `role` / `order` from trivial syntactic
   cues (the `+`-delimited part list in Golden Gate, "clone X into Y" for
   restriction, etc.).

2. **LLM enrichment** (optional, best-effort). When OPENAI_API_KEY is set,
   a single JSON call fills in the fields the deterministic layer misses:
   constraints ("avoid internal BsaI", "keep native Kozak"), multi-mutation
   lists, role refinement, and clarification need. The LLM call is structured
   and capped; failures are swallowed — the deterministic `CloningRequest` is
   always returned.

No handler is required to read from `CloningRequest` yet. This module is
additive: the Step 3.5 hook in `chat.py` builds the request, logs it, and
passes it alongside the existing `intent_result`.
"""
from __future__ import annotations

import json
import logging
import os
import re
from typing import Any, Dict, List, Optional, Tuple

from .canonical_request import (
    AssemblyPlan,
    CloningRequest,
    Constraint,
    Enzyme,
    Mutation,
    Part,
    Vector,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Enzyme classification (used by the deterministic layer)
# ---------------------------------------------------------------------------

_TYPE_II_ENZYMES = {
    "ecori", "hindiii", "bamhi", "xhoi", "nhei", "kpni", "saci", "psti",
    "sali", "agei", "mlui", "bglii", "xbai", "spei", "noti", "asci", "clai",
    "pacI", "sbfi", "fsei", "asisi", "swai", "pmei", "smai", "haeiii",
    "avrii", "ndei", "ncoi", "apai", "bstbi", "dpni", "ssti", "apali",
}
_TYPE_IIS_ENZYMES = {
    "bsai", "bsmbi", "bbsi", "sapi", "esp3i", "paqci", "bsgi", "bpii",
    "btgzi", "eco31i",
}
# Canonical spelling table — keeps whatever case the user typed mapped back
# to the community-standard spelling.
_ENZYME_CANONICAL = {
    "ecori": "EcoRI", "hindiii": "HindIII", "bamhi": "BamHI", "xhoi": "XhoI",
    "nhei": "NheI", "kpni": "KpnI", "saci": "SacI", "psti": "PstI",
    "sali": "SalI", "agei": "AgeI", "mlui": "MluI", "bglii": "BglII",
    "xbai": "XbaI", "spei": "SpeI", "noti": "NotI", "asci": "AscI",
    "clai": "ClaI", "bsai": "BsaI", "bsmbi": "BsmBI", "bbsi": "BbsI",
    "sapi": "SapI", "esp3i": "Esp3I", "paqci": "PaqCI",
    "avrii": "AvrII", "ndei": "NdeI", "ncoi": "NcoI", "apai": "ApaI",
    "pacI": "PacI", "sbfi": "SbfI", "fsei": "FseI", "asisi": "AsiSI",
    "swai": "SwaI", "pmei": "PmeI", "smai": "SmaI",
    "bstbi": "BstBI", "dpni": "DpnI", "ssti": "SstI", "apali": "ApaLI",
}


def _classify_enzyme(name: str) -> str:
    key = name.lower()
    if key in _TYPE_IIS_ENZYMES:
        return "type_iis"
    if key in _TYPE_II_ENZYMES:
        return "type_ii"
    return "unknown"


def _canonical_enzyme_name(name: str) -> str:
    return _ENZYME_CANONICAL.get(name.lower(), name)


# ---------------------------------------------------------------------------
# Deterministic parsers — these mirror what per-workflow handlers already do,
# consolidated so every workflow sees the same answer.
# ---------------------------------------------------------------------------

_CONNECTOR_SPLIT_RE = re.compile(r"\s*[+,]\s*|\s+and\s+|\s+then\s+", re.IGNORECASE)
_PREAMBLE_STRIP_RE = re.compile(
    r"^(?:please\s+)?(?:design|assemble|clone|insert|build|make|create|generate)\b"
    r"(?:\s+(?:a|an|the|primers?|oligos?|gibson|golden\s+gate|restriction|"
    r"assembly|cloning|fragments?|construct)\b)*"
    r"(?:\s+(?:for|of|to|that|which|assemble|clone|build)\b)?"
    r"(?:\s+(?:a|an|the|primers?|oligos?|gibson|golden\s+gate|restriction|"
    r"assembly|cloning|fragments?|construct|plasmid|vector)\b)*\s*",
    re.IGNORECASE,
)
_ROLE_KEYWORDS = {
    "promoter": "promoter",
    "polya": "polyA", "poly-a": "polyA", "poly(a)": "polyA",
    "terminator": "terminator",
    "utr": "utr",
    "cds": "cds", "gene": "cds", "orf": "cds", "protein": "cds",
    "tag": "tag",
    "linker": "linker",
    "mcs": "mcs",
    "origin": "origin", "ori": "origin",
    "backbone": "backbone", "vector": "backbone",
    "selection": "selection_marker", "marker": "selection_marker",
}


def _strip_preamble(text: str) -> str:
    prev = None
    out = text.strip()
    while out != prev:
        prev = out
        out = _PREAMBLE_STRIP_RE.sub("", out).strip()
    return out


def _infer_part_role(token: str) -> Tuple[str, str]:
    """Returns (clean_name, role)."""
    t = token.strip().strip(".").strip()
    lower = t.lower()
    for kw, role in _ROLE_KEYWORDS.items():
        if lower.endswith(" " + kw) or lower.endswith("-" + kw):
            name = t[: -(len(kw))].strip(" -").strip()
            return name, role
        if lower == kw:
            return t, role
    return t, "unknown"


def _parts_from_plus_list(message: str, fragment_sequences: Dict[str, str]) -> List[Part]:
    """
    Extract ordered Part entries from a "+"/","-delimited part list in the
    message. Example: "EF1a promoter + eGFP CDS + bGH polyA" → 3 parts.

    Only runs when the message has at least one "+". Sequence-labeled fragments
    (captured by `extract_sequences`) are merged in as `source="prompt_sequence"`.
    """
    parts: List[Part] = []
    # Restrict to the segment after an "assemble"/"clone"/"design" verb if
    # present, so we don't pick up incidental "+"s inside parameter notation.
    stripped = _strip_preamble(message)
    if "+" not in stripped and "," not in stripped:
        return parts
    # Bound the scope: everything up to the first period or newline.
    scope = re.split(r"[.\n]", stripped, maxsplit=1)[0]
    tokens = [t for t in _CONNECTOR_SPLIT_RE.split(scope) if t.strip()]
    if len(tokens) < 2:
        return parts
    # Reject tokens that look like parameters ("25 bp overlap") or verbs.
    junk_re = re.compile(r"\b(overlap|tm|bp|degree|degrees|primers?|assembly)\b", re.IGNORECASE)
    for idx, tok in enumerate(tokens):
        if junk_re.search(tok):
            continue
        clean, role = _infer_part_role(tok)
        if not clean or len(clean) > 80:
            continue
        # If the fragment's sequence was also pasted (e.g. "Frag1: ATGC..."),
        # keep both: the KB name for identification, the sequence for skip-KB.
        seq = None
        source = "prompt_name"
        for label, s in fragment_sequences.items():
            if label.lower() == clean.lower():
                seq = s
                source = "prompt_sequence"
                break
        parts.append(Part(name=clean, sequence=seq, role=role, order=idx, source=source))
    return parts


def _parts_from_sequences(seq_data: Dict[str, Any], already_named: set) -> List[Part]:
    """Any labeled sequence that wasn't matched by the name-list pass."""
    out: List[Part] = []
    for label, seq in seq_data.get("fragments", {}).items():
        if label.lower() in already_named:
            continue
        out.append(Part(name=label, sequence=seq, role="fragment",
                        source="prompt_sequence"))
    return out


def _enzymes_from_message_and_intent(message: str, intent_result: Dict[str, Any]) -> List[Enzyme]:
    enzymes: List[Enzyme] = []
    seen: set = set()
    msg_lower = message.lower()

    def _add(name: str, preference: str = "required") -> None:
        canonical = _canonical_enzyme_name(name)
        if canonical.lower() in seen:
            return
        seen.add(canonical.lower())
        enzymes.append(Enzyme(name=canonical, kind=_classify_enzyme(canonical),
                              preference=preference))

    # Intent parser already pulled restriction enzymes into a list.
    rc = intent_result.get("restriction_cloning") or {}
    for enz in (rc.get("enzymes") or []):
        _add(enz)
    gg = intent_result.get("golden_gate") or {}
    if gg.get("enzyme"):
        _add(gg["enzyme"])
    sgrna = intent_result.get("sgrna") or {}
    if sgrna.get("enzyme"):
        _add(sgrna["enzyme"])

    # Sweep the message for any other known enzyme names the intent parser
    # didn't surface (e.g. user named three enzymes but parser truncated to 2).
    for key, canonical in _ENZYME_CANONICAL.items():
        # Use word boundary so "psti" doesn't match inside "phosphatidylserine".
        if re.search(rf"\b{re.escape(key)}\b", msg_lower):
            _add(canonical)

    return enzymes


def _mutations_from_intent(intent_result: Dict[str, Any]) -> List[Mutation]:
    sdm = intent_result.get("sdm") or {}
    if not sdm.get("mutation_type") or not sdm.get("target_method"):
        return []
    m = Mutation(
        mutation_type=sdm["mutation_type"],
        target_method=sdm["target_method"],
        target_feature_name=sdm.get("target_feature_name"),
        target_position_start=sdm.get("target_position_start"),
        target_position_end=sdm.get("target_position_end"),
        codon_position=sdm.get("codon_position"),
        codon_from=sdm.get("codon_from"),
        codon_to=sdm.get("codon_to"),
        terminus=sdm.get("terminus"),
        old_sequence=sdm.get("old_sequence"),
        new_sequence=sdm.get("new_sequence"),
        description=sdm.get("description"),
    )
    return [m]


_INTO_VECTOR_RE = re.compile(
    r"\b(?:into|in)\s+([A-Za-z][A-Za-z0-9_.\-]{2,40})",
    re.IGNORECASE,
)


def _vector_from_intent(message: str, intent_result: Dict[str, Any],
                         has_target: bool, has_inventory: bool) -> Optional[Vector]:
    rc = intent_result.get("restriction_cloning") or {}
    sgrna = intent_result.get("sgrna") or {}
    if rc.get("vector_name"):
        return Vector(name=rc["vector_name"],
                      source="target_file" if has_target else "prompt_name")
    if sgrna.get("vector_name"):
        return Vector(name=sgrna["vector_name"],
                      source="target_file" if has_target else "prompt_name")
    # Catch "insertion of X into Y" / "clone X into Y" for Gateway & friends.
    m = _INTO_VECTOR_RE.search(message)
    if m:
        candidate = m.group(1).strip()
        # Filter obvious non-vector words.
        if candidate.lower() not in {"the", "this", "that", "a", "an", "my", "our"}:
            return Vector(name=candidate,
                          source="target_file" if has_target else "prompt_name")
    # Backbone uploaded as a file with no explicit name is still a vector.
    if has_target:
        return Vector(source="target_file")
    return None


def _constraints_from_message(message: str) -> List[Constraint]:
    """Cheap keyword sweep; the LLM pass fills in the rest."""
    out: List[Constraint] = []
    msg_lower = message.lower()
    # "avoid internal BsaI"
    for key, canonical in _ENZYME_CANONICAL.items():
        if re.search(rf"\b(?:avoid|no|without|without any)\s+(?:internal\s+)?{re.escape(key)}\b",
                     msg_lower):
            out.append(Constraint(kind="avoid_internal_site", target=canonical,
                                  polarity="negative"))
    if re.search(r"\b(keep|preserve)\s+(?:the\s+)?(?:native\s+)?kozak", msg_lower):
        out.append(Constraint(kind="preserve_feature", target="Kozak"))
    if re.search(r"\bpreserve\s+(?:the\s+)?reading\s+frame\b", msg_lower):
        out.append(Constraint(kind="preserve_reading_frame"))
    if re.search(r"\bmoclo\b", msg_lower):
        out.append(Constraint(kind="use_standard", target="moclo"))
    if re.search(r"\bloop\s+assembly\b", msg_lower):
        out.append(Constraint(kind="use_standard", target="loop"))
    return out


def _assembly_from_intent(intent: str, intent_result: Dict[str, Any]) -> AssemblyPlan:
    method_map = {
        "gibson_design": "gibson",
        "inv_gib": "gibson",
        "golden_gate_primer_design": "golden_gate",
        "sgrna_golden_gate": "golden_gate",
        "restriction_cloning": "restriction",
        "gateway_cloning": "gateway",
        "sdm_design": "sdm",
        "pcr_design": "pcr",
        "multi_pcr_design": "pcr",
        "annotate_gb": "annotate",
        "plasmid_design": "plasmid_design",
        "repp": "plasmid_design",
    }
    topology = (intent_result.get("gibson_design") or {}).get("assembly")
    return AssemblyPlan(method=method_map.get(intent), topology=topology)


# ---------------------------------------------------------------------------
# Deterministic merge
# ---------------------------------------------------------------------------

def _build_deterministic(
    message: str,
    redacted_message: str,
    intent_result: Dict[str, Any],
    seq_data: Dict[str, Any],
    has_target: bool,
    has_inventory: bool,
) -> CloningRequest:
    intent = intent_result.get("intent", "unknown")

    # Parts: try the +/comma list first (Golden Gate / Gateway module chains),
    # then fold in any labeled sequences not already covered.
    parts = _parts_from_plus_list(message, seq_data.get("fragments", {}))
    already_named = {p.name.lower() for p in parts if p.name}
    parts.extend(_parts_from_sequences(seq_data, already_named))

    # Restriction cloning: insert_name / insert_sequence should surface as a Part.
    rc = intent_result.get("restriction_cloning") or {}
    if rc.get("insert_name") or rc.get("insert_sequence"):
        name = rc.get("insert_name")
        if not (name and name.lower() in already_named):
            parts.append(Part(
                name=name,
                sequence=rc.get("insert_sequence"),
                role="insert",
                source="prompt_sequence" if rc.get("insert_sequence") else "prompt_name",
            ))

    # sgRNA guide — record as a Part with role="sgrna_guide" so handlers that
    # migrate can consume it uniformly.
    sgrna = intent_result.get("sgrna") or {}
    if sgrna.get("grna_sequence"):
        parts.append(Part(
            name=None,
            sequence=sgrna["grna_sequence"],
            role="sgrna_guide",
            source="prompt_sequence",
            metadata={"guide_length": len(sgrna["grna_sequence"])},
        ))

    vector = _vector_from_intent(message, intent_result, has_target, has_inventory)
    enzymes = _enzymes_from_message_and_intent(message, intent_result)
    mutations = _mutations_from_intent(intent_result)
    constraints = _constraints_from_message(message)
    assembly = _assembly_from_intent(intent, intent_result)

    # Pass-through pockets
    return CloningRequest(
        intent=intent,
        parts=parts,
        vector=vector,
        enzymes=enzymes,
        mutations=mutations,
        constraints=constraints,
        assembly=assembly,
        gibson_params=intent_result.get("gibson_design"),
        pcr_params=intent_result.get("pcr"),
        sdm_params=intent_result.get("sdm"),
        sgrna_params=intent_result.get("sgrna"),
        golden_gate_params=intent_result.get("golden_gate"),
        restriction_params=intent_result.get("restriction_cloning"),
        gateway_params=intent_result.get("gateway"),
        confidence=(intent_result.get("router_notes") or {}).get("confidence", 1.0),
        raw_message=message,
        redacted_message=redacted_message,
    )


# ---------------------------------------------------------------------------
# LLM enrichment (best-effort, non-blocking)
# ---------------------------------------------------------------------------

_ENRICHMENT_SYSTEM_PROMPT = """You are enriching a canonical cloning request.
You will receive:
- A redacted user message (DNA replaced with placeholders)
- The intent label chosen by the router
- A preliminary canonical request (JSON) already populated by deterministic parsers

Your job: return a JSON patch that ADDS OR CORRECTS fields the deterministic
parser would have missed. Do not rewrite fields that are already correct.
Never output DNA sequences. Output exactly one JSON object.

Focus on:
- constraints: ["avoid_internal_site", "preserve_feature", "preserve_reading_frame",
  "junction_at", "require_domestication", "use_standard", "max_fragments",
  "selection"]. Each constraint has {kind, target, polarity, notes}.
- multi-mutation requests: if the user asked for multiple mutations, return
  them as a list of mutation objects with the full Mutation schema.
- part role refinement: if a part's role is "unknown" but the message implies
  it's e.g. a promoter/cds/terminator, return a refined role.
- needs_clarification: true ONLY if the request is genuinely ambiguous in a
  way that would cause the wrong workflow to run. Include a single short
  clarification_question if so.

OUTPUT SCHEMA:
{
  "constraints": [{"kind": str, "target": str|null, "polarity": "positive"|"negative"|"preferred"|"dispreferred", "notes": str|null}],
  "mutations": [...Mutation...],
  "part_role_refinements": [{"name": str, "role": str}],
  "needs_clarification": bool,
  "clarification_question": str|null,
  "notes": [str]
}
Return {} with no keys if you have nothing to add."""


async def _enrich_with_llm(
    redacted_message: str,
    preliminary: CloningRequest,
    api_key: str,
) -> Optional[Dict[str, Any]]:
    try:
        from openai import AsyncOpenAI
    except Exception:
        return None
    client = AsyncOpenAI(api_key=api_key)
    user_content = (
        f"INTENT: {preliminary.intent}\n"
        f"MESSAGE: {redacted_message}\n"
        f"PRELIMINARY_REQUEST: {json.dumps(preliminary.to_dict(), default=str)}"
    )
    try:
        response = await client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {"role": "system", "content": _ENRICHMENT_SYSTEM_PROMPT},
                {"role": "user", "content": user_content},
            ],
            response_format={"type": "json_object"},
            max_tokens=600,
        )
        return json.loads(response.choices[0].message.content)
    except Exception as exc:
        logger.info("request_normalizer LLM enrichment skipped: %s", exc)
        return None


def _apply_enrichment(req: CloningRequest, patch: Dict[str, Any]) -> None:
    for c in patch.get("constraints") or []:
        if not isinstance(c, dict) or not c.get("kind"):
            continue
        if req.has_constraint(c["kind"], c.get("target")):
            continue
        req.constraints.append(Constraint(
            kind=c["kind"], target=c.get("target"),
            polarity=c.get("polarity", "positive"), notes=c.get("notes"),
        ))
    for m in patch.get("mutations") or []:
        if not isinstance(m, dict) or not m.get("mutation_type") or not m.get("target_method"):
            continue
        if any(
            existing.mutation_type == m["mutation_type"]
            and existing.target_method == m["target_method"]
            and existing.codon_position == m.get("codon_position")
            and existing.target_feature_name == m.get("target_feature_name")
            for existing in req.mutations
        ):
            continue
        req.mutations.append(Mutation(
            mutation_type=m["mutation_type"],
            target_method=m["target_method"],
            target_feature_name=m.get("target_feature_name"),
            target_position_start=m.get("target_position_start"),
            target_position_end=m.get("target_position_end"),
            codon_position=m.get("codon_position"),
            codon_from=m.get("codon_from"),
            codon_to=m.get("codon_to"),
            terminus=m.get("terminus"),
            old_sequence=m.get("old_sequence"),
            new_sequence=m.get("new_sequence"),
            truncation_length=m.get("truncation_length"),
            description=m.get("description"),
        ))
    for refinement in patch.get("part_role_refinements") or []:
        if not isinstance(refinement, dict):
            continue
        name = (refinement.get("name") or "").lower()
        role = refinement.get("role")
        if not name or not role:
            continue
        for p in req.parts:
            if p.name and p.name.lower() == name and p.role == "unknown":
                p.role = role
    if patch.get("needs_clarification"):
        req.needs_clarification = True
        if patch.get("clarification_question"):
            req.clarification_question = patch["clarification_question"]
    for note in patch.get("notes") or []:
        if isinstance(note, str):
            req.normalizer_notes.append(note)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

async def normalize_request(
    message: str,
    redacted_message: str,
    intent_result: Dict[str, Any],
    seq_data: Dict[str, Any],
    has_target: bool,
    has_inventory: bool,
    enable_llm_enrichment: bool = True,
) -> CloningRequest:
    """
    Build a CloningRequest from raw inputs + existing intent-parser output.

    Always returns a CloningRequest; LLM enrichment is best-effort and silently
    skipped on failure.
    """
    req = _build_deterministic(
        message=message,
        redacted_message=redacted_message,
        intent_result=intent_result,
        seq_data=seq_data,
        has_target=has_target,
        has_inventory=has_inventory,
    )

    api_key = os.getenv("OPENAI_API_KEY") if enable_llm_enrichment else None
    if api_key:
        patch = await _enrich_with_llm(redacted_message, req, api_key)
        if patch:
            _apply_enrichment(req, patch)

    return req
