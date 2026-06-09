"""Output-emitter tool schemas + dispatch chain.

Layered on top of v1's `dispatch_tool`: emitter names route to local
handlers in `agent_v2.outputs.*`; everything else falls through to v1.
The Main Agent's `tools` list = v1's AIPLASMIDDESIGN_TOOLS + EMITTER_TOOLS.
"""
from __future__ import annotations

from typing import Any, Optional

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.outputs.assembled_gb import emit_assembled_gb
from agent_v2.outputs.parts_order_csv import emit_parts_order
from agent_v2.outputs.protocol_csv import emit_protocol
from agent_v2.outputs.workflow_trace import emit_workflow_trace
from agent_v2.outputs.guides_csv import emit_guides_csv
from agent_v2.outputs.guides_gb import emit_guides_gb
from agent_v2.feature_resolver import resolve_feature_position
from agent_v2.crispr_tools import (
    design_guides_tool,
    design_pegrnas_tool,
    design_primers_tool,
    find_genomic_record_tool,
)


EMIT_ASSEMBLED_GB_TOOL: dict[str, Any] = {
    "name": "emit_assembled_gb",
    "description": (
        "Emit the final assembled plasmid as assembled.gb. Call this AFTER all "
        "assembly + verification is complete, with the attachment_id of the "
        "assembled product (typically att_product_N where N >= 2 if a single "
        "input was uploaded). Returns a file envelope and the product length."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id": {
                "type": "string",
                "description": "att_product_N id of the assembled plasmid.",
            },
        },
        "required": ["attachment_id"],
    },
}


EMIT_PARTS_ORDER_TOOL: dict[str, Any] = {
    "name": "emit_parts_order",
    "description": (
        "Emit parts_order.csv listing every part required for the assembly. "
        "Classify each into inventory (on hand) / order_addgene / order_synthesis "
        "(>=40 bp gBlock) / order_oligo (<40 bp). Include rough vendor + cost + "
        "lead-time defaults (overridable per part)."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "parts": {
                "type": "array",
                "description": (
                    "List of parts. Each: part_id, name, length_bp (required); "
                    "optional: origin (inventory|user_supplied|input|designed_oligo|knowledge_base), "
                    "addgene_id, vendor_hint, cost_usd, lead_time_days, notes."
                ),
                "items": {
                    "type": "object",
                    "properties": {
                        "part_id": {"type": "string"},
                        "name": {"type": "string"},
                        "length_bp": {"type": "integer"},
                        "origin": {"type": "string"},
                        "addgene_id": {"type": "string"},
                        "vendor_hint": {"type": "string"},
                        "cost_usd": {"type": "number"},
                        "lead_time_days": {"type": "integer"},
                        "notes": {"type": "string"},
                    },
                    "required": ["part_id", "name", "length_bp"],
                },
            }
        },
        "required": ["parts"],
    },
}


EMIT_PROTOCOL_TOOL: dict[str, Any] = {
    "name": "emit_protocol",
    "description": (
        "Emit protocol.csv: a wet-lab step-by-step playbook for the chosen "
        "assembly method. Methods supported: gibson, gateway, restriction, "
        "sdm, sgrna_gg, golden_gate. Pass `custom_steps` to override or "
        "append individual steps."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "assembly_method": {
                "type": "string",
                "description": "gibson | gateway | restriction | sdm | sgrna_gg | golden_gate | none",
            },
            "custom_steps": {
                "type": "array",
                "description": (
                    "Optional list of partial step dicts. If `step_num` is set "
                    "(1-indexed) the matching step is updated in place; "
                    "otherwise the step is appended."
                ),
                "items": {
                    "type": "object",
                    "properties": {
                        "step_num": {"type": "integer"},
                        "category": {"type": "string"},
                        "inputs": {"type": "string"},
                        "output": {"type": "string"},
                        "instrument": {"type": "string"},
                        "time_min": {"type": "integer"},
                        "temp_C": {"type": "string"},
                        "reagents": {"type": "string"},
                        "notes": {"type": "string"},
                    },
                },
            },
        },
        "required": ["assembly_method"],
    },
}


EMIT_WORKFLOW_TRACE_TOOL: dict[str, Any] = {
    "name": "emit_workflow_trace",
    "description": (
        "Emit workflow_trace.txt: a plain-text aggregate of the agent run. "
        "Sections: AGENT TRACE (tool calls), PLAN.MD (final checked plan), "
        "DECISIONS LEDGER (chosen vs alternative), optional DESIGN "
        "VERIFICATION + EXPLORE FINDINGS. Pass everything you have - all "
        "fields are optional except this is the LAST tool you call."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "session_id": {"type": "string"},
            "turn_id": {"type": "string"},
            "assembly_method": {"type": "string"},
            "product_attachment_id": {"type": "string"},
            "timestamp": {"type": "string"},
            "agent_trace": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "iteration": {"type": "integer"},
                        "tool": {"type": "string"},
                        "args_summary": {"type": "string"},
                        "result_keys": {"type": "array", "items": {"type": "string"}},
                    },
                },
            },
            "plan_md": {"type": "string"},
            "decisions": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "choice": {"type": "string"},
                        "alternative": {"type": "string"},
                        "reason": {"type": "string"},
                    },
                },
            },
            "verifier": {
                "type": "object",
                "properties": {
                    "passed": {"type": "boolean"},
                    "warnings": {"type": "array"},
                },
            },
            "findings": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "role": {"type": "string"},
                        "summary_md": {"type": "string"},
                        "key_facts": {"type": "object"},
                    },
                },
            },
        },
    },
}


RESOLVE_FEATURE_POSITION_TOOL: dict[str, Any] = {
    "name": "resolve_feature_position",
    "description": (
        "Map a user reference like 'D10A in Cas9' or '-35 box of CMV promoter' "
        "to a deterministic plasmid coordinate. Call this BEFORE simulate_assembly / "
        "score_sanger_primer / digest_plasmid whenever the user names a feature "
        "instead of giving a numeric position. Returns the bp position plus the "
        "codon and amino acid (for kind='aa_residue'). Use kind='aa_residue' for "
        "amino-acid substitutions, 'bp_offset' for upstream/downstream offsets, "
        "'feature_start' / 'feature_end' for boundary references."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id": {
                "type": "string",
                "description": "Registered attachment id (att_input_N or att_product_N).",
            },
            "feature_name": {
                "type": "string",
                "description": (
                    "Case-insensitive substring match against annotated feature names "
                    "(e.g. 'Cas9', 'GFP', 'AmpR', 'CMV promoter')."
                ),
            },
            "kind": {
                "type": "string",
                "enum": ["aa_residue", "bp_offset", "feature_start", "feature_end"],
                "description": (
                    "aa_residue: 1-indexed amino acid number within the CDS. "
                    "bp_offset: 0-indexed bp offset from feature start (or end if - strand). "
                    "feature_start / feature_end: just the boundary."
                ),
            },
            "offset": {
                "type": "integer",
                "description": (
                    "For kind=aa_residue: 1-indexed residue number (e.g. 10 for D10A). "
                    "For kind=bp_offset: bp offset (negative allowed for upstream)."
                ),
            },
        },
        "required": ["attachment_id", "feature_name", "kind"],
    },
}


DESIGN_GUIDES_TOOL: dict[str, Any] = {
    "name": "design_guides",
    "description": (
        "Design Cas9 (or Cas12a) sgRNAs against a region of a registered "
        "plasmid. PAM-scans both strands, scores with Doench 2014 Rule Set 1 "
        "(default) or a deterministic heuristic, counts on-plasmid off-targets. "
        "Use this when the user asks for an sgRNA targeting a specific base or "
        "amino-acid residue: first resolve the residue with "
        "resolve_feature_position, then call design_guides with a ~50 bp window "
        "centred on the resolved plasmid position."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id": {"type": "string", "description": "att_input_N / att_product_N id."},
            "region_start": {"type": "integer", "description": "1-indexed inclusive. **Pass a single integer**, not a coordinate dict from a prior tool call."},
            "region_end":   {"type": "integer", "description": "1-indexed inclusive. **Pass a single integer.**"},
            "pam":           {"type": "string", "description": "IUPAC PAM. NGG (Cas9), TTTV (Cas12a)."},
            "guide_length":  {"type": "integer", "description": "Default 20."},
            "pam_position":  {"type": "string", "description": "3prime (Cas9) | 5prime (Cas12a)."},
            "max_guides":    {"type": "integer", "description": "Default 5; cap top-N after sorting."},
            "min_score":     {"type": "number"},
            "score_method":  {"type": "string", "description": "doench2014 (default) | heuristic"},
        },
        "required": ["attachment_id", "region_start", "region_end"],
    },
}


DESIGN_PEGRNAS_TOOL: dict[str, Any] = {
    "name": "design_pegrnas",
    "description": (
        "Design pegRNAs for prime editing against a substitution / insertion / "
        "deletion at a 1-based plasmid position. Returns top-N pegRNAs with "
        "spacer / scaffold / RTT / PBS components and the full assembled pegRNA "
        "(5' -> 3': spacer-scaffold-RTT-PBS). For PE3, an opposite-strand "
        "nicking gRNA is selected within +/-100 nt of the pegRNA nick; PE3b is "
        "flagged when the ngRNA overlaps the edit (substitution edits only). "
        "Resolve the edit position with resolve_feature_position first if the "
        "user named a residue (e.g. 'D10A in Cas9')."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id": {"type": "string"},
            "edit_start":    {"type": "integer", "description": "1-indexed inclusive + strand position. **If you called resolve_feature_position, pass its `edit_start_1based` field directly.** Do NOT pass `plasmid_position` (which is 0-indexed) or the whole result object."},
            "edit_end":      {"type": "integer", "description": "1-indexed inclusive + strand position. **Pass `edit_end_1based` from resolve_feature_position.**"},
            "alt":           {"type": "string", "description": "+ strand replacement / inserted sequence. For + strand features, alt = desired codon. For - strand features (feature_strand='-' in resolve_feature_position output), alt = revcomp(desired_sense_codon). Required for substitution + insertion. Example: KEAP1 R15C is on the - strand; sense Cys codon TGT -> alt='ACA' (= revcomp('TGT'))."},
            "edit_type":     {"type": "string", "description": "substitution (default) | insertion | deletion"},
            "n_results":     {"type": "integer", "description": "Default 3."},
            "use_pe3":       {"type": "boolean", "description": "Default true."},
        },
        "required": ["attachment_id", "edit_start", "edit_end"],
    },
}


DESIGN_PRIMERS_TOOL: dict[str, Any] = {
    "name": "design_primers",
    "description": (
        "Primer3-driven primer design with application-aware ranking. "
        "application='fragment' (PCR cloning), 'sanger' (banded design + 10-axis "
        "Sanger quality scorer), or 'illumina' (Nextera read-1 / read-2 "
        "adapters prepended for NGS amplicon prep). For NGS / Sanger readouts "
        "around an sgRNA cut, set region = sgRNA +/-250 bp and excluded = "
        "sgRNA +/-75 bp. Region coords are 1-indexed inclusive on the "
        "PLASMID. The tool slices the template internally and translates "
        "returned primer positions back to plasmid coordinates."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id":     {"type": "string"},
            "region_start":      {"type": "integer", "description": "1-indexed inclusive on plasmid. **Pass a single integer.** Same convention as design_guides — do NOT pass the entire resolve_feature_position result dict."},
            "region_end":        {"type": "integer", "description": "1-indexed inclusive on plasmid. **Pass a single integer.**"},
            "excluded_start":    {"type": "integer", "description": "Optional. 1-indexed inclusive on plasmid. **Pass a single integer.**"},
            "excluded_end":      {"type": "integer", "description": "Optional. 1-indexed inclusive on plasmid. **Pass a single integer.**"},
            "application":       {"type": "string", "description": "fragment | sanger (default) | illumina"},
            "product_size_min":  {"type": "integer"},
            "product_size_max":  {"type": "integer"},
            "primer_opt_tm":     {"type": "number", "description": "Default 60 C."},
            "num_return":        {"type": "integer", "description": "primer3 candidate pairs. Default 5."},
        },
        "required": ["attachment_id", "region_start", "region_end"],
    },
}


EMIT_GUIDES_CSV_TOOL: dict[str, Any] = {
    "name": "emit_guides_csv",
    "description": (
        "Emit guides.csv — one wide row per designed sgRNA, pegRNA, ngRNA, "
        "Sanger / Illumina primer, and (optionally) sgRNA cloning oligo pair. "
        "Pass digested entries from design_guides / design_pegrnas / "
        "design_primers directly; the emitter does not re-run any design. "
        "Cloning oligos (for sgRNA Golden Gate into pX330 / lentiCRISPR / etc.) "
        "are passed in pre-computed — the agent picks the destination vector "
        "and synthesises the BbsI/BsmBI-overhanged top + bottom oligos."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "guides":          {"type": "array", "description": "sgRNA dicts from design_guides; each carries name, spacer, pam, start, end, direction, score, score_method, gc_fraction, n_offtargets. Add `target_attachment_id` per row."},
            "pegrnas":         {"type": "array", "description": "pegRNA dicts from design_pegrnas; each carries name, rank, spacer, pam, spacer_start, spacer_end, direction, cas9_score, rtt, pbs, scaffold, full_pegrna, predicted_efficiency, is_dpam, is_pe3b, ngrna (nested), edit_type, edit_ref, edit_alt. Add `target_attachment_id` per row."},
            "primers":         {"type": "array", "description": "Primer-pair dicts from design_primers; each carries application, left_primer / left_annealing / left_adapter, right_primer / right_annealing / right_adapter, left_pos_plasmid, right_pos_plasmid, product_size, pair_label, etc. Add `target_attachment_id` per row."},
            "cloning_oligos":  {"type": "array", "description": "Optional. List of {name, spacer, oligo_top, oligo_bottom, notes} pairs the agent synthesised for sgRNA Golden Gate cloning."},
            "descriptor":      {"type": "string", "description": "Short tag prefixed to ALL CRISPR output filenames (guides.csv, guides.gb, parts_order.csv, protocol.csv, workflow_trace.txt). Use the gene + edit identifier, e.g. \"KEAP1_R15C\" or \"CGAS_D431S_K479E\". Set ONCE on emit_guides_csv and it propagates to the others."},
        },
    },
}


EMIT_GUIDES_GB_TOOL: dict[str, Any] = {
    "name": "emit_guides_gb",
    "description": (
        "Emit guides.gb — the input attachment's GenBank file (plasmid or "
        "genomic) with the chosen sgRNAs / pegRNAs / ngRNAs / Sanger + Illumina "
        "primers appended as misc_RNA / primer_bind features. The frontend's "
        "CircularPlasmidViewer renders this via /agent_v2/annotate-on-upload, "
        "so the user sees the guide cut sites + flanking amplicons overlaid on "
        "their original sequence. Coordinates are in the target's plasmid coord "
        "space (0-indexed). Reads the source GenBank text stashed by the router "
        "at upload time; falls back to a minimal record if not cached."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "target_attachment_id": {"type": "string"},
            "guides":   {"type": "array", "description": "Same shape as emit_guides_csv.guides. **If omitted, the dispatcher reuses what was passed to the preceding emit_guides_csv call** — but pass it explicitly when in doubt."},
            "pegrnas":  {"type": "array", "description": "Same shape as emit_guides_csv.pegrnas. Same fallback rule as guides."},
            "primers":  {"type": "array", "description": "Same shape as emit_guides_csv.primers. Requires left_pos_plasmid / right_pos_plasmid for annotation placement. Same fallback rule."},
            "descriptor": {"type": "string", "description": "Optional filename prefix; inherited from emit_guides_csv if omitted."},
        },
        "required": ["target_attachment_id"],
    },
}


# emit_workflow_trace is INTENTIONALLY excluded — it's emitted server-side
# by the orchestrator right after main_agent returns, so the LLM doesn't
# have to compose the big args dict + end with a final text response in
# the same turn (a known cause of the run pushing past the 5-minute
# Vercel SSE-proxy timeout). The schema + handler remain in this file so
# tests + the dispatch table can still reference them.
EMITTER_TOOLS: list[dict[str, Any]] = [
    EMIT_ASSEMBLED_GB_TOOL,
    EMIT_PARTS_ORDER_TOOL,
    EMIT_PROTOCOL_TOOL,
    EMIT_GUIDES_CSV_TOOL,
    EMIT_GUIDES_GB_TOOL,
]

# Non-emitter local tools (helpers like coordinate resolvers).
RESOLVER_TOOLS: list[dict[str, Any]] = [RESOLVE_FEATURE_POSITION_TOOL]

FIND_GENOMIC_RECORD_TOOL: dict[str, Any] = {
    "name": "find_genomic_record",
    "description": (
        "Look up an external gene by symbol on NCBI, download its .gb file, "
        "annotate it as a genomic record, and register it in the "
        "AttachmentRegistry. Returns the new attachment_id along with the "
        "RefSeqGene / mRNA accession, organism, gene symbol, and CDS "
        "protein length. Use this when the user references a gene "
        "(e.g. 'CGAS', 'KEAP1', 'NRF2') the agent does NOT already have a "
        "local .gb for — the user's hint is usually phrases like 'human "
        "CGAS gene (look up .gb)', 'fetch from GenBank', or just naming a "
        "gene with no upload. Prefers RefSeqGene (NG_*); falls back to "
        "RefSeq mRNA (NM_*) when no RefSeqGene exists — the response's "
        "`db_source` field tells you which. mRNA fallback has no introns, "
        "so intronic or splice-junction edits are out of scope on that "
        "kind of attachment."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "gene_symbol": {"type": "string", "description": "Gene symbol, e.g. 'CGAS'."},
            "organism": {"type": "string", "description": "Default 'Homo sapiens'. Free-text Linnaean name, e.g. 'Mus musculus'."},
        },
        "required": ["gene_symbol"],
    },
}


# CRISPR + primer design wrappers (delegate to v1 designers, registry-aware).
CRISPR_TOOLS: list[dict[str, Any]] = [
    DESIGN_GUIDES_TOOL,
    DESIGN_PEGRNAS_TOOL,
    DESIGN_PRIMERS_TOOL,
    FIND_GENOMIC_RECORD_TOOL,
]



# Handler table. Append new emitters here as they ship.



# ─────────────────────────────────────────────────────────────────────
# Builder-driven simulate_assembly: parses instruction → IntentSpec,
# annotates inventory attachments to extract role-tagged Parts, runs
# the agent_v2.builder revise loop (max 3 iterations), and returns
# either a verified product attachment or the verifier's diagnostics.
# ─────────────────────────────────────────────────────────────────────
SIMULATE_ASSEMBLY_MAX_ITERS = 3


# Map annotation feature_class / name keywords → builder role taxonomy
# (see agent_v2.builder.intent_spec.KNOWN_ROLES).
_FEATURE_CLASS_TO_ROLE = {
    "promoter": "promoter",
    "cds": "cds",
    "protein_coding": "cds",
    "polya": "polya",
    "polya_signal": "polya",
    "polyadenylation_signal": "polya",
    "terminator": "polya",
    "selection_marker": "selection_marker",
    "antibiotic_resistance": "selection_marker",
    "resistance_gene": "selection_marker",
    "origin": "origin",
    "rep_origin": "origin",
    "ori": "origin",
    "ltr": "ltr",
    "lentiviral_cis": "lentiviral_cis",
    "wpre": "wpre",
    "scaffold": "scaffold",
    "grna_scaffold": "scaffold",
    "spacer": "spacer",
    "stuffer": "stuffer",
    "tag": "tag",
    "nls": "tag",
    "kozak": "kozak",
    "polylinker": "polylinker",
    "mcs": "polylinker",
    "att_site": "att_site",
    "enhancer": "enhancer",
}


def _annotation_to_role(ann: dict) -> Optional[str]:
    """Best-effort mapping from an annotated feature to a builder role.
    Returns None when the feature has no recognised role (e.g. random
    misc_feature); those features are skipped."""
    fc = (ann.get("feature_class") or "").lower()
    if fc in _FEATURE_CLASS_TO_ROLE:
        return _FEATURE_CLASS_TO_ROLE[fc]
    kb = ann.get("kb_data") or {}
    fc2 = (kb.get("feature_class") or "").lower()
    if fc2 in _FEATURE_CLASS_TO_ROLE:
        return _FEATURE_CLASS_TO_ROLE[fc2]
    name = (ann.get("name") or "").lower()
    if any(k in name for k in ("promoter", " cmv", "ef1", "u6 ", "h1 ", "cag ",
                                "cmv ie", "pgk", "ubc", "sv40 prom")):
        return "promoter"
    if any(k in name for k in ("polya", "poly a", "bgh", "rbglob")):
        return "polya"
    if any(k in name for k in ("origin", "pmb1", "cole1", "f1 ori", " ori",
                                "puc ori")):
        return "origin"
    if any(k in name for k in ("ampr", "kanr", "puror", "neor", "blastr",
                                "hygror", "bla(", " bla ", " cat ", "tetr")):
        return "selection_marker"
    if " ltr" in name or name.endswith("ltr"):
        return "ltr"
    if "scaffold" in name:
        return "scaffold"
    if any(k in name for k in (" nls", "flag-", "ha tag", "myc tag", "his tag")):
        return "tag"
    if "kozak" in name:
        return "kozak"
    if any(k in name for k in ("mcs", "polylinker")):
        return "polylinker"
    if any(k in name for k in ("cas9", "cas12", "mcherry", "egfp", "gfp", "rfp",
                                "bfp", "luciferase", "cre recombinase", "abe8",
                                "be4", "p2a", "t2a", "f2a", "e2a")):
        return "cds"
    return None


def _parse_intent(instruction: str):
    """Parse the user's free-text assembly instruction into an
    IntentSpec the builder + verifier can check against."""
    from agent_v2.builder.intent_spec import IntentSpec
    il = (instruction or "").lower()
    host = "mammalian"
    if any(k in il for k in ("e. coli", "escherichia", "bacterial expression")):
        host = "bacterial"
    elif any(k in il for k in ("yeast", "saccharomyces", "s. cerevisiae")):
        host = "yeast"
    if any(k in il for k in ("sgrna cloning", "guide cloning",
                              "pol3-gg", "stuffer-scaffold")):
        return IntentSpec.for_sgrna_cloning_backbone()
    preferred: list[str] = []
    for kw in ("CMV", "EF1a", "EF1α", "U6", "H1", "CAG", "Kozak",
                "Cas9", "Cas12", "mCherry", "GFP", "EGFP", "BGH polyA",
                "SV40 polyA", "WPRE", "P2A", "T2A", "F2A", "E2A", "NLS",
                "AmpR", "PuroR", "NeoR", "BlastR", "HygroR", "pMB1",
                "ColE1", "f1 ori", "pUC57", "pcDNA3.1", "lentiCRISPR"):
        if kw.lower() in il:
            preferred.append(kw)
    return IntentSpec.for_expression_cassette(host_scope=host, preferred=preferred)


async def simulate_assembly_fast(
    args,
    registry,
    *,
    output_dir=None,
):
    """Build a verified VirtualConstruct from inventory attachments.

    Pipeline:
      1. Parse `instruction` → IntentSpec (function + required roles).
      2. For each inventory / target attachment, annotate the source
         plasmid, then slice each role-tagged annotation into a Part
         (with 50 bp upstream + downstream junctions).
      3. Hand the resulting PartSet to agent_v2.builder.build() with
         max_iters=SIMULATE_ASSEMBLY_MAX_ITERS. The builder lays out
         slots in DEFAULT_ORDER, then iteratively reorders / reorients
         / swaps parts in response to verifier diagnostics.
      4. On success: register the materialized sequence as a product
         attachment and return method + slot summary.
      5. On failure: return ok=False with the unresolved diagnostics so
         the caller can see exactly which roles / modules / interactions
         are missing or wrong.

    No empty-args fallback. If no source attachments are supplied, the
    builder still runs against an empty PartSet so the diagnostics
    surface the missing-role errors explicitly.
    """
    from splicify_api.annotation_cache import annotate_llm_cached
    from agent_v2.builder.builder import build
    from agent_v2.builder.part_set import Part, PartSet, extract_part

    target_id = args.get("target_attachment_id")
    inv_ids = args.get("inventory_attachment_ids") or []
    instruction = (args.get("instruction") or "")

    src_ids: list[str] = []
    if target_id:
        src_ids.append(target_id)
    for i in inv_ids:
        if i and i not in src_ids:
            src_ids.append(i)

    intent = _parse_intent(instruction)

    # Annotate each source plasmid and extract role-tagged Parts.
    extracted_parts: list[Part] = []
    annotation_summaries: list[dict] = []
    for aid in src_ids:
        att = registry.get(aid)
        if not att or not att.sequence:
            annotation_summaries.append(
                {"attachment_id": aid, "error": "attachment not found or empty"}
            )
            continue
        try:
            env = await annotate_llm_cached(att.sequence, circular=att.circular)
        except Exception as e:
            annotation_summaries.append(
                {"attachment_id": aid,
                 "error": f"annotation failed: {type(e).__name__}: {e}"}
            )
            continue
        anns = env.get("annotations") or env.get("features") or []
        parts_from_src: list[Part] = []
        for ann in anns:
            role = _annotation_to_role(ann)
            if not role:
                continue
            p = extract_part(aid, att.sequence, ann, role=role,
                              circular=att.circular)
            if p:
                parts_from_src.append(p)
                extracted_parts.append(p)
        annotation_summaries.append({
            "attachment_id": aid,
            "name": att.name,
            "n_annotations": len(anns),
            "n_extracted_parts": len(parts_from_src),
            "roles": sorted({p.role for p in parts_from_src}),
        })

    part_set = PartSet(parts=extracted_parts)

    # Run the builder revise loop. The builder itself enforces
    # max_iters internally; on exhaustion or no-mechanical-fix it
    # returns success=False with the unresolved diagnostics.
    try:
        result = await build(part_set, intent,
                              max_iters=SIMULATE_ASSEMBLY_MAX_ITERS)
    except Exception as e:
        return {
            "ok": False,
            "error": f"builder crashed: {type(e).__name__}: {e}",
            "intent": intent.to_dict(),
            "part_set_summary": part_set.to_dict(),
            "annotation_summaries": annotation_summaries,
        }

    if not result.success:
        diags = [d.to_dict() for d in (result.unresolved_diagnostics or [])]
        missing_roles = sorted({
            (d.get("structured_action") or {}).get("missing_role")
            for d in diags
            if (d.get("structured_action") or {}).get("missing_role")
        })
        missing_modules = sorted({
            (d.get("structured_action") or {}).get("missing_module")
            for d in diags
            if (d.get("structured_action") or {}).get("missing_module")
        })
        return {
            "ok": False,
            "error": (
                f"builder failed after {SIMULATE_ASSEMBLY_MAX_ITERS} iterations: "
                f"{len(diags)} unresolved diagnostic(s)"
            ),
            "intent": intent.to_dict(),
            "part_set_summary": part_set.to_dict(),
            "annotation_summaries": annotation_summaries,
            "diagnostics": diags,
            "missing_roles": missing_roles,
            "missing_modules": missing_modules,
            "method_pick": result.method_pick,
            "method_assessment": (result.method_assessment.to_dict()
                                    if result.method_assessment else None),
            "journal": [{"iteration": e.iteration, "action": e.action,
                          "detail": e.detail} for e in result.journal],
            "hint": (
                "Each diagnostic.structured_action names the role/module "
                "that's missing. Register more source plasmids "
                "(lookup_kb_part / find_external_part / find_genomic_record) "
                "to fill them, then retry simulate_assembly."
            ),
        }

    # Success: materialize and register the product.
    sequence = (result.final_construct.materialize() if result.final_construct else "")
    method_pick = result.method_pick or "gibson"
    label = "_".join(s[:20].replace(" ", "_") for s in (intent.preferred_features[:3]
                                                          or [intent.function]))
    product_name = f"assembled_{label}"
    aid = registry.register_product(name=product_name, sequence=sequence,
                                      circular=intent.topology == "circular")

    return {
        "ok": True,
        "method_hint": method_pick,
        "product_attachment_id": aid,
        "product_name": product_name,
        "length_bp": len(sequence),
        "topology": intent.topology,
        "fragment_count": len(result.final_construct.slots),
        "fragments": [
            {"name": s.part.name, "role": s.part.role,
              "length_bp": s.part.length_bp,
              "orientation": s.orientation,
              "source_plasmid_id": s.part.source_plasmid_id}
            for s in result.final_construct.slots
        ],
        "intent": intent.to_dict(),
        "verifier_summary": (result.final_verification.annotation_summary
                              if result.final_verification else {}),
        "method_assessment": (result.method_assessment.to_dict()
                                if result.method_assessment else None),
        "annotation_summaries": annotation_summaries,
        "notes": (
            "Built via agent_v2.builder revise-loop (verified against "
            "IntentSpec). Pass product_attachment_id to verify_assembly + "
            "annotate_attachment + emit_assembled_gb."
        ),
    }


FIND_EXTERNAL_PART_MAX_PER_CALL = 2
# Lowered from 6 → 3 (2026-05-20): with the coverage-score / stop_searching
# signal returned by find_external_part_cloning, the agent should converge
# on a parent vector in 1-2 calls. The smaller budget forces it to commit
# instead of forever fishing for "a better match".
FIND_EXTERNAL_PART_MAX_PER_SESSION = 3


async def find_external_part_cloning(
    args,
    registry,
    *,
    output_dir=None,
):
    """Cloning-side find_external_part: search Addgene + auto-download
    + annotate + register in the AttachmentRegistry so the Main agent
    can pass attachment_id to annotate_attachment / simulate_assembly."""
    from splicify_api.external_search import (
        search_addgene, fetch_addgene_entry, download_addgene_gb,
    )
    from splicify_api.annotation_cache import annotate_llm_cached
    from splicify_api.agent.agent_tools import extract_seq_from_genbank

    description = (args.get("description") or "").strip()
    if not description:
        return {"ok": False, "error": "description is required"}
    max_results = max(1, min(FIND_EXTERNAL_PART_MAX_PER_CALL,
                              int(args.get("max_results") or 1)))
    existing = int(getattr(registry, "_external_lookups_done", 0) or 0)
    remaining = FIND_EXTERNAL_PART_MAX_PER_SESSION - existing
    if remaining <= 0:
        return {"ok": False, "candidates": [], "n_candidates": 0,
                "error": (f"External-lookup budget exhausted "
                            f"({FIND_EXTERNAL_PART_MAX_PER_SESSION} calls); "
                            "ask the user to upload the .gb directly.")}
    max_results = min(max_results, remaining)
    setattr(registry, "_external_lookups_done", existing + max_results)

    # Optional `required_features`: a list of keyword phrases the user
    # wants the parent to contain (extracted from the design brief).
    # When supplied, we compute per-candidate coverage (case-insensitive
    # substring against module_types + annotation feature names) and
    # include a `coverage_score` ∈ [0,1] in each candidate dict so the
    # LLM can stop searching as soon as one candidate scores ≥ 0.75.
    required_features = args.get("required_features") or []
    required_features = [
        s.strip().lower() for s in required_features
        if isinstance(s, str) and s.strip()
    ]

    def _coverage(envelope: dict) -> tuple[float, list[str], list[str]]:
        if not required_features:
            return (0.0, [], [])
        modules = envelope.get("modules") or []
        feats = envelope.get("annotations") or envelope.get("features") or []
        haystack = " | ".join(
            (m.get("module_type") or "").lower() for m in modules
        ) + " | " + " | ".join(
            (f.get("name") or "").lower() for f in feats
        )
        matched = [kw for kw in required_features if kw in haystack]
        missing = [kw for kw in required_features if kw not in haystack]
        score = len(matched) / len(required_features)
        return (round(score, 3), matched, missing)

    # Addgene's search punishes long descriptive queries (e.g.
    # "third-generation lentiviral all-in-one CRISPR vector" → 0 hits).
    # First try the caller's query; if it returns 0 candidates,
    # synthesize shorter fallback queries from `required_features` and
    # retry until we find SOMETHING (or run out of fallbacks). The
    # fallback budget is internal — does not consume the per-session
    # external-lookup cap.
    candidates = await search_addgene(description, max_results=max_results)
    queries_tried = [description]
    if not candidates and required_features:
        # Build short fallback queries from required_features. Strip
        # underscored module names (e.g. "lentiviral_payload") and
        # tokenize. Prefer compact 2-3 word combos.
        kw_tokens: list[str] = []
        for kw in required_features:
            for token in kw.replace("_", " ").split():
                t = token.strip().lower()
                if len(t) >= 3 and t not in {"the", "and", "of", "for"} \
                        and t not in kw_tokens:
                    kw_tokens.append(t)
        # Generate a small ladder of progressively shorter queries.
        candidate_queries = []
        if len(kw_tokens) >= 2:
            candidate_queries.append(" ".join(kw_tokens[:3]))
        if len(kw_tokens) >= 1:
            candidate_queries.append(kw_tokens[0])
        # Also try the LAST descriptive token in description (often the
        # plasmid name itself).
        descr_tokens = [t for t in description.split() if len(t) >= 4]
        if descr_tokens:
            candidate_queries.append(descr_tokens[-1])
            if len(descr_tokens) >= 2:
                candidate_queries.append(" ".join(descr_tokens[-2:]))
        # Dedup + run
        seen_q = {description.lower()}
        for q in candidate_queries:
            if q.lower() in seen_q:
                continue
            seen_q.add(q.lower())
            queries_tried.append(q)
            candidates = await search_addgene(q, max_results=max_results)
            if candidates:
                break
    out = []
    best_score = 0.0
    for cand in candidates:
        entry = await fetch_addgene_entry(cand.addgene_id)
        d = cand.to_dict()
        if entry:
            d.update({
                "description": entry.description, "depositor": entry.depositor,
                "pmid": entry.pmid, "doi": entry.doi,
                "paper_title": entry.paper_title,
                "sequences_page": f"https://www.addgene.org/{cand.addgene_id}/sequences/",
            })
        try:
            gb_text = await download_addgene_gb(cand.addgene_id)
            if gb_text:
                seq = extract_seq_from_genbank(gb_text)
                if seq and len(seq) >= 50:
                    envelope = await annotate_llm_cached(seq, circular=True)
                    envelope.setdefault("sequence", seq)
                    aid = registry.register_product(
                        name=(entry.name if entry else cand.name), sequence=seq, circular=True,
                    )
                    _ANNOTATION_CACHE[aid] = envelope
                    d["attachment_id"] = aid
                    d["length_bp"] = len(seq)
                    d["module_types"] = sorted({
                        m.get("module_type") for m in (envelope.get("modules") or [])
                        if m.get("module_type")
                    })[:10]
                    score, matched, missing = _coverage(envelope)
                    d["coverage_score"] = score
                    d["coverage_matched"] = matched
                    d["coverage_missing"] = missing
                    best_score = max(best_score, score)
        except Exception as e:
            d["download_error"] = f"{type(e).__name__}: {e}"
        out.append(d)

    result = {
        "ok": True, "source": "addgene", "query": description,
        "queries_tried": queries_tried,
        "n_candidates": len(out), "candidates": out,
    }
    if required_features:
        result["best_coverage_score"] = best_score
        result["coverage_threshold"] = 0.75
        if best_score >= 0.75:
            result["stop_searching"] = True
            # Identify the winning candidate (highest coverage with an
            # attachment_id) and stash it on the registry so cassette_swap
            # can pick it as the default target — even if the agent forgets
            # to thread the id through.
            winners = [c for c in out
                       if c.get("attachment_id")
                       and (c.get("coverage_score") or 0) >= 0.75]
            if winners:
                top = max(winners, key=lambda c: c.get("coverage_score") or 0)
                setattr(registry, "_locked_parent_attachment_id", top["attachment_id"])
                result["locked_parent_attachment_id"] = top["attachment_id"]
            result["stop_reason"] = (
                f"A candidate covered {int(best_score * 100)}% of the "
                f"requested features (≥75% threshold). Stop calling "
                f"find_external_part for this parent and proceed to "
                f"cassette_swap (or replace_region for explicit-bounds "
                f"edits)."
            )
        else:
            result["stop_searching"] = False
    return result


FIND_EXTERNAL_PART_CLONING_TOOL = {
    "name": "find_external_part",
    "description": (
        "Search Addgene + auto-download + annotate + register a candidate "
        "plasmid in the local AttachmentRegistry. Returns the new "
        "attachment_id you can pass to find_features / graft_parts. Use "
        "when the user described a specific external plasmid (by name, "
        "Addgene ID, or distinguishing-feature description) AND "
        "lookup_kb_part returned no high-confidence local match.\n\n"
        "CONVERGENCE RULE — call this AT MOST ONCE per parent vector. "
        "Pass `required_features` (a list of keyword phrases pulled "
        "from the user's design brief: e.g. ['lentiviral_payload', "
        "'guide_expression_cassette', 'cas9', 'puror']). The tool "
        "computes per-candidate coverage against the downloaded "
        "annotation and returns `best_coverage_score` + `stop_searching: "
        "true` once any candidate covers ≥75% of the keywords. When "
        "`stop_searching` is true, IMMEDIATELY stop calling "
        "find_external_part and proceed to find_features + graft_parts "
        "with the returned attachment_id. Do not call again to look "
        "for 'a better match' — the budget is small and the cassette "
        "swap doesn't need a perfect parent."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "description": {"type": "string",
                              "description": ("Free-text query for Addgene. "
                                               "Short, distinguishing phrases "
                                               "work best (e.g. 'lentiCRISPR "
                                               "v2 all-in-one'). Avoid pasting "
                                               "the entire user prompt.")},
            "max_results": {"type": "integer", "default": 1},
            "required_features": {
                "type": "array",
                "items": {"type": "string"},
                "description": ("Keyword phrases the parent vector must "
                                 "contain — used to compute coverage_score "
                                 "and emit the stop_searching signal. "
                                 "Examples: ['cas9', 'puror', 'u6 promoter', "
                                 "'lentiviral_payload', 'wpre']. "
                                 "Case-insensitive substring match against "
                                 "module_types + feature names."),
            },
        },
        "required": ["description"],
    },
}




async def design_primers_batch_tool(
    args,
    registry,
    *,
    output_dir=None,
):
    """Run N design_primers requests in a single LLM turn. Each
    request reuses the underlying primer3 call from the per-region
    tool; batching saves the LLM round-trip per request (~10-15s each)."""
    from agent_v2.crispr_tools import design_primers_tool as _single

    requests = args.get("requests") or []
    if not isinstance(requests, list) or not requests:
        return {"ok": False, "error": "requests must be a non-empty list"}

    results = []
    for i, req in enumerate(requests):
        if not isinstance(req, dict):
            results.append({"ok": False, "error": f"request {i} not a dict"})
            continue
        # Forward each request through the existing single-region handler.
        out = await _single(req, registry, output_dir=output_dir)
        if isinstance(out, dict):
            out["batch_index"] = i
            label = req.get("pair_label") or req.get("application") or f"req_{i}"
            out.setdefault("pair_label", label)
        results.append(out)

    return {
        "ok": True,
        "n_requests": len(requests),
        "n_successful": sum(1 for r in results if isinstance(r, dict) and r.get("ok", True) is not False),
        "results": results,
    }


DESIGN_PRIMERS_BATCH_TOOL = {
    "name": "design_primers_batch",
    "description": (
        "Batch wrapper around design_primers — accepts a list of "
        "per-region requests and returns a list of results. Use this "
        "instead of N separate design_primers calls when you need "
        "primers for multiple regions on the same construct (typical: "
        "NGS amplicon + Sanger amplicon for the same sgRNA cut site, "
        "or one primer pair per part in a Gibson assembly). Saves one "
        "LLM turn per request (~10-15s each).\n\n"
        "sgRNA cut-site defaults: for cut at position P, both NGS "
        "(application='illumina') and Sanger (application='sanger') "
        "use region_start=P-250, region_end=P+250, excluded_start=P-75, "
        "excluded_end=P+75. Product sizes are auto-selected per "
        "application (150-290 for illumina, 250-500 for sanger)."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "requests": {
                "type": "array",
                "description": "List of design_primers requests. Each item is the same shape as design_primers's input: {attachment_id, region_start, region_end, excluded_start?, excluded_end?, application?, product_size_min?, product_size_max?, primer_opt_tm?, num_return?, pair_label?}.",
                "items": {
                    "type": "object",
                    "properties": {
                        "attachment_id":   {"type": "string"},
                        "region_start":    {"type": "integer"},
                        "region_end":      {"type": "integer"},
                        "excluded_start":  {"type": "integer"},
                        "excluded_end":    {"type": "integer"},
                        "application":     {"type": "string", "description": "fragment | sanger (default) | illumina"},
                        "product_size_min": {"type": "integer"},
                        "product_size_max": {"type": "integer"},
                        "primer_opt_tm":   {"type": "number"},
                        "num_return":      {"type": "integer"},
                        "pair_label":      {"type": "string", "description": "Human-readable label for this request (e.g. 'EMX1 cut-site NGS')."},
                    },
                    "required": ["attachment_id", "region_start", "region_end"],
                },
            },
        },
        "required": ["requests"],
    },
}


# ─────────────────────────────────────────────────────────────────────
# graft_parts: deterministic ordered concatenation with dual-orientation
# verification. Smaller hammer than simulate_assembly_fast — caller
# specifies the EXACT part order; this tool just glues them together
# and checks the annotated result. Designed for cassette swap / graft
# workflows where the user has already decided the layout.
# ─────────────────────────────────────────────────────────────────────
import re


def _rc_dna(seq: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
            "a": "t", "t": "a", "g": "c", "c": "g", "n": "n"}
    return "".join(comp.get(b, "N") for b in reversed(seq or ""))


async def _resolve_graft_part(
    part: dict, registry, *,
    kb_organism: str = "h_sapiens",
) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Return (sequence, name, source_label) for one graft part spec.

    Accepted spec shapes:
      {"attachment_id": "att_X"}                     — use registered sequence
      {"kb_part_name": "mCherry"}                    — look up + back-translate
      {"sequence": "ACGT..."}                        — caller-supplied DNA
    """
    aid = part.get("attachment_id")
    if aid:
        att = registry.get(aid)
        if not att:
            return None, None, f"unknown attachment_id={aid!r}"
        seq = att.sequence
        # Sub-region extraction: when caller provides start/end (1-based,
        # inclusive — matches GenBank convention), slice that region.
        # Supports circular wrap (start > end) by concatenating the end
        # of the plasmid with the beginning.
        start = part.get("start")
        end = part.get("end")
        if start is not None and end is not None:
            L = len(seq)
            try:
                s = max(1, int(start)) - 1
                e = int(end)
            except (TypeError, ValueError):
                return None, None, f"bad start/end on attachment {aid}"
            if s < 0 or e <= 0 or e > L:
                if att.circular and 1 <= int(start) <= L and 1 <= int(end) <= L \
                        and int(start) > int(end):
                    # Circular wrap: take [start..L] + [1..end]
                    seq = seq[s:L] + seq[:e]
                else:
                    return None, None, (f"start/end out of bounds for "
                                          f"attachment {aid} (len={L})")
            elif s >= e:
                if att.circular:
                    seq = seq[s:L] + seq[:e]
                else:
                    return None, None, (f"start ({start}) must be < end ({end}) "
                                          f"for linear attachment {aid}")
            else:
                seq = seq[s:e]
        return seq, (part.get("name") or att.name), \
               f"attachment:{aid}" + (f"[{start}..{end}]" if start and end else "")

    raw = part.get("sequence")
    if raw:
        clean = "".join(c for c in raw.upper() if c in "ACGTN")
        if not clean:
            return None, None, "sequence had no DNA characters"
        return clean, (part.get("name") or "raw_sequence"), "literal_sequence"

    kb_name = part.get("kb_part_name") or part.get("kb_name")
    if kb_name:
        try:
            from splicify_api.agent.agent_tools import tool_lookup_kb_part
        except Exception as e:
            return None, None, f"kb tool unavailable: {type(e).__name__}: {e}"
        res = await tool_lookup_kb_part(
            {"name": kb_name, "register_attachment": True,
              "organism": part.get("organism") or kb_organism},
            registry,
        )
        matches = res.get("matches") or []
        top = matches[0] if matches else None
        if not top or not top.get("attachment_id"):
            return None, None, f"kb lookup found no DNA for {kb_name!r}"
        att = registry.get(top["attachment_id"])
        if not att:
            return None, None, f"kb attachment {top['attachment_id']} missing"
        return att.sequence, (part.get("name") or top.get("name") or kb_name), \
               f"kb:{top.get('name') or kb_name} ({top.get('provenance')})"

    return None, None, "part spec needs attachment_id, sequence, or kb_part_name"


def _detect_cassette_swap_signature(parts: list[dict]) -> Optional[dict]:
    """Inspect a graft_parts `parts` list to see if it matches an
    in-place edit pattern (parent_5' + new_content + parent_3' with
    the same parent attachment_id on both flanks). When detected,
    returns a dict with the inferred replace_region args. Otherwise
    returns None.

    Signature requirements:
      - len(parts) >= 3
      - parts[0]  has `attachment_id` + `start` + `end`
      - parts[-1] has `attachment_id` + `start` + `end`
      - parts[0].attachment_id == parts[-1].attachment_id
      - parts[0].end + 1 <= parts[-1].start (slices don't overlap or invert)
      - At least one middle part is NOT a slice of the same parent
        (otherwise this is just a 3-piece parent reassembly, no edit)
    """
    if not isinstance(parts, list) or len(parts) < 3:
        return None
    head, tail = parts[0], parts[-1]
    middle = parts[1:-1]
    if not (isinstance(head, dict) and isinstance(tail, dict)):
        return None
    head_aid = head.get("attachment_id")
    tail_aid = tail.get("attachment_id")
    if not head_aid or head_aid != tail_aid:
        return None
    head_end = head.get("end")
    tail_start = tail.get("start")
    if not isinstance(head_end, int) or not isinstance(tail_start, int):
        return None
    if head_end + 1 > tail_start:
        return None
    # At least one middle part must NOT be a same-parent slice.
    middle_is_all_same_parent = all(
        isinstance(m, dict)
        and m.get("attachment_id") == head_aid
        and isinstance(m.get("start"), int)
        and isinstance(m.get("end"), int)
        for m in middle
    )
    if middle_is_all_same_parent:
        return None
    return {
        "parent_attachment_id": head_aid,
        "region_start": head_end + 1,
        "region_end": tail_start - 1,
        "replacement_parts": middle,
    }


async def graft_parts(
    args,
    registry,
    *,
    output_dir=None,
):
    """Concatenate an ordered list of parts (in the order given), run
    the annotator on the result, and — when the assembled construct
    fails annotation expectations — retry with every part marked
    `is_insert: true` flipped to its reverse complement.

    AUTO-ROUTE: when the parts list matches a cassette-swap signature
    (parent_5' slice + new content + parent_3' slice of the SAME
    attachment), this function delegates to `replace_region` so the
    parent stays intact and only the conflicting region is swapped.
    The redirection is logged in the returned `auto_routed_to_replace_region`
    field so callers can see it.

    Args:
      parts: ordered list of part specs (see _resolve_graft_part).
        Each may also carry `orientation: 1|-1` (default 1) and
        `is_insert: bool` (default false). When `try_both_orientations`
        is true (default), insert-flagged parts are flipped if the
        first attempt fails verification.
      topology: "circular" | "linear" (default "circular").
      try_both_orientations: bool (default true).
      required_roles: list[str] — roles the annotated product must
        contain (e.g. ["promoter", "cds", "polya"]). Defaults to
        ["promoter", "cds"] if any insert part name includes "cds"
        keywords, else no role check.
      product_name: optional label for the registered attachment.
    """
    # Cassette-swap auto-detect: if the parts list is parent_5' + new
    # content + parent_3' of the same attachment, delegate to
    # replace_region. This implements rule C from the replacement
    # rule set: in-place edits stay in-place; graft_parts is reserved
    # for true de-novo construction. The `_skip_auto_route` flag
    # prevents infinite recursion when replace_region itself calls
    # graft_parts internally.
    parts_in = args.get("parts") or []
    swap_sig = (None if args.get("_skip_auto_route")
                else _detect_cassette_swap_signature(parts_in))
    if swap_sig:
        rr_args = {
            "target_attachment_id": swap_sig["parent_attachment_id"],
            "region_start": swap_sig["region_start"],
            "region_end": swap_sig["region_end"],
            "replacement_parts": swap_sig["replacement_parts"],
            "topology": args.get("topology"),
            "try_both_orientations": args.get("try_both_orientations", True),
            "required_roles": args.get("required_roles") or [],
            "product_name": args.get("product_name"),
            "organism": args.get("organism"),
        }
        # No preserve_module — the caller didn't tell us which module
        # to protect, so just run the edit and let required_roles catch
        # anything important. (When the agent calls replace_region
        # directly, it should pass preserve_module.)
        result = await replace_region(rr_args, registry, output_dir=output_dir)
        result["auto_routed_to_replace_region"] = True
        result["auto_route_reason"] = (
            f"graft_parts received a cassette-swap signature "
            f"(parent {swap_sig['parent_attachment_id']} on both flanks, "
            f"region {swap_sig['region_start']}..{swap_sig['region_end']}). "
            f"Delegated to replace_region. Pass preserve_module on the "
            f"next call to also enforce a module-survival check."
        )
        return result
    from splicify_api.annotation_cache import annotate_llm_cached
    from splicify_api.agent.agent_tools import tool_lookup_kb_part  # noqa: F401

    parts_in = args.get("parts") or []
    if not isinstance(parts_in, list) or not parts_in:
        return {"ok": False, "error": "parts must be a non-empty ordered list"}
    topology = args.get("topology", "circular")
    try_both = bool(args.get("try_both_orientations", True))
    required_roles = args.get("required_roles") or []
    kb_organism = args.get("organism") or "h_sapiens"

    # Resolve each part to a concrete sequence.
    resolved: list[dict] = []
    for i, spec in enumerate(parts_in):
        if not isinstance(spec, dict):
            return {"ok": False, "error": f"part #{i} is not a dict"}
        seq, name, source = await _resolve_graft_part(
            spec, registry, kb_organism=kb_organism,
        )
        if not seq:
            return {
                "ok": False,
                "error": f"part #{i} could not be resolved: {source}",
                "parts_resolved_so_far": [r["name"] for r in resolved],
            }
        resolved.append({
            "index": i,
            "name": name,
            "source": source,
            "length_bp": len(seq),
            "orientation": int(spec.get("orientation") or 1),
            "is_insert": bool(spec.get("is_insert", False)),
            "_sequence": seq,
        })

    def _materialize(plan: list[dict]) -> str:
        out = []
        for r in plan:
            s = r["_sequence"]
            if r["orientation"] == -1:
                s = _rc_dna(s)
            out.append(s)
        return "".join(out).upper()

    async def _attempt(plan: list[dict], label: str) -> dict:
        """Materialize, annotate, evaluate against required_roles."""
        seq = _materialize(plan)
        try:
            env = await annotate_llm_cached(seq, circular=topology == "circular")
        except Exception as e:
            return {"label": label, "ok": False, "annotation_error":
                    f"{type(e).__name__}: {e}", "length_bp": len(seq),
                    "sequence": seq}
        feats = env.get("annotations") or env.get("features") or []
        modules = env.get("modules") or []
        roles_seen = set()
        # Feature-class aliases — multiple annotators use slightly
        # different vocabularies for the same role.
        _CLASS_ALIASES = {
            "polya_signal": "polya", "polyadenylation_signal": "polya",
            "terminator": "polya",
            "rep_origin": "origin", "ori": "origin",
            "selection_marker": "selection_marker",
            "antibiotic_resistance": "selection_marker",
            "protein_coding": "cds",
        }
        for f in feats:
            fc = (f.get("feature_class") or "").lower()
            if fc:
                roles_seen.add(_CLASS_ALIASES.get(fc, fc))
            name = (f.get("name") or "").lower()
            # Token-based name matching (tolerates "poly(A)", "bGH polyA", etc).
            n_collapsed = re.sub(r"[^a-z0-9]", "", name)
            for needle, role in (
                ("promoter", "promoter"), ("polya", "polya"),
                ("polyasignal", "polya"), ("polyadenylation", "polya"),
                ("terminator", "polya"),
                ("ori", "origin"), ("origin", "origin"),
                ("cas9", "cds"), ("mcherry", "cds"), ("gfp", "cds"),
                ("rfp", "cds"), ("bfp", "cds"), ("yfp", "cds"),
                ("ampr", "selection_marker"), ("kanr", "selection_marker"),
                ("puror", "selection_marker"), ("bleor", "selection_marker"),
                ("neor", "selection_marker"), ("hygror", "selection_marker"),
            ):
                if needle in n_collapsed:
                    roles_seen.add(role)
        missing = [r for r in required_roles if r not in roles_seen]
        return {
            "label": label, "ok": not missing, "length_bp": len(seq),
            "sequence": seq, "missing_roles": missing,
            "roles_seen": sorted(roles_seen),
            "module_types": sorted({m.get("module_type") for m in modules
                                      if m.get("module_type")}),
            "n_features": len(feats),
            "annotation_envelope": env,  # kept so the caller can stash it in _ANNOTATION_CACHE
        }

    attempts: list[dict] = [await _attempt(resolved, "forward")]
    chosen = attempts[0] if attempts[0]["ok"] else None

    if not chosen and try_both and any(r["is_insert"] for r in resolved):
        flipped = [dict(r) for r in resolved]
        for r in flipped:
            if r["is_insert"]:
                r["orientation"] = -r["orientation"]
        attempts.append(await _attempt(flipped, "insert_flipped"))
        if attempts[-1]["ok"]:
            chosen = attempts[-1]
            resolved = flipped

    # If still no pass and required_roles wasn't supplied, accept the
    # forward attempt as long as we got any annotation back. The caller
    # can validate further via verify_assembly.
    if not chosen and not required_roles and attempts[0].get("length_bp"):
        chosen = attempts[0]

    if not chosen:
        return {
            "ok": False,
            "error": "graft did not pass required_roles check in any orientation",
            "attempts": [{k: v for k, v in a.items() if k != "sequence"}
                          for a in attempts],
            "parts_resolved": [{k: v for k, v in r.items() if not k.startswith("_")}
                                 for r in resolved],
        }

    seq = chosen["sequence"]
    label = args.get("product_name") or "graft_" + "_".join(
        r["name"].replace(" ", "_")[:16] for r in resolved[:3]
    )
    aid = registry.register_product(
        name=label, sequence=seq, circular=topology == "circular",
    )

    # Stash the verification-time annotation under the new attachment_id
    # so find_features / replace_region / annotate_attachment can read
    # from cache instead of paying the annotation cost again.
    env = chosen.get("annotation_envelope")
    if isinstance(env, dict):
        env["sequence"] = seq
        _ANNOTATION_CACHE[aid] = env

    return {
        "ok": True,
        "product_attachment_id": aid,
        "product_name": label,
        "length_bp": len(seq),
        "topology": topology,
        "orientation_used": chosen["label"],
        "fragments": [{"name": r["name"], "length_bp": r["length_bp"],
                         "orientation": r["orientation"],
                         "is_insert": r["is_insert"],
                         "source": r["source"]} for r in resolved],
        "roles_seen": chosen.get("roles_seen", []),
        "module_types": chosen.get("module_types", []),
        "n_features_annotated": chosen.get("n_features", 0),
        "attempts": [{k: v for k, v in a.items() if k not in ("sequence", "annotation_envelope")}
                      for a in attempts],
        "notes": (
            "Sequential ordered concatenation; product registered as a new "
            "attachment. Annotation cached so downstream tools (find_features, "
            "emit_assembled_gb, replace_region's module check) hit O(1)."
        ),
    }


GRAFT_PARTS_TOOL = {
    "name": "graft_parts",
    "description": (
        "Deterministic ordered-concatenation assembler. Use this when "
        "you already know the exact part order for a cassette swap, "
        "graft, or de-novo build, and want the tool to glue them "
        "together + annotate + return a product attachment_id. "
        "Each part can be {attachment_id} (registered plasmid or "
        "fragment), {kb_part_name} (KB lookup + automatic protein "
        "back-translation if needed), or {sequence} (caller-supplied "
        "DNA). Mark the cassette to insert with `is_insert: true`; "
        "if the first orientation fails the required-roles check, the "
        "tool retries with all insert parts reverse-complemented. "
        "Prefer this over simulate_assembly for cassette swaps where "
        "you've already decided the layout — it's faster and "
        "deterministic."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "parts": {
                "type": "array",
                "description": (
                    "Ordered list of parts. Each item: at least one of "
                    "{attachment_id, kb_part_name, sequence}, plus "
                    "optional {name, orientation: 1|-1, is_insert: bool, "
                    "organism (for back-translation when kb_part_name "
                    "is protein-only; default h_sapiens), start, end "
                    "(1-based inclusive sub-region of an attachment — "
                    "use to splice in a cassette region; circular wrap "
                    "is supported when start > end)}."
                ),
                "items": {"type": "object"},
                "minItems": 1,
            },
            "topology": {"type": "string", "enum": ["circular", "linear"],
                          "default": "circular"},
            "try_both_orientations": {"type": "boolean", "default": True},
            "required_roles": {
                "type": "array", "items": {"type": "string"},
                "description": (
                    "Roles the annotated product must contain "
                    "(e.g. ['promoter','cds','polya']). When empty, "
                    "no role check — only annotation success is needed."
                ),
            },
            "product_name": {"type": "string"},
            "organism": {"type": "string",
                          "description": "Default organism for KB back-translation."},
        },
        "required": ["parts"],
    },
}


# ─────────────────────────────────────────────────────────────────────
# replace_region: in-place edit of an existing plasmid. Treats the
# parent as a single unit that stays intact except for ONE region,
# which gets excised and replaced with new content. Wraps graft_parts
# internally so orientation-flip + role verification are reused, then
# checks that the module containing the original region survives.
# ─────────────────────────────────────────────────────────────────────

# Module-type → expected role checklist. Used when the caller asks to
# preserve a module after the edit but doesn't enumerate the roles.
_MODULE_REQUIRED_ROLES = {
    "mammalian_pol2_expression_cassette": ["promoter", "cds", "polya"],
    "bacterial_expression_cassette":      ["promoter", "cds"],
    "guide_expression_cassette":          ["promoter", "scaffold"],
    "bacterial_selection_cassette":       ["selection_marker", "origin"],
    "mammalian_selection_cassette":       ["promoter", "selection_marker", "polya"],
    "lentiviral_payload":                 ["promoter", "cds", "polya"],
}


async def replace_region(args, registry, *, output_dir=None):
    """In-place edit: keep the parent plasmid intact and swap ONE region
    for new content.

    Use this for cassette swaps where the user wants to modify a part of
    an existing plasmid (e.g. replace PuroR with mCherry inside the
    EF-1α-driven Cas9 cassette of lentiCRISPR v2). The tool builds the
    final construct as:

        parent[1 .. region_start-1]  +  <replacement_parts>  +  parent[region_end+1 .. end]

    and delegates to graft_parts internally, so the existing dual-
    orientation retry + annotation-based role check still apply.

    After the graft succeeds, the tool re-checks the product's annotation
    to confirm the module that originally contained the excised region
    still appears in the product (i.e. the replacement integrated
    without breaking the surrounding cassette). On failure it returns
    ok=False with module_preserved=false so the caller can adjust.

    Args:
        target_attachment_id (str): plasmid to edit (must be registered).
        region_start (int, 1-based inclusive): start of the region to excise.
        region_end (int, 1-based inclusive): end of the region to excise.
        replacement_parts (list): same shape as graft_parts.parts. Each
            entry is auto-marked is_insert=true if not already set.
        preserve_module (str, optional): module_type the product MUST
            still contain after the edit (e.g.
            "mammalian_pol2_expression_cassette"). If not given,
            required_roles is empty and only annotation success is
            checked.
        required_roles (list, optional): explicit role list, overrides
            the preserve_module → roles mapping.
        try_both_orientations (bool, default True): forwarded to graft_parts.
        product_name (str, optional): label for the registered product
            attachment.
    """
    aid = args.get("target_attachment_id")
    rs = args.get("region_start")
    re_end = args.get("region_end")
    replacement = args.get("replacement_parts") or []

    if not aid:
        return {"ok": False, "error": "target_attachment_id is required"}
    if not isinstance(rs, int) or not isinstance(re_end, int):
        return {"ok": False,
                "error": "region_start and region_end must be 1-based integers"}
    if not replacement or not isinstance(replacement, list):
        return {"ok": False, "error": "replacement_parts must be a non-empty list"}

    att = registry.get(aid)
    if not att:
        return {"ok": False, "error": f"unknown attachment_id={aid!r}"}
    L = len(att.sequence)

    # Guard #4 — target_attachment_id must be a parent backbone, not
    # a KB-registered fragment. KB hits (named "kb_…") and small
    # linear sequences are inputs to splice IN, not targets to edit.
    is_kb_stub = (att.name or "").startswith("kb_")
    is_small_linear = (not att.circular) and L < 3000
    if is_kb_stub or is_small_linear:
        # Suggest the largest circular attachment as the right target.
        candidates = []
        for other_aid, other_att in (
            (registry.items or {}).items()
            if hasattr(registry, "items") else []
        ):
            if other_aid == aid:
                continue
            if (other_att.name or "").startswith("kb_"):
                continue
            if not other_att.circular:
                continue
            candidates.append({
                "attachment_id": other_aid,
                "name": other_att.name,
                "length_bp": len(other_att.sequence),
            })
        candidates.sort(key=lambda c: -c["length_bp"])
        reason = ("looks like a KB-registered fragment (name starts with 'kb_')"
                  if is_kb_stub else
                  f"is linear and only {L} bp")
        return {
            "ok": False,
            "error": (
                f"target_attachment_id={aid!r} {reason}. Use the parent "
                f"plasmid backbone as target instead — typically the "
                f"largest circular attachment registered (the find_external_part "
                f"download, or the user's upload). replacement_parts is where "
                f"KB hits / fragments go; target is the plasmid being edited."
            ),
            "suggested_targets": candidates[:5],
            "hint": (
                "Re-call replace_region with target_attachment_id = the "
                "largest circular attachment from suggested_targets."
            ),
        }

    if rs < 1 or re_end > L or rs > re_end:
        return {"ok": False,
                "error": (f"invalid region: start={rs}, end={re_end}, "
                           f"parent length={L}")}

    # Guard #1 — rebuild anti-pattern. When the excision is zero-bp
    # (region_end == region_start - 1 ≡ insert at point) AND
    # replacement_parts contain slices of the SAME parent, the agent is
    # almost certainly trying to rebuild the parent inside the
    # replacement instead of editing it. Refuse with diagnostics.
    excised_bp = re_end - rs + 1
    parent_slices_in_replacement = [
        r for r in replacement
        if isinstance(r, dict)
        and r.get("attachment_id") == aid
        and isinstance(r.get("start"), int)
        and isinstance(r.get("end"), int)
    ]
    if excised_bp <= 1 and parent_slices_in_replacement:
        return {
            "ok": False,
            "error": (
                f"rebuild anti-pattern: excision is only {excised_bp} bp but "
                f"replacement_parts contains {len(parent_slices_in_replacement)} "
                f"slice(s) of the same parent ({aid}). This means you're trying "
                f"to reconstruct the parent vector inside the replacement region "
                f"instead of editing it in place. Use find_cassette_for + "
                f"find_features to get real cassette boundaries, then pass them "
                f"as region_start..region_end (the actual span to excise) and "
                f"keep only NEW content (literal sequences or KB hits without "
                f"parent_match) in replacement_parts. The parent flanks are "
                f"supplied automatically by replace_region."
            ),
            "replaced_region": {
                "parent_attachment_id": aid,
                "region_start": rs, "region_end": re_end,
                "excised_bp": excised_bp,
                "parent_slices_in_replacement": [
                    {"start": s.get("start"), "end": s.get("end")}
                    for s in parent_slices_in_replacement
                ],
            },
            "hint": "Call find_cassette_for(target_attachment_id, query=<cassette anchor>) before retrying.",
        }

    # Guard #2 — annotation-anchored bounds. Cassette swaps should
    # cut AT feature / submodule boundaries (within ±25 bp). If both
    # region_start AND region_end are deep interior coordinates (far
    # from any annotated boundary) AND the excision is non-trivial
    # (>50 bp), the agent likely guessed the bounds. Refuse with the
    # nearest valid boundary candidates surfaced so it can retry.
    if excised_bp > 50:
        env = await _get_or_annotate(att)
        if env:
            boundary_positions: list[tuple[int, str, str]] = []
            for f in (env.get("annotations") or env.get("features") or []):
                fs, fe = f.get("start"), f.get("end")
                if isinstance(fs, int):
                    boundary_positions.append((fs, "start", f.get("name") or "?"))
                if isinstance(fe, int):
                    boundary_positions.append((fe, "end", f.get("name") or "?"))
            for m in env.get("modules") or []:
                ms, me = m.get("start"), m.get("end")
                mn = m.get("module_type") or m.get("name") or "?"
                if isinstance(ms, int):
                    boundary_positions.append((ms, "start", f"[module] {mn}"))
                if isinstance(me, int):
                    boundary_positions.append((me, "end", f"[module] {mn}"))
            for sm in (m.get("submodules") or []) if env.get("modules") else []:
                pass  # submodule boundaries are covered by find_features

            if boundary_positions:
                def _nearest(pos: int) -> tuple[int, str, str]:
                    return min(boundary_positions, key=lambda x: abs(x[0] - pos))
                near_start = _nearest(rs)
                near_end = _nearest(re_end)
                # Allow ±25 bp tolerance to either boundary
                start_ok = abs(near_start[0] - rs) <= 25
                end_ok = abs(near_end[0] - re_end) <= 25
                if not start_ok and not end_ok:
                    # Suggest a few nearby valid boundaries
                    nearby = sorted(boundary_positions,
                                     key=lambda x: min(abs(x[0] - rs),
                                                         abs(x[0] - re_end)))[:8]
                    return {
                        "ok": False,
                        "error": (
                            f"region_start ({rs}) and region_end ({re_end}) "
                            f"are not within 25 bp of any annotated feature or "
                            f"module boundary in parent {aid}. The nearest "
                            f"boundary to region_start is "
                            f"{near_start[1]} of {near_start[2]!r} at "
                            f"position {near_start[0]} (Δ={abs(near_start[0]-rs)} bp); "
                            f"to region_end it's {near_end[1]} of "
                            f"{near_end[2]!r} at position {near_end[0]} "
                            f"(Δ={abs(near_end[0]-re_end)} bp). "
                            f"Cassette swaps should cut at real feature "
                            f"boundaries. Call find_cassette_for then use "
                            f"submodule start/end values directly."
                        ),
                        "replaced_region": {
                            "parent_attachment_id": aid,
                            "region_start": rs, "region_end": re_end,
                            "excised_bp": excised_bp,
                        },
                        "nearby_boundaries": [
                            {"position": p, "side": side, "name": name}
                            for (p, side, name) in nearby
                        ],
                        "hint": ("find_cassette_for(target_attachment_id, "
                                   "query=<feature inside the cassette>) returns "
                                   "submodule start/end you can use directly."),
                    }

    preserve = args.get("preserve_module")
    required_roles = list(args.get("required_roles") or [])
    if preserve and not required_roles:
        required_roles = list(_MODULE_REQUIRED_ROLES.get(preserve, []))

    # Build the 3-region graft: parent_5' + replacement... + parent_3'
    parts: list[dict] = []
    if rs > 1:
        parts.append({
            "attachment_id": aid, "start": 1, "end": rs - 1,
            "name": f"parent_5prime[1..{rs - 1}]",
        })
    # Auto-mark replacement parts as is_insert (so the orientation-flip
    # retry only flips the insert, leaving the parent slices fixed).
    for r in replacement:
        if not isinstance(r, dict):
            return {"ok": False,
                    "error": "every replacement_parts entry must be a dict"}
        if "is_insert" not in r:
            r = {**r, "is_insert": True}
        parts.append(r)
    if re_end < L:
        parts.append({
            "attachment_id": aid, "start": re_end + 1, "end": L,
            "name": f"parent_3prime[{re_end + 1}..{L}]",
        })

    topology = args.get("topology", "circular" if att.circular else "linear")
    product_name = (args.get("product_name")
                    or f"edited_{att.name}_at_{rs}-{re_end}")

    g = await graft_parts({
        "parts": parts,
        "topology": topology,
        "try_both_orientations": args.get("try_both_orientations", True),
        "required_roles": required_roles,
        "product_name": product_name,
        "organism": args.get("organism"),
        # Skip the cassette-swap auto-detect — replace_region IS the
        # cassette-swap path, the parts list we just built would
        # re-trigger detection and recurse forever.
        "_skip_auto_route": True,
    }, registry, output_dir=output_dir)

    edit_info = {
        "parent_attachment_id": aid,
        "parent_length_bp": L,
        "region_start": rs,
        "region_end": re_end,
        "excised_bp": re_end - rs + 1,
        "parent_5prime_bp": rs - 1,
        "parent_3prime_bp": L - re_end,
        "n_replacement_parts": len(replacement),
        "preserve_module": preserve,
        "required_roles_checked": required_roles,
    }

    if not g.get("ok"):
        return {**g, "replaced_region": edit_info}

    # Module-preservation check: read the product annotation (cached by
    # graft_parts) and confirm preserve_module still appears.
    module_preserved = None
    parent_modules_at_region: list[str] = []
    if preserve:
        product_env = _ANNOTATION_CACHE.get(g.get("product_attachment_id")) or {}
        product_mods = product_env.get("modules") or []
        module_types = {m.get("module_type") for m in product_mods}
        module_preserved = preserve in module_types

        # For context, list which modules in the parent overlapped the
        # excised region — helps the caller see what was supposed to be
        # preserved vs what was inadvertently disrupted.
        parent_env = _ANNOTATION_CACHE.get(aid) or {}
        for m in parent_env.get("modules") or []:
            ms = m.get("start") or 0
            me = m.get("end") or 0
            if max(ms, rs) <= min(me, re_end):
                mt = m.get("module_type")
                if mt and mt not in parent_modules_at_region:
                    parent_modules_at_region.append(mt)

    return {
        **g,
        "replaced_region": edit_info,
        "module_preserved": module_preserved,
        "parent_modules_overlapping_region": parent_modules_at_region,
    }


# ─────────────────────────────────────────────────────────────────────
# cassette_swap: high-level deterministic workflow. The LLM extracts
# intent ("keep Cas9, replace what comes after it with P2A-mCherry");
# this tool does all the mechanical work — pick the parent, run
# find_cassette_for, compute boundaries from submodule coords, call
# replace_region. One LLM call replaces a 15-step ReAct chain that
# kept picking the wrong attachment_id or the wrong cut points.
# ─────────────────────────────────────────────────────────────────────
def _auto_pick_parent(registry) -> Optional[str]:
    """Pick the largest circular non-KB attachment as the parent. Used
    when the caller doesn't specify target_attachment_id explicitly."""
    items = (registry.items if hasattr(registry, "items") else {}) or {}
    best_aid = None
    best_len = 0
    for aid, att in items.items():
        if not getattr(att, "circular", False):
            continue
        if (att.name or "").startswith("kb_"):
            continue
        n = len(att.sequence or "")
        if n > best_len:
            best_len = n
            best_aid = aid
    return best_aid


async def cassette_swap(args, registry, *, output_dir=None):
    """End-to-end deterministic cassette swap.

    Hardcodes the mechanical decisions the LLM kept getting wrong:
      1. Pick the parent (largest circular plasmid if not specified).
      2. Run find_cassette_for(parent, cassette_anchor) → cassette + submodules.
      3. Resolve keep_through_feature.end → region_start.
         Resolve excise_through_feature.end → region_end.
      4. Call replace_region with computed bounds + the caller's
         replacement_parts and preserve_module derived from the cassette.

    Returns the replace_region envelope plus a `decisions` field
    documenting every coordinate / boundary the tool resolved, so the
    caller can see exactly what was computed.

    Args:
      target_attachment_id (str, optional): explicit parent. If omitted,
        the largest circular non-KB attachment is auto-picked.
      cassette_anchor (str, required): feature name that identifies which
        cassette to edit (e.g. 'Cas9', 'PuroR', 'EGFP').
      keep_through (str, required): feature inside the cassette to
        KEEP up to and including. Excision starts at this feature's
        end + 1. Pass the cassette_anchor itself when the swap is
        immediately after the anchor (e.g. 'Cas9' for a Cas9-X swap).
      excise_through (str, required): feature inside the cassette to
        DELETE up to and including. Excision ends at this feature's
        end. Pass the LAST feature you want gone (e.g. 'PuroR').
      replacement_parts (list, required): same shape as replace_region —
        ordered list of new content to splice in.
      product_name (str, optional): label for the registered product.
    """
    cassette_anchor = (args.get("cassette_anchor") or "").strip()
    keep_through = (args.get("keep_through") or "").strip()
    excise_through = (args.get("excise_through") or "").strip()
    replacement_parts = args.get("replacement_parts") or []
    if not (cassette_anchor and keep_through and excise_through):
        return {"ok": False,
                "error": "cassette_anchor, keep_through, and excise_through are required"}
    if not replacement_parts:
        return {"ok": False, "error": "replacement_parts must be non-empty"}

    # Step 1 — pick parent. Priority order:
    #   (a) explicit target_attachment_id from the caller
    #   (b) the parent locked by find_external_part's coverage check
    #       (stashed on registry._locked_parent_attachment_id)
    #   (c) auto-pick the largest circular non-KB attachment
    aid = args.get("target_attachment_id")
    auto_source = "explicit"
    if not aid:
        locked = getattr(registry, "_locked_parent_attachment_id", None)
        if locked and registry.get(locked):
            aid = locked
            auto_source = "locked_by_find_external_part"
    if not aid:
        aid = _auto_pick_parent(registry)
        auto_source = "largest_circular"
    if not aid:
        return {"ok": False,
                "error": "no parent attachment found; upload a backbone or call find_external_part first"}
    att = registry.get(aid)
    if not att:
        return {"ok": False, "error": f"unknown attachment_id={aid!r}"}

    # Step 2 — cassette + submodules via find_cassette_for
    cassette_info = await find_cassette_for_tool({
        "target_attachment_id": aid, "query": cassette_anchor,
    }, registry)
    if not cassette_info.get("ok"):
        return {"ok": False,
                "error": (f"find_cassette_for({cassette_anchor!r}) failed: "
                           f"{cassette_info.get('error')}"),
                "parent_attachment_id": aid}
    cassettes = cassette_info.get("cassettes") or []
    if not cassettes:
        return {"ok": False,
                "error": (f"no cassette containing {cassette_anchor!r} found in "
                           f"parent {aid}"),
                "find_cassette_result": cassette_info}
    # Prefer mammalian_pol2_expression_cassette over broader lentiviral_payload etc.
    cassettes.sort(key=lambda c: (
        0 if "expression_cassette" in (c.get("module_type") or "") else 1,
        (c.get("cassette_end") or 0) - (c.get("cassette_start") or 0),
    ))
    cassette = cassettes[0]
    preserve_module = cassette.get("module_type")

    # Step 3 — resolve boundaries from the cassette's submodules + the
    # broader feature annotation envelope (find_features-style lookup).
    env = await _get_or_annotate(att)
    all_feats = (env.get("annotations") or env.get("features") or []) if env else []
    sub_feats = cassette.get("submodules") or []

    def _find_by_name(needle: str) -> Optional[dict]:
        q = needle.lower()
        # Try cassette submodules first (more specific), then full annotation list.
        for src in (sub_feats, all_feats):
            for f in src:
                name = (f.get("name") or "").lower()
                if q in name:
                    return f
        return None

    keep_f = _find_by_name(keep_through)
    excise_f = _find_by_name(excise_through)
    decisions = {
        "parent_attachment_id": aid,
        "parent_pick_source": auto_source,
        "parent_name": att.name,
        "parent_length_bp": len(att.sequence),
        "cassette_anchor": cassette_anchor,
        "cassette_module_type": preserve_module,
        "cassette_bounds": (cassette.get("cassette_start"), cassette.get("cassette_end")),
        "keep_through_feature": keep_f,
        "excise_through_feature": excise_f,
    }
    if not keep_f or not excise_f:
        missing = []
        if not keep_f:
            missing.append(f"keep_through={keep_through!r}")
        if not excise_f:
            missing.append(f"excise_through={excise_through!r}")
        return {"ok": False,
                "error": (f"could not resolve feature(s): {', '.join(missing)} "
                           f"in parent {aid}. Submodules in cassette: "
                           f"{[s.get('name') for s in sub_feats]}"),
                "decisions": decisions,
                "available_features": [f.get("name") for f in all_feats][:30]}

    region_start = (keep_f.get("end") or 0) + 1
    region_end = excise_f.get("end") or 0
    decisions["computed_region_start"] = region_start
    decisions["computed_region_end"] = region_end
    decisions["excised_bp"] = region_end - region_start + 1

    if region_end < region_start:
        return {"ok": False,
                "error": (f"computed region_end ({region_end}) < region_start "
                           f"({region_start}); keep_through is downstream of "
                           f"excise_through — check your feature ordering."),
                "decisions": decisions}

    # Step 4 — call replace_region with computed bounds
    result = await replace_region({
        "target_attachment_id": aid,
        "region_start": region_start,
        "region_end": region_end,
        "replacement_parts": replacement_parts,
        "preserve_module": preserve_module,
        "product_name": (args.get("product_name")
                           or f"swapped_{att.name}_{cassette_anchor}-to-{excise_through}"),
        "try_both_orientations": args.get("try_both_orientations", True),
    }, registry, output_dir=output_dir)

    result["decisions"] = decisions
    return result


CASSETTE_SWAP_TOOL = {
    "name": "cassette_swap",
    "description": (
        "Deterministic end-to-end cassette swap. Use this AS THE "
        "PRIMARY assembly tool for any 'replace X with Y inside an "
        "existing plasmid' workflow — it hardcodes the mechanical "
        "decisions (parent selection, find_cassette_for, boundary "
        "computation, replace_region call) so you only have to extract "
        "the design intent.\n\n"
        "Workflow:\n"
        "  1. The tool auto-picks the parent (largest circular non-KB "
        "attachment) unless you pass target_attachment_id explicitly.\n"
        "  2. It calls find_cassette_for(parent, cassette_anchor) to "
        "locate the cassette containing the anchor feature.\n"
        "  3. It computes region_start = (keep_through feature).end + 1 "
        "and region_end = (excise_through feature).end — both resolved "
        "against the cassette's submodule list (then the full feature "
        "annotation as fallback).\n"
        "  4. It calls replace_region with those bounds + your "
        "replacement_parts + the cassette's module_type as "
        "preserve_module.\n\n"
        "Pass cassette_anchor = the feature that identifies the "
        "cassette (typically the CDS you want to KEEP intact, e.g. "
        "'Cas9'). Pass keep_through = the last feature you want "
        "preserved (often the same as cassette_anchor when the swap "
        "is immediately after it). Pass excise_through = the last "
        "feature to DELETE (e.g. 'PuroR' for a PuroR → mCherry "
        "swap)."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "target_attachment_id": {"type": "string",
                                       "description": ("Optional. Largest circular "
                                                         "non-KB attachment is auto-picked "
                                                         "when omitted.")},
            "cassette_anchor": {"type": "string",
                                  "description": ("Feature name that identifies "
                                                    "which cassette to edit "
                                                    "(e.g. 'Cas9').")},
            "keep_through": {"type": "string",
                              "description": ("Feature inside the cassette to "
                                                "KEEP up to and including. "
                                                "Excision starts at this "
                                                "feature's end + 1.")},
            "excise_through": {"type": "string",
                                "description": ("Feature inside the cassette to "
                                                  "DELETE up to and including. "
                                                  "Excision ends at this feature's "
                                                  "end.")},
            "replacement_parts": {
                "type": "array",
                "description": ("Ordered list of NEW content to splice in. "
                                  "Same shape as replace_region's replacement_parts."),
                "items": {"type": "object"},
                "minItems": 1,
            },
            "product_name": {"type": "string"},
            "try_both_orientations": {"type": "boolean", "default": True},
        },
        "required": ["cassette_anchor", "keep_through",
                       "excise_through", "replacement_parts"],
    },
}


REPLACE_REGION_TOOL = {
    "name": "replace_region",
    "description": (
        "In-place edit of an existing plasmid. PREFER this over "
        "graft_parts for cassette swaps and any 'replace X with Y inside "
        "an existing vector' workflow. The parent plasmid is treated as "
        "ONE unit that stays intact except for the supplied region — "
        "this avoids accidentally dropping intergenic / regulatory "
        "sequences that aren't named as features.\n\n"
        "Workflow:\n"
        "  1. find_features on the parent to locate the boundaries of "
        "the conflicting region (e.g. the start of NLS-FLAG-PuroR "
        "and the end of PuroR / start of polyA).\n"
        "  2. Build replacement_parts as the ordered list of NEW "
        "content (literal sequences, KB-looked-up CDSes, etc.). "
        "Existing parent content is supplied AUTOMATICALLY by this "
        "tool as parent_5prime / parent_3prime flanks — do NOT add "
        "the parent backbone to replacement_parts.\n"
        "  3. Pass preserve_module so the tool verifies the "
        "surrounding cassette survives the edit (e.g. "
        "'mammalian_pol2_expression_cassette' for a Cas9-cassette swap).\n"
        "  4. The tool runs graft_parts internally, retries with "
        "flipped insert(s) if the first orientation fails, and "
        "returns module_preserved=true/false."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "target_attachment_id": {"type": "string"},
            "region_start": {"type": "integer",
                              "description": "1-based inclusive start of the region to excise."},
            "region_end": {"type": "integer",
                            "description": "1-based inclusive end of the region to excise."},
            "replacement_parts": {
                "type": "array",
                "description": (
                    "Ordered list of NEW content to insert in place of "
                    "the excised region. Each item is one of "
                    "{attachment_id}, {kb_part_name, organism?}, "
                    "{sequence}; optional {name, orientation: 1|-1}. "
                    "Parts are auto-marked is_insert=true so the "
                    "dual-orientation retry only flips them, leaving "
                    "the parent flanks fixed."
                ),
                "items": {"type": "object"},
                "minItems": 1,
            },
            "preserve_module": {
                "type": "string",
                "description": (
                    "Module type that should still be detected in the "
                    "product (e.g. 'mammalian_pol2_expression_cassette'). "
                    "Drives the required_roles check via a built-in "
                    "module → roles mapping."
                ),
            },
            "required_roles": {
                "type": "array", "items": {"type": "string"},
                "description": ("Explicit role list; overrides the "
                                  "preserve_module → roles mapping."),
            },
            "try_both_orientations": {"type": "boolean", "default": True},
            "topology": {"type": "string", "enum": ["circular", "linear"]},
            "product_name": {"type": "string"},
            "organism": {"type": "string",
                          "description": "Default organism for KB back-translation."},
        },
        "required": ["target_attachment_id", "region_start",
                      "region_end", "replacement_parts"],
    },
}


# ─────────────────────────────────────────────────────────────────────
# Per-attachment annotation cache. Used by find_features / graft_parts
# / find_external_part / upload paths so a plasmid is annotated AT
# MOST ONCE per process lifetime; subsequent boundary lookups are
# instant. Populated by:
#   - router._add_genbank_to_registry (upload path, eager)
#   - find_external_part_cloning (Addgene download, eager)
#   - _get_or_annotate (lazy fallback when a tool needs annotation)
# ─────────────────────────────────────────────────────────────────────
_ANNOTATION_CACHE: dict[str, dict] = {}


def cache_annotation(attachment_id: str, envelope: dict) -> None:
    """Public helper for non-agent_v2 code (router, orchestrator) to
    seed the annotation cache when they already have an envelope.
    Idempotent — calling twice with the same id is safe."""
    if attachment_id and isinstance(envelope, dict):
        envelope.setdefault("sequence", "")  # downstream lookups expect the key
        _ANNOTATION_CACHE[attachment_id] = envelope


async def _get_or_annotate(att) -> dict:
    """Return the cached annotation envelope for an Attachment, running
    annotate_llm_cached on first request and stashing the result. Keyed
    on attachment_id; safe to call from any tool that needs feature
    coordinates."""
    if att is None or not att.sequence:
        return {}
    cached = _ANNOTATION_CACHE.get(att.attachment_id)
    if cached is not None:
        return cached
    from splicify_api.annotation_cache import annotate_llm_cached
    env = await annotate_llm_cached(att.sequence, circular=att.circular)
    if isinstance(env, dict):
        # Make sequence retrievable for the interpreter's PlasmidIndex.
        env.setdefault("sequence", att.sequence)
        _ANNOTATION_CACHE[att.attachment_id] = env
    return env or {}


async def find_features_tool(args, registry, *, output_dir=None):
    """Locate a feature by name within a registered attachment and
    return its 1-based start/end coordinates (plus direction and KB
    metadata). Wraps the interpreter pipeline's PlasmidIndex.find_features.

    Use the returned coordinates as `start`/`end` in a follow-up
    `graft_parts` call when you need to splice a region in or out of an
    existing plasmid (cassette swap, deletion, insertion).
    """
    aid = args.get("attachment_id")
    query = (args.get("query") or "").strip()
    if not aid or not query:
        return {"ok": False, "error": "attachment_id and query are required"}
    att = registry.get(aid)
    if att is None:
        return {"ok": False, "error": f"unknown attachment_id={aid!r}"}

    env = await _get_or_annotate(att)
    if not env:
        return {"ok": False, "error": "annotation failed or returned empty"}

    from agent_v2.interpreter.plasmid_index import PlasmidIndex
    idx = PlasmidIndex.from_envelope(aid, env, name=att.name)

    # Try features first (CDS / promoter / etc.), then modules
    # (cassettes), then cloning features (restriction sites).
    feats = idx.find_features(query)
    mods = idx.find_modules(query)
    cf = idx.find_cloning_features(query)

    def _slim(d: dict) -> dict:
        return {k: d.get(k) for k in
                ("name", "start", "end", "direction", "strand",
                  "module_type", "rule_id", "type",
                  "feature_family", "subtype", "description")
                if d.get(k) is not None}

    return {
        "ok": True,
        "attachment_id": aid,
        "attachment_name": att.name,
        "attachment_length_bp": len(att.sequence),
        "query": query,
        "features": [_slim(f) for f in feats[:10]],
        "modules":  [_slim(m) for m in mods[:10]],
        "cloning_features": [_slim(c) for c in cf[:10]],
        "n_features": len(feats),
        "n_modules": len(mods),
        "n_cloning_features": len(cf),
    }


async def find_cassette_for_tool(args, registry, *, output_dir=None):
    """Return the functional cassette (module) that contains a named
    feature, with its boundaries + the ordered list of submodules
    inside it. Use this BEFORE replace_region to pick the right
    excision region: replace the WHOLE cassette, replace only one
    submodule, or replace a region between two submodules.

    Wraps the interpreter pipeline's PlasmidIndex (find_features +
    find_modules + expression_cassette_for) so a single call answers
    "what's the cassette I need to edit?"
    """
    aid = args.get("target_attachment_id") or args.get("attachment_id")
    query = (args.get("query") or args.get("feature_name") or "").strip()
    if not aid or not query:
        return {"ok": False,
                "error": "target_attachment_id and query are required"}
    att = registry.get(aid)
    if att is None:
        return {"ok": False, "error": f"unknown attachment_id={aid!r}"}

    env = await _get_or_annotate(att)
    if not env:
        return {"ok": False, "error": "annotation failed or returned empty"}

    from agent_v2.interpreter.plasmid_index import PlasmidIndex
    idx = PlasmidIndex.from_envelope(aid, env, name=att.name)

    feats = idx.find_features(query)
    if not feats:
        # Fall back to module-name search
        mods_direct = idx.find_modules(query)
        if not mods_direct:
            return {"ok": False, "attachment_id": aid,
                    "query": query,
                    "error": f"no feature or module matching {query!r}"}
        return {
            "ok": True, "attachment_id": aid, "attachment_name": att.name,
            "query": query,
            "matched_via": "module_name",
            "cassettes": [{
                "module_type": m.get("module_type"),
                "rule_id": m.get("rule_id"),
                "module_name": m.get("name"),
                "cassette_start": m.get("start"),
                "cassette_end": m.get("end"),
                "strand": m.get("strand"),
                "submodules": [
                    {"name": s.get("name"), "start": s.get("start"),
                      "end": s.get("end"),
                      "submodule_type": s.get("submodule_type")}
                    for s in (m.get("submodules") or [])
                ],
            } for m in mods_direct[:5]],
        }

    # Build cassette entries: for each matching feature, find the
    # module(s) that fully contain it.
    cassettes: list[dict] = []
    seen: set[tuple] = set()
    for feat in feats[:5]:
        fs = feat.get("start")
        fe = feat.get("end")
        if not isinstance(fs, int) or not isinstance(fe, int):
            continue
        containing = []
        for m in idx.modules():
            ms = m.get("start")
            me = m.get("end")
            if not isinstance(ms, int) or not isinstance(me, int):
                continue
            if ms <= fs and me >= fe:
                containing.append(m)
        if not containing:
            continue
        # Pick the SMALLEST (most-specific) module first
        containing.sort(key=lambda m: (m.get("end") or 0) - (m.get("start") or 0))
        for m in containing[:2]:
            key = (m.get("module_type"), m.get("start"), m.get("end"))
            if key in seen:
                continue
            seen.add(key)
            cassettes.append({
                "trigger_feature": feat.get("name"),
                "trigger_feature_start": fs,
                "trigger_feature_end": fe,
                "module_type": m.get("module_type"),
                "rule_id": m.get("rule_id"),
                "module_name": m.get("name"),
                "cassette_start": m.get("start"),
                "cassette_end": m.get("end"),
                "strand": m.get("strand"),
                "submodules": [
                    {"name": s.get("name"), "start": s.get("start"),
                      "end": s.get("end"),
                      "submodule_type": s.get("submodule_type")}
                    for s in (m.get("submodules") or [])
                ],
            })

    # Also surface expression_cassette_for() when applicable — it
    # explicitly identifies promoter + polyA boundaries, which is the
    # most common cassette-swap shape.
    cassette_pair = idx.expression_cassette_for(query)

    return {
        "ok": True,
        "attachment_id": aid,
        "attachment_name": att.name,
        "query": query,
        "matched_via": "feature",
        "n_matching_features": len(feats),
        "matching_features": [
            {"name": f.get("name"), "start": f.get("start"),
              "end": f.get("end"),
              "direction": f.get("direction") or f.get("strand")}
            for f in feats[:5]
        ],
        "cassettes": cassettes,
        "expression_cassette": cassette_pair,
        "guidance": (
            "For a WHOLE-cassette swap, pass cassette_start..cassette_end "
            "to replace_region. For a partial swap (e.g. keep Cas9 but "
            "replace the post-Cas9 cassette tail), pick the trigger_feature.end + 1 "
            "as region_start and the start of the next downstream submodule "
            "(or the cassette_end) as region_end. expression_cassette.promoter.end and "
            "expression_cassette.polyA.start give clean upstream and downstream "
            "anchors when the cassette is a Pol II expression unit."
        ),
    }


FIND_CASSETTE_FOR_TOOL = {
    "name": "find_cassette_for",
    "description": (
        "Return the functional cassette (module) containing a named "
        "feature in a registered plasmid, plus its sub-feature list "
        "and the explicit promoter / polyA boundaries when the "
        "cassette is a Pol II expression unit. PREFER this over raw "
        "find_features when you need to decide WHERE to cut for a "
        "replace_region call — find_features gives you a single "
        "feature's coords, but cassette swaps usually need the "
        "boundaries of the surrounding functional unit. The response "
        "includes a `guidance` field explaining how to map cassette / "
        "submodule coordinates onto replace_region's region_start / "
        "region_end inputs."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "target_attachment_id": {"type": "string"},
            "query": {"type": "string",
                       "description": ("Case-insensitive substring against "
                                        "feature names (e.g. 'Cas9') or "
                                        "module names (e.g. 'mammalian_pol2_expression_cassette').")},
        },
        "required": ["target_attachment_id", "query"],
    },
}


async def ask_plasmid_tool(args, registry, *, output_dir=None):
    """Free-form natural-language question over one or more registered
    plasmids. Delegates to the interpreter ReAct loop (run_interpreter)
    which has 8 deterministic plasmid-analysis tools. Returns the
    natural-language answer plus the trace so the caller can extract
    coordinates from the tool calls.

    Use this as a fallback when find_cassette_for / find_features
    can't resolve a design question (e.g. "where should I cut to add a
    second sgRNA cassette without disrupting the existing Cas9
    expression?"). Slower than the deterministic tools because it
    runs its own LLM loop.
    """
    from agent_v2.interpreter.plasmid_registry import PlasmidRegistry
    from agent_v2.interpreter.agent import run_interpreter

    question = (args.get("question") or "").strip()
    if not question:
        return {"ok": False, "error": "question is required"}
    aids = args.get("attachment_ids")
    if aids is None and args.get("target_attachment_id"):
        aids = [args["target_attachment_id"]]
    if not aids:
        # Default: every registered plasmid that looks like a real
        # backbone (>= 1 kb). Skips KB stubs.
        aids = [a.attachment_id for a in (registry.items.values()
                                            if hasattr(registry, "items") else [])
                if len(a.sequence) >= 1000]
    if not aids:
        return {"ok": False, "error": "no registered plasmids to query"}

    pr = PlasmidRegistry()
    for aid in aids:
        att = registry.get(aid)
        if att is None:
            continue
        env = await _get_or_annotate(att)
        if not env:
            continue
        env["sequence"] = att.sequence
        pr.register(aid, env, name=att.name)

    if pr.n() == 0:
        return {"ok": False, "error": "no annotations available for the requested attachments"}

    try:
        result = await run_interpreter(question, pr,
                                         max_iters=int(args.get("max_iters") or 6))
    except Exception as e:
        return {"ok": False,
                "error": f"interpreter crashed: {type(e).__name__}: {e}"}

    # Extract any coordinates surfaced by the interpreter's tool calls.
    extracted_coords: list[dict] = []
    for entry in (result.trace or []):
        res = entry.get("result") or {}
        for k in ("matches", "features", "modules", "cloning_features", "cassettes"):
            for item in (res.get(k) or []):
                if isinstance(item, dict) and "start" in item and "end" in item:
                    extracted_coords.append({
                        "from_tool": entry.get("tool"),
                        "name": item.get("name"),
                        "start": item.get("start"),
                        "end": item.get("end"),
                        "module_type": item.get("module_type"),
                    })

    return {
        "ok": True,
        "answer": result.answer,
        "citations": result.citations,
        "n_tool_calls": result.n_tool_calls,
        "extracted_coords": extracted_coords[:25],
        "plasmid_ids_queried": list(pr.items.keys()) if hasattr(pr, "items") else aids,
    }


ASK_PLASMID_TOOL = {
    "name": "ask_plasmid",
    "description": (
        "Free-form natural-language question over registered plasmids. "
        "Delegates to the interpreter ReAct loop with 8 plasmid-analysis "
        "sub-tools. Use this as a fallback when find_cassette_for / "
        "find_features don't resolve a design question — e.g. "
        "'which submodule of the Cas9 cassette should I excise to "
        "replace PuroR with mCherry without disrupting the bGH polyA?' "
        "Returns the natural-language answer plus extracted_coords for "
        "any start/end pairs surfaced during the loop. SLOWER than the "
        "deterministic tools (runs its own LLM iterations), so reach for "
        "find_cassette_for first."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "question": {"type": "string"},
            "target_attachment_id": {"type": "string",
                                       "description": "Scope the question to one plasmid."},
            "attachment_ids": {
                "type": "array", "items": {"type": "string"},
                "description": ("Scope the question to a specific list. "
                                  "When omitted and target_attachment_id is "
                                  "also omitted, every registered plasmid "
                                  ">= 1 kb is included."),
            },
            "max_iters": {"type": "integer", "default": 6,
                            "description": "Interpreter loop iteration cap."},
        },
        "required": ["question"],
    },
}


FIND_FEATURES_TOOL = {
    "name": "find_features",
    "description": (
        "Look up the 1-based start/end coordinates of a named feature, "
        "module, or restriction site inside a registered attachment. "
        "Use this to discover the boundaries of a cassette you want to "
        "swap out (e.g. find Cas9, PuroR, EF-1-alpha) before passing "
        "the coords as start/end inputs to graft_parts. The first call "
        "for an attachment lazily annotates it via annotate_llm_cached; "
        "subsequent calls hit a process-local cache."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id": {"type": "string"},
            "query": {"type": "string",
                       "description": ("Case-insensitive substring against "
                                        "feature/module/cloning-feature names "
                                        "(e.g. 'Cas9', 'PuroR', 'WPRE', "
                                        "'EF-1-alpha', 'BsmBI').")},
        },
        "required": ["attachment_id", "query"],
    },
}


# simulate_assembly_fast remains importable for direct programmatic use
# (and so old runs / tests don't break), but it is NOT exposed to the
# main agent anymore. graft_parts is the only assembly tool the LLM
# sees. To re-enable simulate_assembly later, re-add it to
# EMITTER_HANDLERS and stop filtering it out of make_full_tool_roster.
_ARCHIVED_TOOL_NAMES = {"simulate_assembly", "graft_parts"}


# ─────────────────────────────────────────────────────────────────────
# Parent-aware KB lookup: wraps v1's tool_lookup_kb_part to surface
# when the queried feature ALREADY appears in a registered parent
# plasmid's annotation. KB hits are generic / codon-optimized
# sequences that won't match the parent's exact bases — for cassette
# swaps the agent should slice the parent via find_features +
# replace_region rather than registering a parallel KB copy.
# ─────────────────────────────────────────────────────────────────────
async def lookup_kb_part_with_parent_check(args, registry, *, output_dir=None):
    """v2 wrapper around v1's tool_lookup_kb_part. After the v1 call,
    cross-references every registered attachment's cached annotation
    (via `_ANNOTATION_CACHE`) for features whose name contains the
    query token. When a match is found, the response gets a
    `parent_match` list naming the attachment + coordinates plus a
    `parent_match_note` instructing the LLM to slice the parent
    instead of using the KB hit.

    The underlying KB call is unchanged — back-translation, attachment
    auto-registration, and the ranked match list still happen. The
    parent-match info is purely additive context for the agent.
    """
    from splicify_api.agent.agent_tools import tool_lookup_kb_part
    result = await tool_lookup_kb_part(args, registry)
    if not isinstance(result, dict):
        return result
    query_raw = (args.get("name") or "").strip()
    query = query_raw.lower()
    if not query:
        return result

    # Tokenize the query so multi-word features like "BGH polyA" or
    # "EF-1-alpha promoter" still match feature names containing any
    # significant token.
    _STOPWORDS = {"the", "a", "of", "and", "or", "to", "from"}
    tokens = [t for t in re.split(r"[^A-Za-z0-9]+", query)
              if len(t) >= 3 and t not in _STOPWORDS]
    if not tokens:
        tokens = [query]

    parent_matches = []
    seen_aids: set[str] = set()
    items = getattr(registry, "items", None) or {}
    for aid, att in list(items.items()):
        # Skip just-registered KB attachments (the v1 call may have
        # registered the top hit as att_product_X; we don't want it to
        # match itself).
        name = getattr(att, "name", "") or ""
        role = getattr(att, "role", "")
        if role == "product" and name.startswith("kb_"):
            continue
        if aid in seen_aids:
            continue
        env = _ANNOTATION_CACHE.get(aid)
        if not env:
            continue
        feats = env.get("annotations") or env.get("features") or []
        modules = env.get("modules") or []

        def _matches(hay: str) -> bool:
            if not hay:
                return False
            hay_low = hay.lower()
            return any(t in hay_low for t in tokens)

        hit = None
        for f in feats:
            haystack = " | ".join(filter(None, [
                f.get("name"),
                f.get("sseqid"),
                (f.get("kb_data") or {}).get("gene_name"),
                (f.get("kb_data") or {}).get("protein_name"),
            ]))
            if _matches(haystack):
                hit = {
                    "feature_name": f.get("name"),
                    "feature_class": f.get("feature_class")
                                     or (f.get("kb_data") or {}).get("feature_class"),
                    "start": f.get("start"),
                    "end": f.get("end"),
                    "direction": f.get("direction") or f.get("strand"),
                    "source": "feature",
                }
                break
        if not hit:
            for m in modules:
                if _matches(m.get("module_type") or "") or _matches(m.get("name") or ""):
                    hit = {
                        "module_type": m.get("module_type"),
                        "module_name": m.get("name"),
                        "start": m.get("start"),
                        "end": m.get("end"),
                        "source": "module",
                    }
                    break
        if hit:
            hit["parent_attachment_id"] = aid
            hit["parent_name"] = name
            parent_matches.append(hit)
            seen_aids.add(aid)

    if parent_matches:
        result["parent_match"] = parent_matches
        result["parent_match_note"] = (
            f"'{query_raw}' already appears in {len(parent_matches)} "
            f"registered parent attachment(s). "
            f"Prefer find_features + replace_region to slice the parent's "
            f"existing copy at the listed coordinates — the KB hit is a "
            f"generic / codon-optimized version that won't match the "
            f"parent's exact bases. Use the KB hit only if the parent's "
            f"version is functionally wrong (e.g. different organism, "
            f"different isoform) and the design requires substitution."
        )
    return result


EMITTER_HANDLERS: dict[str, Any] = {
    "emit_assembled_gb": emit_assembled_gb,
    "emit_parts_order": emit_parts_order,
    "emit_protocol": emit_protocol,
    "emit_workflow_trace": emit_workflow_trace,
    "emit_guides_csv": emit_guides_csv,
    "emit_guides_gb": emit_guides_gb,
    "resolve_feature_position": resolve_feature_position,
    "design_guides": design_guides_tool,
    "design_pegrnas": design_pegrnas_tool,
    "design_primers": design_primers_tool,
    "find_external_part": find_external_part_cloning,
    "design_primers_batch": design_primers_batch_tool,
    "find_genomic_record": find_genomic_record_tool,
    # NOTE: graft_parts is INTERNAL only — used by replace_region under
    # the hood. The LLM no longer sees it in the tool roster; cassette
    # swaps go through replace_region exclusively. Keeping the function
    # exported (`graft_parts` in module scope) so replace_region can
    # call it; just not in EMITTER_HANDLERS.
    "find_features": find_features_tool,
    "find_cassette_for": find_cassette_for_tool,
    "ask_plasmid": ask_plasmid_tool,
    "replace_region": replace_region,
    "cassette_swap": cassette_swap,
    "lookup_kb_part": lookup_kb_part_with_parent_check,
}


def make_full_tool_roster() -> list[dict[str, Any]]:
    """Main agent's full tool roster: v1's tools minus any in
    _ARCHIVED_TOOL_NAMES, plus the agent_v2 emitters, find_features,
    and replace_region. graft_parts is INTERNAL (no schema exposed)
    — all assembly goes through replace_region."""
    from splicify_api.agent.tool_schemas import AIPLASMIDDESIGN_TOOLS
    v1 = [t for t in AIPLASMIDDESIGN_TOOLS
          if t.get("name") not in _ARCHIVED_TOOL_NAMES]
    return v1 + list(EMITTER_TOOLS) + list(RESOLVER_TOOLS) + list(CRISPR_TOOLS) + [
        FIND_EXTERNAL_PART_CLONING_TOOL,
        DESIGN_PRIMERS_BATCH_TOOL,
        FIND_FEATURES_TOOL,
        FIND_CASSETTE_FOR_TOOL,
        ASK_PLASMID_TOOL,
        CASSETTE_SWAP_TOOL,
        REPLACE_REGION_TOOL,
    ]


# Process-local cache of emit_guides_csv args, keyed by either an explicit
# target_attachment_id or, if absent, by the first per-row target_attachment_id
# found inside one of the arrays. Used to back-fill emit_guides_gb /
# emit_parts_order / emit_protocol when the LLM forgets to repeat the data.
_EMITTER_ARG_CACHE: dict[str, dict[str, Any]] = {}
_EMITTER_CACHE_KEYS = ("guides", "pegrnas", "primers", "cloning_oligos", "descriptor")
# emit_parts_order / emit_protocol have no target_attachment_id in their
# args, so the per-target cache lookup misses. Fall back to the LAST
# descriptor cached by any emit_guides_csv call this process has seen,
# so the filename prefix applies to those emitters too.
_LAST_DESCRIPTOR: dict[str, Any] = {"value": None}

def _cache_key_for(args: dict[str, Any]) -> Optional[str]:
    target = args.get("target_attachment_id")
    if isinstance(target, str) and target:
        return target
    for key in ("pegrnas", "guides", "primers"):
        for item in (args.get(key) or []):
            if isinstance(item, dict) and isinstance(item.get("target_attachment_id"), str):
                return item["target_attachment_id"]
    return None

def _cache_emit_args(args: dict[str, Any]) -> None:
    # Stash the descriptor regardless of cache key — emit_parts_order /
    # emit_protocol have no target_attachment_id to key by, so they fall
    # back to this last-seen descriptor for their filename prefix.
    if isinstance(args.get("descriptor"), str) and args["descriptor"].strip():
        _LAST_DESCRIPTOR["value"] = args["descriptor"]
    key = _cache_key_for(args)
    if not key:
        return
    bucket = _EMITTER_ARG_CACHE.setdefault(key, {})
    for k in _EMITTER_CACHE_KEYS:
        v = args.get(k)
        if v:
            bucket[k] = v

def _fill_from_cache(args: dict[str, Any]) -> dict[str, Any]:
    key = _cache_key_for(args)
    merged = dict(args)
    bucket: Optional[dict[str, Any]] = None
    if key and key in _EMITTER_ARG_CACHE:
        bucket = _EMITTER_ARG_CACHE[key]
    elif _EMITTER_ARG_CACHE:
        # No per-target key matched (emit_parts_order / emit_protocol carry
        # no target_attachment_id in their args). Fall back to the MOST
        # RECENT cached bucket so the auto-derive path in those emitters
        # can still see pegrnas / primers / cloning_oligos / descriptor.
        bucket = list(_EMITTER_ARG_CACHE.values())[-1]
    if bucket:
        for k in _EMITTER_CACHE_KEYS:
            if not merged.get(k) and bucket.get(k):
                merged[k] = bucket[k]
    # Descriptor falls back to the last-cached value when nothing in the
    # per-target bucket matched — needed for emit_parts_order / emit_protocol
    # whose args carry no target_attachment_id.
    if not merged.get("descriptor") and _LAST_DESCRIPTOR.get("value"):
        merged["descriptor"] = _LAST_DESCRIPTOR["value"]
    return merged

_FILL_FROM_CACHE_TOOLS = {"emit_guides_gb", "emit_parts_order", "emit_protocol"}

async def dispatch_with_emitters(
    name: str,
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    """Two-tier dispatch: emitters first, then v1's dispatch_tool.
    Caches emit_guides_csv args and back-fills downstream emitters."""
    if name == "emit_guides_csv":
        _cache_emit_args(args)
    elif name in _FILL_FROM_CACHE_TOOLS:
        args = _fill_from_cache(args)
    handler = EMITTER_HANDLERS.get(name)
    if handler is not None:
        return await handler(args, registry, output_dir=output_dir)
    from splicify_api.agent.agent_tools import dispatch_tool
    return await dispatch_tool(name, args, registry)
