"""
Anthropic tool schemas for the AIPlasmidDesign agent (v2).

The agent never sees raw DNA. Plasmids and assembly products are referenced
by `attachment_id`; tools fetch the underlying sequence server-side and
return high-level properties only (features, modules, cut maps, primer
scores, design verdicts, choice-match hashes).
"""
from __future__ import annotations

# --- AIPlasmidDesign tools (preferred) -----------------------------------

ANNOTATE_ATTACHMENT_TOOL = {
    "name": "annotate_attachment",
    "description": (
        "Annotate an attached plasmid by attachment_id using AIPlasmidDesign's "
        "deterministic annotation pipeline. Returns features (name, type, "
        "start, end, strand), high-level modules, and cloning features "
        "(restriction sites, Gateway att sites). Use this FIRST for any "
        "question about what is on a plasmid."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id": {"type": "string"},
        },
        "required": ["attachment_id"],
    },
}

SIMULATE_ASSEMBLY_TOOL = {
    "name": "simulate_assembly",
    "description": (
        "Simulate a cloning reaction (Gibson, Golden Gate, Gateway, "
        "restriction, SDM) using AIPlasmidDesign's per-intent design "
        "handlers. The assembled product is registered as a NEW attachment "
        "so you can subsequently call annotate_attachment, verify_assembly, "
        "or compare_to_choice on its attachment_id. Returns the new "
        "product_attachment_id, intent classification, and predesign verdict."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "instruction":              {"type": "string"},
            "target_attachment_id":     {"type": "string"},
            "inventory_attachment_ids": {"type": "array", "items": {"type": "string"}},
        },
        "required": ["instruction"],
    },
}

DIGEST_PLASMID_TOOL = {
    "name": "digest_plasmid",
    "description": (
        "Compute restriction-digest fragments for an attached plasmid. "
        "Returns cut positions and fragment lengths. Use for restriction "
        "screens."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id": {"type": "string"},
            "enzymes":       {"type": "array", "items": {"type": "string"}},
            "circular":      {"type": "boolean", "default": True},
        },
        "required": ["attachment_id", "enzymes"],
    },
}

FIND_PRIMER_BINDING_TOOL = {
    "name": "find_primer_binding_sites",
    "description": (
        "Locate where primers bind on an attached template. Returns "
        "binding positions and strand. Use this BEFORE score_sanger_primer "
        "if you don't already know primer binding positions."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "template_attachment_id": {"type": "string"},
            "primers":                {"type": "array", "items": {"type": "string"}},
            "circular":               {"type": "boolean", "default": True},
        },
        "required": ["template_attachment_id", "primers"],
    },
}

SCORE_SANGER_PRIMER_TOOL = {
    "name": "score_sanger_primer",
    "description": (
        "Score one or more candidate primers for Sanger-sequencing quality "
        "across a named target feature. Wraps the application-aware "
        "Sanger scorer (Tm, self-dimer, hairpin, length, GC, GC clamp, "
        "homopolymer, primer-to-target distance, mispriming, template "
        "secondary structure). Returns per-primer overall_score [0-100], "
        "rating, warnings, and breakdown. USE THIS for any 'which primer "
        "for Sanger sequencing across X' question rather than reasoning "
        "about distances yourself."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "template_attachment_id": {"type": "string"},
            "primers": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Candidate primer sequences (forward orientation as written).",
            },
            "target_feature_name": {
                "type": "string",
                "description": (
                    "Name of the feature the read should cover (e.g. 'gRNA TU', "
                    "'sgRNA scaffold', 'mScarlet'). The tool resolves it "
                    "against the annotated features. Provide either this OR "
                    "target_position."
                ),
            },
            "target_position": {
                "type": "integer",
                "description": "0-indexed bp on the template that should fall inside the readable window. Use if you know the position.",
            },
        },
        "required": ["template_attachment_id", "primers"],
    },
}

ANALYZE_DESIGN_INTENT_TOOL = {
    "name": "analyze_design_intent",
    "description": (
        "Run AIPlasmidDesign's design-completeness analyzer over a free-text "
        "cloning request. Returns the inferred intent (e.g. gibson, "
        "golden_gate, sdm), the parts the user named, what expression-feature "
        "kinds are present (promoter, CDS, polyA, tag, etc.), what kinds are "
        "missing, and design_warnings (e.g. 'CDS without promoter'). Use "
        "this for 'what is the purpose of this cloning?' questions or to "
        "sanity-check whether a request will work at all."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "user_message": {
                "type": "string",
                "description": "The cloning request text — the same text the user describing their experiment would write.",
            },
        },
        "required": ["user_message"],
    },
}

VERIFY_ASSEMBLY_TOOL = {
    "name": "verify_assembly",
    "description": (
        "Run AIPlasmidDesign's interaction-driven design verifier over an "
        "attached plasmid (or a simulate_assembly product). Annotates the "
        "sequence, builds the SBO interaction graph, and checks every "
        "expression feature against its expected role (promoter has a "
        "downstream CDS, polyA flanks the CDS, tag sits inside a CDS, "
        "att site has a partner, enhancer pairs with the right promoter). "
        "Returns {passed, warnings, suggestions, summary}. Use this to "
        "answer 'is this assembly correct / what is wrong with this design'."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id": {"type": "string"},
        },
        "required": ["attachment_id"],
    },
}

COMPARE_TO_CHOICE_TOOL = {
    "name": "compare_to_choice",
    "description": (
        "For multiple-choice questions whose options are full DNA sequences "
        "(presented to you as '[redacted DNA, N nt]'), check whether an "
        "attachment matches a particular choice letter. The tool hashes the "
        "attachment's sequence (in every rotation, both strands for circular "
        "plasmids) and compares against the unredacted choice text held "
        "server-side. Returns {match: bool, length_matches: bool, "
        "best_rotation_offset, hamming_distance_estimate}. Use this AFTER "
        "simulate_assembly when the question asks for the 'resulting "
        "plasmid sequence' — call it once per choice letter (A, B, C, D)."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_id": {
                "type": "string",
                "description": "Usually the product attachment from simulate_assembly (e.g. 'att_product_1').",
            },
            "choice_letter": {
                "type": "string",
                "description": "A single letter A-H matching the option in the question.",
            },
        },
        "required": ["attachment_id", "choice_letter"],
    },
}

GOLDEN_GATE_ASSEMBLE_TOOL = {
    "name": "golden_gate_assemble",
    "description": (
        "DETERMINISTIC in-silico Golden Gate assembly. Digests each input "
        "plasmid with the named Type-IIS enzyme (Esp3I/BsmBI/BsaI/BbsI/SapI), "
        "finds compatible 4-nt overhangs, walks the unique cyclic ligation "
        "order, and registers the assembled product as a new attachment. "
        "Use this for 'what is the resulting plasmid sequence?' multi-part "
        "Golden Gate questions BEFORE compare_to_choice. Returns "
        "product_attachment_id (or feasible=false + reason)."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "attachment_ids": {
                "type": "array",
                "items": {"type": "string"},
                "description": "All input plasmid attachment_ids (target + inventory).",
            },
            "enzyme": {
                "type": "string",
                "enum": ["Esp3I", "BsmBI", "BsaI", "BbsI", "SapI"],
            },
        },
        "required": ["attachment_ids", "enzyme"],
    },
}

ROUTE_WORKFLOW_TOOL = {
    "name": "route_workflow",
    "description": (
        "Score every cloning method (Gateway, Gibson, Golden Gate, "
        "restriction, SDM, sgRNA-Golden-Gate, PCR-extension Gibson, "
        "synthesis fallback) against an attached target + inventory. "
        "Returns the winner plus full per-method reports (feasible, score, "
        "rationale). USE THIS for any 'which approach' / 'which method' / "
        "'which protocol' question where the answer is a cloning method, "
        "BEFORE answering. Do not reason about method feasibility from "
        "general biology heuristics — this tool runs the actual feasibility "
        "assessors."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "target_attachment_id": {"type": "string"},
            "inventory_attachment_ids": {
                "type": "array",
                "items": {"type": "string"},
            },
        },
        "required": ["target_attachment_id"],
    },
}

LOOKUP_KB_PART_TOOL = {
    "name": "lookup_kb_part",
    "description": (
        "Search AIPlasmidDesign's knowledge base for a feature or part by "
        "name (e.g. 'MCP', 'MS2 coat protein', 'eGFP', 'CMV promoter'). "
        "If `attachment_id` is supplied, also aligns the KB sequence "
        "against the attachment (both strands, circular-aware) and "
        "reports the fraction of the KB part present plus a qualitative "
        "label ('about 1/3', 'about 2/3', 'essentially complete'). "
        "USE THIS when a question references a feature/protein name that "
        "the annotation pipeline does not recognise (e.g. MCP, NLS variants, "
        "less-common tags), or asks how much of a named protein is on a "
        "fragment."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "name":          {"type": "string"},
            "attachment_id": {"type": "string"},
        },
        "required": ["name"],
    },
}

# --- Anthropic server-side web search (fallback) -------------------------

WEB_SEARCH_TOOL = {
    "type": "web_search_20250305",
    "name": "web_search",
    "max_uses": 5,
}

AIPLASMIDDESIGN_TOOLS = [
    ANNOTATE_ATTACHMENT_TOOL,
    SIMULATE_ASSEMBLY_TOOL,
    DIGEST_PLASMID_TOOL,
    FIND_PRIMER_BINDING_TOOL,
    SCORE_SANGER_PRIMER_TOOL,
    ANALYZE_DESIGN_INTENT_TOOL,
    VERIFY_ASSEMBLY_TOOL,
    COMPARE_TO_CHOICE_TOOL,
    GOLDEN_GATE_ASSEMBLE_TOOL,
    ROUTE_WORKFLOW_TOOL,
    LOOKUP_KB_PART_TOOL,
]

ALL_TOOLS = AIPLASMIDDESIGN_TOOLS + [WEB_SEARCH_TOOL]

SYSTEM_PROMPT = """You are a molecular cloning expert agent backed by AIPlasmidDesign's
deterministic toolchain.

CRITICAL RULES:
1. You will NEVER see raw DNA. Plasmids are referenced by `attachment_id`.
   The tools resolve sequences server-side and return high-level properties
   only.
2. NEVER quote, transcribe, or reproduce DNA. For DNA-valued multiple-choice
   options (shown as '[redacted DNA, N nt]'), use compare_to_choice to
   match an attachment against each option letter.
3. ALWAYS call at least one AIPlasmidDesign tool before answering any
   question with attachments.
4. Use web_search ONLY for external context (Addgene metadata, vendor specs,
   citations).
5. POST-ASSEMBLY GEOMETRY for Sanger questions: when the question asks
   about a primer for sequencing AFTER a cloning step (oligo annealing,
   Golden Gate insert, Gibson assembly), you MUST first simulate that
   assembly (golden_gate_assemble for Type-IIS Golden Gate, simulate_assembly
   otherwise) to register a product attachment, then call score_sanger_primer
   on the PRODUCT attachment, not the pre-assembly template.
6. METHOD/ENZYME ENUMERATION: when the multiple-choice options are cloning
   methods (Gibson / Golden Gate / Gateway / restriction / SDM / USER) or
   restriction enzymes (BsmBI / BsaI / BbsI / SapI / etc.), you MUST call
   route_workflow once (for method choice) OR digest_plasmid once per
   candidate enzyme (for enzyme choice) BEFORE answering. Do not pick a
   method or enzyme from general biology heuristics. The unified-predesign
   router is the source of truth on method feasibility.

TOOL-SELECTION GUIDE:
- "What is on this plasmid / where is feature X?"     → annotate_attachment
- "Is this design correct / what's wrong with it?"    → verify_assembly
- "What's the purpose of this cloning request?"       → analyze_design_intent
- "Which primer for Sanger sequencing across X?"      → score_sanger_primer
  (call annotate_attachment first if you need the position of X)
- "Which fragment lengths confirm a correct clone?"   → digest_plasmid
- "What does this assembly produce?"                  → golden_gate_assemble for Type-IIS Golden Gate;
                                                          simulate_assembly for Gibson/restriction/Gateway/SDM
- "Which choice matches the assembly product?"        → compare_to_choice
  (after simulate_assembly; call once per option letter)
- "Where do these primers bind on this template?"     → find_primer_binding_sites
- "Which cloning method/approach should I use?"       → route_workflow (REQUIRED before answering)
- "Which enzyme(s) for this Golden Gate?"             → digest_plasmid for each candidate (REQUIRED)
- "How much of protein X is on this fragment?"        → lookup_kb_part(name=X, attachment_id=...)
- annotate_attachment did not return feature X        → lookup_kb_part(name=X) — KB has parts the annotator may miss

ANSWER FORMAT:
- 1-3 sentences of reasoning grounded in tool output.
- End with a single line: "ANSWER: <letter>" (e.g. "ANSWER: B").
"""
