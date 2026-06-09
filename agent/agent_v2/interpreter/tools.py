"""Anthropic tool schemas + handlers for the DNA-annotation
interpreter. Each tool wraps a PlasmidRegistry method, takes an
optional plasmid_id (omit to fan out across the inventory), and
returns a JSON-serialisable dict that flows straight back into the
agent's tool-result message.
"""
from __future__ import annotations

from typing import Any, Optional

from agent_v2.interpreter.plasmid_registry import PlasmidRegistry


# ──────────────────────────────────────────────────────────────────
# Tool schemas (Anthropic spec)
# ──────────────────────────────────────────────────────────────────
INTERPRETER_TOOLS: list[dict[str, Any]] = [
    {
        "name": "resolve_plasmid",
        "description": (
            "Map a user reference (a filename, a gene name, or a description like "
            "'the lentivirus with EGFP and PuroR') to one or more concrete plasmid_ids "
            "in the registry. Call this FIRST whenever the user names a specific "
            "plasmid. Returns {ok, method, matches: [{plasmid_id, name, score, reason}]} "
            "or ok=false when nothing matches."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {"type": "string", "description": "Free-text plasmid reference."},
            },
            "required": ["query"],
        },
    },
    {
        "name": "plasmid_summary",
        "description": (
            "Get a high-level overview of one plasmid (or every plasmid in the "
            "registry if plasmid_id is omitted). Use this for orientation before "
            "drilling in."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "plasmid_id": {"type": "string",
                                "description": "Omit to fan out across the inventory."},
            },
        },
    },
    {
        "name": "find_modules",
        "description": (
            "Search modules by module_type, rule_id, or name (case-insensitive "
            "substring). Returns coords, submodules, and metadata. Examples: "
            "'guide_expression_cassette', 'sgrna', 'golden_gate', "
            "'lac_alpha_disrupted', 'POL3-GG-01', 'lentiviral_payload'."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {"type": "string"},
                "plasmid_id": {"type": "string", "description": "Omit to fan out."},
            },
            "required": ["query"],
        },
    },
    {
        "name": "find_features",
        "description": (
            "Search pLannotate feature annotations by name / sseqid / gene_name / "
            "protein_name. Returns coords + KB metadata (gene_name, protein_name, "
            "organism)."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {"type": "string"},
                "plasmid_id": {"type": "string", "description": "Omit to fan out."},
            },
            "required": ["query"],
        },
    },
    {
        "name": "find_cloning_features",
        "description": (
            "Search cloning-feature annotations: Type II / Type IIs restriction "
            "sites, Gateway att sites, PCR design warnings. Query matches name, "
            "subtype, or feature_family."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {"type": "string"},
                "plasmid_id": {"type": "string", "description": "Omit to fan out."},
            },
            "required": ["query"],
        },
    },
    {
        "name": "lookup_amino_acid",
        "description": (
            "Resolve a 1-based amino-acid position within a feature (or the parent "
            "ORF directly). Returns the residue letter, full name, codon DNA "
            "sequence, the AA's position within the feature and within the parent "
            "ORF, and the ORF coordinates. Use this for 'What is the Nth amino "
            "acid in Cas9?' style questions."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "feature_name": {"type": "string",
                                  "description": "Sub-feature name (e.g. 'Cas9', 'PuroR') "
                                                 "or an ORF name."},
                "aa_index": {"type": "integer", "minimum": 1,
                              "description": "1-based AA position WITHIN the feature."},
                "plasmid_id": {"type": "string", "description": "Omit to fan out."},
            },
            "required": ["feature_name", "aa_index"],
        },
    },
    {
        "name": "expression_cassette_for",
        "description": (
            "Find the expression cassette (upstream regulatory promoter + CDS + "
            "downstream regulatory polyA) that drives the given CDS. Walks the "
            "interaction_builder output first, then falls back to scanning "
            "mammalian_pol2_expression_cassette modules. Returns promoter + polyA "
            "names and coordinates."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "cds_name": {"type": "string",
                              "description": "CDS feature name (e.g. 'Cas9', 'EGFP')."},
                "plasmid_id": {"type": "string", "description": "Omit to fan out."},
            },
            "required": ["cds_name"],
        },
    },
    {
        "name": "infer_application",
        "description": (
            "Pattern-match the module composition of a plasmid against a small "
            "rule table to suggest its most likely application (e.g. 'CRISPR-Cas9 "
            "lentiviral knockout vector', 'Blue/white cloning vector', 'Mammalian "
            "expression vector'). Use this for 'What is this plasmid used for?'."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "plasmid_id": {"type": "string",
                                "description": "Omit to fan out across the inventory."},
            },
        },
    },
    {
        "name": "find_external_part",
        "description": (
            "ESCALATION: when the local registry has no match for a specific "
            "part the user named (a promoter, sgRNA backbone, fluorescent "
            "reporter, etc.), search Addgene and AUTO-DOWNLOAD + ANNOTATE the "
            "top candidate's GenBank file so subsequent tool calls "
            "(find_features, find_modules, lookup_amino_acid, "
            "expression_cassette_for, infer_application) can answer detail "
            "questions about it without a second turn. The result includes "
            "the candidate's name, Addgene ID, URL, description, depositor, "
            "paper title + PMID + DOI, plus a registered_as: plasmid_id you "
            "can pass straight back into the other tools. Use ONLY after the "
            "local tools have come up empty AND the question names a "
            "specific part — never as a shotgun search."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "description": {
                    "type": "string",
                    "description": (
                        "A short description of the part the user is asking "
                        "about (e.g. 'lentiCRISPR v2', 'Tet-On 3G inducible "
                        "promoter', 'pSpCas9(BB)-2A-GFP')."
                    ),
                },
                "max_results": {
                    "type": "integer",
                    "description": "How many top candidates to return (1..5).",
                    "default": 3,
                },
            },
            "required": ["description"],
        },
    },
   {
        "name": "redirect_to_workflow",
        "description": (
            "ESCAPE HATCH: if the user's request is fundamentally better "
            "served by a different workflow, call this to redirect the "
            "orchestrator instead of producing a half-useful answer. Use it "
            "sparingly — only when the local tools cannot answer at all. "
            "Targets:\n"
            "  - PLASMID_CLONING: the user wants to BUILD / DESIGN / ASSEMBLE / "
            "MODIFY a plasmid (verbs: build, design, clone, assemble, insert, "
            "mutate). Cannot be satisfied by lookup alone.\n"
            "  - CRISPR_GUIDE: the user wants sgRNA design / scoring / PAM "
            "choice. Lookup of an existing plasmid is NOT this.\n"
            "  - REJECT: outside molecular biology / cloning / guide design.\n"
            "Reason should be one short sentence explaining WHY the current "
            "workflow can't satisfy the request."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target_intent": {
                    "type": "string",
                    "enum": ["PLASMID_CLONING", "CRISPR_GUIDE", "REJECT"],
                },
                "reason": {"type": "string"},
            },
            "required": ["target_intent", "reason"],
        },
    },
]


# ──────────────────────────────────────────────────────────────────
# Dispatch
# ──────────────────────────────────────────────────────────────────
def dispatch_interpreter_tool(name: str, args: dict[str, Any], registry: PlasmidRegistry) -> dict[str, Any]:
    """Route a tool name to its registry method. Returns a normalised
    result dict. Unknown tool names get a {"ok": False, "error": ...}.
    """
    if name == "resolve_plasmid":
        return registry.resolve_plasmid(args.get("query", ""))

    if name == "plasmid_summary":
        return registry.fan_out("summary", plasmid_id=args.get("plasmid_id"))

    if name == "find_modules":
        return registry.fan_out("find_modules",
                                 plasmid_id=args.get("plasmid_id"),
                                 query=args.get("query", ""))

    if name == "find_features":
        return registry.fan_out("find_features",
                                 plasmid_id=args.get("plasmid_id"),
                                 query=args.get("query", ""))

    if name == "find_cloning_features":
        return registry.fan_out("find_cloning_features",
                                 plasmid_id=args.get("plasmid_id"),
                                 query=args.get("query", ""))

    if name == "lookup_amino_acid":
        return registry.fan_out("lookup_amino_acid",
                                 plasmid_id=args.get("plasmid_id"),
                                 feature_name=args.get("feature_name", ""),
                                 aa_index=int(args.get("aa_index", 0)))

    if name == "expression_cassette_for":
        return registry.fan_out("expression_cassette_for",
                                 plasmid_id=args.get("plasmid_id"),
                                 cds_name=args.get("cds_name", ""))

    if name == "infer_application":
        return registry.fan_out("infer_application",
                                 plasmid_id=args.get("plasmid_id"))

    if name == "find_external_part":
        return _dispatch_find_external_part(args, registry)

    if name == "redirect_to_workflow":
        return {
            "ok": True,
            "_redirect": {
                "target_intent": args.get("target_intent"),
                "reason": args.get("reason", ""),
            },
        }

    return {"ok": False, "error": f"Unknown interpreter tool: {name}"}


def _dispatch_find_external_part(args: dict[str, Any], registry: "PlasmidRegistry") -> dict[str, Any]:
    """Run the Addgene search + auto-download + annotate + register
    synchronously from the interpreter's tool-dispatch context."""
    import asyncio
    from splicify_api.external_search import (
        search_addgene, fetch_addgene_entry, download_addgene_gb,
    )
    from splicify_api.annotation_cache import annotate_llm_cached
    from splicify_api.agent.agent_tools import extract_seq_from_genbank

    description = (args.get("description") or "").strip()
    if not description:
        return {"ok": False, "error": "description is required"}
    MAX_PER_CALL = 2
    MAX_PER_SESSION = 6
    max_results = max(1, min(MAX_PER_CALL, int(args.get("max_results") or 1)))
    existing = int(getattr(registry, "_external_lookups_done", 0) or 0)
    remaining = MAX_PER_SESSION - existing
    if remaining <= 0:
        return {"ok": False, "candidates": [], "n_candidates": 0,
                "error": (f"External-lookup budget exhausted "
                            f"({MAX_PER_SESSION} Addgene fetches); "
                            "ask the user to upload directly.")}
    max_results = min(max_results, remaining)
    setattr(registry, "_external_lookups_done", existing + max_results)

    async def _go():
        candidates = await search_addgene(description, max_results=max_results)
        out = []
        for cand in candidates:
            entry = await fetch_addgene_entry(cand.addgene_id)
            d = cand.to_dict()
            registered_pid: str | None = None
            if entry:
                d.update({
                    "description": entry.description,
                    "depositor": entry.depositor,
                    "pmid": entry.pmid,
                    "doi": entry.doi,
                    "paper_title": entry.paper_title,
                    "sequences_page": f"https://www.addgene.org/{cand.addgene_id}/sequences/",
                })
            # Download + annotate + register. If anything fails we
            # still return the metadata + a download_note so the agent
            # can fall back to asking the user to upload manually.
            try:
                gb_text = await download_addgene_gb(cand.addgene_id)
                if gb_text:
                    seq = extract_seq_from_genbank(gb_text)
                    if seq and len(seq) >= 50:
                        envelope = await annotate_llm_cached(seq, circular=True)
                        envelope.setdefault("sequence", seq)
                        pid = f"addgene_{cand.addgene_id}"
                        registry.register(
                            pid, envelope,
                            name=entry.name if entry else cand.name,
                        )
                        registered_pid = pid
                        d["length_bp"] = len(seq)
                        idx = registry.get(pid)
                        if idx is not None:
                            summary = idx.summary()
                            d["module_types"] = summary.get("module_types", [])[:10]
                            d["n_modules"] = summary.get("n_modules", 0)
                            d["n_annotations"] = summary.get("n_annotations", 0)
            except Exception as e:
                d["download_error"] = f"{type(e).__name__}: {e}"
            d["registered_as"] = registered_pid
            if not registered_pid:
                d["download_note"] = (
                    "Could not auto-download the .gb. Ask the user to "
                    "upload it from the sequences page if detailed feature "
                    "/ AA / module lookups are needed."
                )
            out.append(d)
        return out

    try:
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            loop = None
        if loop is None:
            results = asyncio.run(_go())
        else:
            import concurrent.futures
            with concurrent.futures.ThreadPoolExecutor(max_workers=1) as ex:
                fut = ex.submit(asyncio.run, _go())
                results = fut.result(timeout=120)
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}

    return {
        "ok": True,
        "source": "addgene",
        "query": description,
        "n_candidates": len(results),
        "n_registered": sum(1 for r in results if r.get("registered_as")),
        "candidates": results,
    }

