"""Interpreter agent — Sonnet-4.6 ReAct loop over the deterministic
PlasmidRegistry tools. Answers natural-language questions about an
uploaded plasmid or an inventory of plasmids.

Usage:
    result = await run_interpreter(
        question="Where is the gRNA cloning cassette?",
        registry=registry,
    )
    print(result.answer)
    for c in result.citations:
        print(c)

Mocking pattern for tests: pass a fake `client` whose `messages.create`
returns a synthetic response; pass a fake `dispatch_fn` to short-
circuit registry method calls.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

from agent_v2.interpreter.plasmid_registry import PlasmidRegistry
from agent_v2.interpreter.tools import INTERPRETER_TOOLS, dispatch_interpreter_tool


SYSTEM_PROMPT = """You are the DNA Annotation Interpreter — a deterministic reader of plasmid annotation envelopes. You answer questions about modules, features, interactions, and per-amino-acid positions by calling the registry tools. You do NOT speculate beyond what the tools return.

You have an inventory of one or more uploaded plasmids. Every tool accepts an optional `plasmid_id`; omit it to fan out across the whole inventory, supply it to scope to one plasmid.

Discipline:
1. If the user names a specific plasmid (filename, gene name, "the lentivirus with EGFP"), call `resolve_plasmid` FIRST. If `ok=false`, tell the user no registered plasmid matches and stop — do not invent one.
2. For "Where is X?" questions, use `find_modules` (modules), `find_features` (annotations), or `find_cloning_features` (restriction / Gateway). Return start–end coordinates and a one-line description of WHY the rule fired (rule_id, submodules, golden_gate enzyme, etc.).
3. For "What is the Nth amino acid in Y?" use `lookup_amino_acid` with feature_name=Y, aa_index=N. Report all four fields: residue letter + full name, codon DNA, position in feature, position in parent ORF.
4. For "What promoter / polyA drives Y?" use `expression_cassette_for(cds_name=Y)`. Report both regulatory submodules with their coordinates.
5. For "What is this plasmid for?" use `infer_application`. If the confidence is low, say so and list the modules that ARE present.
6. When NO plasmid in the inventory has what was asked for and the question NAMES a specific part (e.g. "lentiCRISPR v2", "Tet-On 3G", "pSpCas9(BB)-2A-GFP", "pegRNA backbone with tevopreQ1"), call `find_external_part(description=...)`. This now AUTO-DOWNLOADS each candidate's GenBank from Addgene, annotates it through the standard pipeline, and registers it in the local registry with a plasmid_id like `addgene_174038`. After the call, you can call any of the local-lookup tools (find_features, find_modules, lookup_amino_acid, expression_cassette_for, infer_application) with that plasmid_id to answer follow-up questions about the candidate in the same turn — no need to ask the user to upload anything when the download succeeded. When the download fails (rare), fall back to surfacing the sequences-page URL and asking the user to upload.
7. When the inventory is empty AND the question is generic ("Where is X?" without naming a particular plasmid), say so explicitly: "None of the 3 uploaded plasmids contain a guide expression cassette." Do not paper over it. Do NOT trigger find_external_part for a generic absence — only when the user named a specific part.
8. The registry may be EMPTY (the user asked about a plasmid without uploading one). In that case skip the local lookups and go straight to `find_external_part` if the question names a specific part; otherwise say something like "I don't have any plasmids loaded — drop a `.gb` file to inspect features, or name a specific plasmid (e.g. 'lentiCRISPR v2') and I'll search Addgene."
9. WORKFLOW REDIRECT: if the user is actually asking to BUILD / DESIGN / ASSEMBLE / MODIFY a plasmid (verbs: build, design, clone, assemble, insert, mutate, knock out, introduce, replace, swap, optimize, generate a vector), call `redirect_to_workflow(target_intent="PLASMID_CLONING", reason=...)` IMMEDIATELY — do NOT run the QA tools first. The cloning pipeline will spin up the design subagents instead. Same for CRISPR guide design (target_intent="CRISPR_GUIDE"): if the user wants sgRNA design / off-target scoring / PAM choice / pegRNA design, redirect. Bias toward redirecting when in doubt — completing the wrong workflow wastes 30-60 s of pipeline time before the user discovers the mismatch.

Output format: a concise plain-language answer (Markdown OK). Always include the relevant coordinates (start–end, in bp) when describing a module or feature. When fanning out across the inventory, group your answer by plasmid (only mention plasmids that had results unless the user asked about all of them).

Hard limits:
- Never invent feature names, coordinates, or AA positions. Every claim must come from a tool result.
- Never return raw DNA past a single codon (which `lookup_amino_acid` already includes).
- ≤ 350 words.
"""


@dataclass
class InterpreterResult:
    answer: str
    citations: list[dict[str, Any]] = field(default_factory=list)
    trace: list[dict[str, Any]] = field(default_factory=list)
    n_tool_calls: int = 0
    redirect: Any = None        # WorkflowRedirect or None


def _block_to_dict(block: Any) -> dict[str, Any]:
    if isinstance(block, dict):
        return block
    out: dict[str, Any] = {"type": getattr(block, "type", "")}
    if out["type"] == "text":
        out["text"] = getattr(block, "text", "")
    elif out["type"] == "tool_use":
        out["id"] = getattr(block, "id", "")
        out["name"] = getattr(block, "name", "")
        out["input"] = getattr(block, "input", {}) or {}
    return out


async def run_interpreter(
    question: str,
    registry: PlasmidRegistry,
    *,
    client: Any = None,
    dispatch_fn: Optional[Callable[[str, dict[str, Any], PlasmidRegistry], dict[str, Any]]] = None,
    on_tool_event: Optional[Callable[[dict[str, Any]], Any]] = None,
    model: Optional[str] = None,
    max_iters: int = 8,
) -> InterpreterResult:
    """Run the interpreter ReAct loop. Returns the final answer plus
    a trace of every tool call + result. Empty registry produces a
    canned response without consulting the LLM."""

    if client is None:
        import anthropic
        client = anthropic.Anthropic()
    if dispatch_fn is None:
        dispatch_fn = dispatch_interpreter_tool
    model = model or os.getenv("AGENT_MODEL", "claude-sonnet-4-6")

    # Initial registry snapshot so the agent knows what it has access to.
    inventory_preview = {
        "n_plasmids": registry.n(),
        "plasmids": [
            {"plasmid_id": p.plasmid_id, "name": p.name,
              "length_bp": len(p.sequence),
              "module_types": sorted({m.get("module_type") for m in p.modules() if m.get("module_type")})[:10]}
            for p in registry.all()
        ],
    }

    messages: list[dict[str, Any]] = [
        {
            "role": "user",
            "content": (
                f"Registry contents (preview):\n{json.dumps(inventory_preview, indent=2)}\n\n"
                f"User question:\n{question}"
            ),
        },
    ]

    trace: list[dict[str, Any]] = []
    citations: list[dict[str, Any]] = []
    n_tool_calls = 0
    final_text = ""

    # Prompt-caching prefix: system + tools + first user message.
    cached_system = [{"type": "text", "text": SYSTEM_PROMPT,
                       "cache_control": {"type": "ephemeral"}}]
    cached_tools = [{**t} for t in INTERPRETER_TOOLS]
    if cached_tools:
        cached_tools[-1] = {**cached_tools[-1],
                              "cache_control": {"type": "ephemeral"}}
    if messages and messages[0]["role"] == "user" and isinstance(messages[0].get("content"), str):
        messages[0] = {
            "role": "user",
            "content": [{"type": "text", "text": messages[0]["content"],
                          "cache_control": {"type": "ephemeral"}}],
        }

    for _ in range(max_iters):
        resp = client.messages.create(
            model=model,
            max_tokens=2048,
            system=cached_system,
            tools=cached_tools,
            messages=messages,
        )

        # Normalise the response so test mocks can pass dicts or SDK objects.
        blocks = [_block_to_dict(b) for b in (resp.content if hasattr(resp, "content") else resp.get("content", []))]
        stop_reason = getattr(resp, "stop_reason", None) or (resp.get("stop_reason") if isinstance(resp, dict) else None)

        # Always extract any text the assistant emitted on this turn.
        text_now = "".join(b.get("text", "") for b in blocks if b.get("type") == "text")
        if text_now.strip():
            final_text = text_now.strip()

        tool_uses = [b for b in blocks if b.get("type") == "tool_use"]

        # Echo the assistant turn back to the conversation.
        messages.append({"role": "assistant", "content": [
            ({"type": "text", "text": b["text"]} if b.get("type") == "text"
             else {"type": "tool_use", "id": b["id"], "name": b["name"], "input": b.get("input", {})})
            for b in blocks
        ]})

        if not tool_uses:
            break  # end_turn

        # Dispatch each tool and append a single tool_result message.
        # Each call is wrapped in timing + start/end event hooks so the
        # parent SSE stream can show live tool-call progress.
        import time as _time_mod
        tool_result_blocks = []
        for tu in tool_uses:
            iteration_idx = n_tool_calls
            n_tool_calls += 1
            t0 = _time_mod.monotonic()
            if on_tool_event is not None:
                try:
                    res_evt = on_tool_event({
                        "phase": "start",
                        "iteration": iteration_idx,
                        "tool": tu["name"],
                        "input": tu.get("input") or {},
                        "t_unix": _time_mod.time(),
                    })
                    if hasattr(res_evt, "__await__"):
                        await res_evt
                except Exception:
                    pass
            try:
                result = dispatch_fn(tu["name"], tu.get("input") or {}, registry)
            except Exception as e:
                result = {"ok": False, "error": f"{type(e).__name__}: {e}"}
            # Redirect short-circuit: if the agent called redirect_to_workflow,
            # stop the loop and surface the redirect to the orchestrator.
            if isinstance(result, dict) and isinstance(result.get("_redirect"), dict):
                from agent_v2.redirect import WorkflowRedirect
                redirect = WorkflowRedirect(
                    target_intent=result["_redirect"].get("target_intent", "PLASMID_CLONING"),
                    reason=result["_redirect"].get("reason", ""),
                )
            elapsed_ms = int((_time_mod.monotonic() - t0) * 1000)
            n_results = (len(result.get("results", []))
                         if isinstance(result, dict) and isinstance(result.get("results"), list)
                         else (result.get("n_candidates") if isinstance(result, dict) else None))
            trace.append({
                "tool": tu["name"],
                "input": tu.get("input"),
                "result_keys": list(result.keys()) if isinstance(result, dict) else None,
                "n_results": len(result.get("results", [])) if isinstance(result, dict) else None,
                "elapsed_ms": elapsed_ms,
            })
            if on_tool_event is not None:
                try:
                    res_evt = on_tool_event({
                        "phase": "end",
                        "iteration": iteration_idx,
                        "tool": tu["name"],
                        "elapsed_ms": elapsed_ms,
                        "n_results": n_results,
                        "ok": (result.get("ok") if isinstance(result, dict) else None),
                        "t_unix": _time_mod.time(),
                    })
                    if hasattr(res_evt, "__await__"):
                        await res_evt
                except Exception:
                    pass
            for r in (result.get("results") or [] if isinstance(result, dict) else []):
                if isinstance(r, dict) and r.get("plasmid_id"):
                    citations.append({
                        "plasmid_id": r["plasmid_id"],
                        "via": tu["name"],
                        "ref": _summarise_citation(tu["name"], r),
                    })
            tool_result_blocks.append({
                "type": "tool_result",
                "tool_use_id": tu["id"],
                "content": json.dumps(result),
            })
            if isinstance(result, dict) and isinstance(result.get("_redirect"), dict):
                # Return immediately with the redirect attached.
                return InterpreterResult(
                    answer=(f"Redirecting to {result['_redirect'].get('target_intent', '?')} "
                            f"workflow: {result['_redirect'].get('reason', '')}"),
                    citations=citations, trace=trace,
                    n_tool_calls=n_tool_calls,
                    redirect=redirect,
                )
        messages.append({"role": "user", "content": tool_result_blocks})

        if stop_reason == "end_turn":
            break

    return InterpreterResult(
        answer=final_text or "(no answer produced)",
        citations=citations,
        trace=trace,
        n_tool_calls=n_tool_calls,
    )


def _summarise_citation(tool: str, r: dict[str, Any]) -> str:
    if tool == "find_modules":
        return f"{r.get('module_type')} {r.get('start')}-{r.get('end')} ({r.get('rule_id')})"
    if tool == "find_features":
        return f"{r.get('name')} {r.get('start')}-{r.get('end')}"
    if tool == "lookup_amino_acid":
        return f"{r.get('feature_name')} aa {r.get('aa_position_in_feature')}={r.get('letter')}"
    if tool == "expression_cassette_for":
        return f"cassette for {r.get('cds_name')} promoter={(r.get('promoter') or {}).get('name')}"
    if tool == "infer_application":
        return r.get("application", "?")
    if tool == "plasmid_summary":
        return f"summary {r.get('name')} ({r.get('n_modules')} modules)"
    if tool == "find_cloning_features":
        return f"{r.get('subtype') or r.get('feature_family')} {r.get('start')}-{r.get('end')}"
    return ""
