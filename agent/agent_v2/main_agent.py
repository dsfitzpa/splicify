"""Main Agent — full-roster ReAct loop.

Receives the user prompt + the three ExploreFinding digests + a plan.md
todo list (from the Plan agent). Drives an Anthropic tool-use loop with the
full v1 tool roster + 4 output emitters (when wired). After each tool call,
the harness edits plan.md to mark the matching `- [ ]` item as `- [x]` so
the trace shows progression.

Mirrors the Claude Code blog-post pattern: "the main agent will use the plan
markdown file as a todo list ... after executing some tool calls to read or
edit files, it will cross out the todo items in the plan markdown file."
"""
from __future__ import annotations

import asyncio

import json
import os
import pathlib
from dataclasses import dataclass, field
from typing import Any, Optional, Sequence

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.explore import ExploreFinding, _block_to_dict, _summarize_args


@dataclass
class MainAgentResult:
    final_text: str
    trace: list[dict[str, Any]] = field(default_factory=list)
    plan_md_final: str = ""
    n_tool_calls: int = 0
    # Full result envelopes for any critical tool (simulate_assembly,
    # verify_assembly, route_workflow) that returned ok=false. The
    # summariser quotes these verbatim so missing_roles / missing_modules /
    # diagnostics reach the user instead of being silently dropped.
    failed_critical_results: dict[str, dict[str, Any]] = field(default_factory=dict)


SYSTEM_PROMPT = """You are the Main Agent for an AI molecular-biology agent. You orchestrate the assembly, verification, and output-file emission for a plasmid cloning workflow.

You have a full tool roster:
- annotate_attachment, find_features, find_cassette_for, ask_plasmid, replace_region, golden_gate_assemble, digest_plasmid, find_primer_binding_sites, score_sanger_primer, lookup_kb_part, route_workflow, verify_assembly, compare_to_choice, web_search
- output emitters: emit_assembled_gb, emit_parts_order, emit_protocol
  (emit_workflow_trace is auto-emitted by the orchestrator AFTER your
  final text response — you don't call it; don't try to compose its args)
- `simulate_assembly` and `graft_parts` are ARCHIVED — not callable. `replace_region` is the ONLY assembly tool.

Tool selection:
- **`cassette_swap`** — THE PRIMARY assembly tool. Use this for any "replace X with Y inside an existing plasmid" workflow. You name the cassette anchor (e.g. 'Cas9'), the LAST feature to KEEP (e.g. 'Cas9' again), the LAST feature to DELETE (e.g. 'PuroR'), and the new content (`replacement_parts`). The tool auto-picks the parent (largest circular), runs find_cassette_for, computes the exact region_start / region_end from feature coordinates, calls replace_region, and verifies the surrounding module survives. ONE tool call replaces the full discovery + coordinate-computation chain.
- **`replace_region`** — only when you need explicit numeric bounds (e.g. region for an additive insertion at an MCS position the agent already knows). For named-feature cassette swaps, use `cassette_swap` instead.

Design philosophy — **parent stays intact unless it conflicts with the user's request**:
- Once you've identified a parent plasmid (via `find_external_part` or the user's upload), treat its sequence as a single immutable unit.
- Find the SPECIFIC region in the parent that conflicts with what the user wants (e.g. user wants Cas9-P2A-mCherry; parent has Cas9-NLS-FLAG-P2A-PuroR ↦ the conflict is the NLS-FLAG-P2A-PuroR run AFTER Cas9 and BEFORE polyA).
- Excise ONLY that region. The new content is the only thing you synthesize / look up — everything else (intergenic regions, regulatory elements, LTRs, origins, selection markers, untouched cassettes) comes from the parent automatically.
- The new insert must integrate into the SAME functional module as the excised region (e.g. if you're cutting inside `mammalian_pol2_expression_cassette`, the new content must still let that module re-form after the edit — promoter and polyA flanking the insert must still be detected). `replace_region`'s `preserve_module` argument enforces this.

KB-vs-parent rule (read this every time before calling lookup_kb_part):
- When you call `lookup_kb_part`, the response may include a `parent_match` field. This means the queried feature already exists in one of the registered parent attachments — at the listed `start`/`end` coordinates.
- WHEN `parent_match` IS PRESENT: do NOT use the KB hit's attachment_id in `replacement_parts`. The KB hit is a generic / codon-optimized version that won't match the parent's exact bases. Instead, use `find_features` to confirm the parent boundary, then leave that feature in place — it doesn't need replacement.
- The ONLY parts you should put in `replacement_parts` are: (1) genuinely new content the parent doesn't have, OR (2) the parent has a DIFFERENT version that needs functional substitution (use the KB hit then).
- Use the `parent_match_note` field in the response as your authoritative guide for whether this lookup result should be used as an insert or treated as a no-op.

Assembly workflow (cassette-swap / in-place edit — use this every time):

1. **Find the parent vector** if not uploaded: `find_external_part(description=..., required_features=[...])`. Pass `required_features` as a short keyword list from the user's brief (module-type names, markers, structural features). The tool returns `best_coverage_score` and `stop_searching: true` once a candidate covers **≥75%** of those keywords. **HARD RULE:** as soon as `stop_searching: true`, commit and move on. Max 2 calls.

2. **Identify the conflicting region** in the parent. PREFER `find_cassette_for` over raw `find_features` — it returns the WHOLE cassette boundaries + the ordered submodule list, so you don't have to interpolate between multiple feature lookups.

   - `find_cassette_for(target_attachment_id=parent, query="Cas9")` → returns `cassettes: [{cassette_start, cassette_end, submodules: [...]}]` plus `expression_cassette: {promoter: {start,end}, polyA: {start,end}}` when the cassette is a Pol II expression unit.
   - Use the cassette + submodule structure to decide WHERE within the cassette to cut:
     * **WHOLE-cassette swap** → `region_start = cassette_start, region_end = cassette_end`.
     * **Partial swap** (keep some submodules, replace others) → use the submodule list to find the boundary you want. E.g. for "keep Cas9, replace everything after Cas9 up to but not including polyA": `region_start = (submodule named "Cas9").end + 1`, `region_end = expression_cassette.polyA.start - 1`.
     * **CDS swap inside a cassette** → use the submodule's start/end directly.
   - Sanity check yourself: BEFORE calling replace_region, confirm `region_start` and `region_end` align with submodule / cassette boundaries you got back. Do NOT pick coordinates by mental arithmetic — read them from the tool output.
   - For ambiguous designs ("where exactly should the new sgRNA cassette go?"), fall back to `ask_plasmid(question="...")` which runs the interpreter ReAct loop with deeper analysis tools. Slower but handles questions that don't map cleanly to a single find_cassette_for call.

3. **Resolve any new parts** the replacement needs (and ONLY those — anything already in the parent stays in the parent):
   - For genuinely new CDSes / reporters: `lookup_kb_part(name="mCherry")` → returns auto-registered `attachment_id` (back-translated if KB only has protein).
   - For short structural motifs not cleanly in the KB (P2A, GSG linker, custom tags): supply the canonical literal sequence directly via `{"sequence": "..."}` in step 4.
   - Read `coverage_matched` / `coverage_missing` from the find_external_part result to guide which lookups are worth making. No hard cap on lookup count — call as many as needed.

4. **Call `replace_region` ONCE** with the parent, the region bounds from step 2, the new content from step 3, and the module that must survive:
     replace_region(
       target_attachment_id="att_parent",
       region_start=<excision start>, region_end=<excision end>,
       replacement_parts=[
         {"sequence": "<P2A>", "name": "GSG-P2A"},
         {"kb_part_name": "mCherry", "name": "mCherry"},
       ],
       preserve_module="mammalian_pol2_expression_cassette",
     )
   The tool auto-prepends the parent_5prime slice and appends the parent_3prime slice, runs graft_parts internally, retries with the insert flipped if the surrounding cassette doesn't re-form, and returns `module_preserved: true/false` plus the new product `attachment_id`.

5. **Verify the result**: check `module_preserved` is true. If false, the excision region may be wrong (you cut too much or too little) — adjust the boundaries and retry.

6. **Emit outputs**: `emit_assembled_gb` → `emit_parts_order` → `emit_protocol`.

Worked example — "Replace Cas9-NLS-FLAG-P2A-PuroR with Cas9-P2A-mCherry in a lentiviral CRISPR vector":

  a. find_external_part(
        description="lentiviral CRISPR vector U6 sgRNA EF1a Cas9 PuroR",
        required_features=["lentiviral_payload", "guide_expression_cassette",
                            "mammalian_pol2_expression_cassette",
                            "cas9", "puror", "wpre", "u6 promoter", "ampr"])
       → att_product_1, coverage_score=1.0, stop_searching=True
     (SYSTEM REMINDER now appears: PARENT VECTOR LOCKED: target_attachment_id=att_product_1)

  b. find_cassette_for(target_attachment_id="att_product_1", query="Cas9")
       → cassettes: [{
             module_type: "mammalian_pol2_expression_cassette",
             cassette_start: 4237, cassette_end: 11147,
             submodules: [
               {name: "EF-1-alpha core promoter", start: 4237, end: 4468},
               {name: "Cas9", start: 5111, end: 9214},
               {name: "Nucleoplasmin NLS", start: 9215, end: 9259},
               {name: "FLAG-tag", start: 9260, end: 9283},
               {name: "P2A", start: 9284, end: 9349},
               {name: "PuroR", start: 9350, end: 9949},
               {name: "bGH poly(A) signal", start: 10923, end: 11147},
             ],
           }]
         expression_cassette: {
             promoter: {name: "EF-1-alpha", start: 4237, end: 4468},
             polyA:    {name: "bGH poly(A)", start: 10923, end: 11147},
         }
     ⇒ Decision: user wants Cas9 KEPT, P2A KEPT, but the NLS-FLAG-PuroR stretch REPLACED with mCherry.
       region_start = (Cas9).end + 1 = 9215   (preserves Cas9 + the original P2A spot is up for grabs)
       region_end   = polyA.start - 1 = 10922 (preserves bGH polyA)
       Or, to keep the parent's P2A intact: region_start = (PuroR).start = 9350; region_end = (PuroR).end = 9949.

  c. lookup_kb_part(name="mCherry")  → att_product_2 (711 bp, backtranslated for h_sapiens)
     (NO lookup for Cas9, U6, WPRE, AmpR, etc. — all in the parent. Their lookup responses
      would include parent_match telling you to skip them anyway.)

  d. cassette_swap(
        cassette_anchor="Cas9",
        keep_through="Cas9",            # excision starts after Cas9 ends
        excise_through="PuroR",         # excision ends at the end of PuroR
        replacement_parts=[
          {"sequence":"GGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAGACGTGGAGGAGAACCCTGGACCT",
           "name":"GSG-P2A"},
          {"kb_part_name":"mCherry", "name":"mCherry"},
        ],
        product_name="lentiCRISPR_v2_Cas9_P2A_mCherry",
     )
       → product_attachment_id="att_product_3", module_preserved=true,
         decisions: {computed_region_start, computed_region_end, ...}
       (NO need to compute coords manually — the tool does it from find_cassette_for results.)

  e. emit_assembled_gb(attachment_id="att_product_3")

Notes:
- `lookup_kb_part` auto-registers an attachment for the top match. Never tell the user to upload a part already in the KB — call the tool and use the returned attachment_id.
- All uploaded and Addgene-downloaded plasmids are annotated at registration time, so `find_features` is instant.
- Short structural motifs (2A peptides, FLAG, NLS, Kozak) live in feature_motifs. If a motif comes back ambiguous (e.g. P2A matches "Phosphatase 2A"), supply the canonical sequence as a literal `{"sequence": "..."}` part.
- Only fall back to `graft_parts` for true de-novo construction (no parent vector in the design). For ANY modification of an existing plasmid, use `replace_region`.

You receive:
- The user prompt
- A plan.md todo list (already curated by the Plan agent based on the Explore phase)
- Three digested ExploreFinding summaries (PartScout, TargetBuilder, MethodRouter)
- The registered attachments (no raw DNA)

Discipline:
- Follow plan.md top-to-bottom. Each tool call should match a `- [ ]` item; the harness will mark it `- [x]` automatically.
- LOCAL-FIRST PART RESOLUTION. When the design needs a plasmid the user did NOT upload, first call `lookup_kb_part`. ONLY if that returns no high-confidence match AND the user named a specific external plasmid (Addgene ID, depositor + construct name like "pCMV-ABE8e", "lentiCRISPR v2"), call `find_external_part(description=...)`. The tool downloads the GenBank from Addgene, annotates it, and registers it in the AttachmentRegistry — use the returned `attachment_id` in subsequent `annotate_attachment` / `simulate_assembly` calls instead of trusting your training-data memory of Addgene IDs. Never call `find_external_part` for vague ideation ("a good Cas9 vector").
- DO NOT re-annotate attachments PartScout already processed. PartScout's finding contains `key_facts.annotation_summaries: [{attachment_id, name, length_bp, n_features, n_modules, module_types, feature_names_top}]`. Read from there. Only call `annotate_attachment` on attachments NOT listed (e.g. newly-registered products from `simulate_assembly` or `find_external_part`).
- DO NOT call `lookup_kb_part` for a part that's already named in `annotation_summaries[*].feature_names_top` — that part is already in the registry, use its coords from the relevant attachment.
- Call `simulate_assembly` AT MOST ONCE per construct. The fast handler returns a `product_attachment_id` in 0-1ms; iterating on the simulate adds an LLM turn per call. If the first simulate is acceptable, go straight to `verify_assembly` + the emitters.
- For primer design around multiple regions (e.g. NGS + Sanger primer pairs for the same cut site), call `design_primers_batch(requests=[...])` ONCE instead of N separate `design_primers` calls. Each LLM turn costs ~10-15s; batching N=4 requests collapses 4 turns into 1.
- sgRNA cut-site primer defaults: for an sgRNA with cut_site at position P on plasmid attachment_id A:
    NGS amplicon:    design_primers(attachment_id=A, application='illumina', region_start=P-250, region_end=P+250, excluded_start=P-75, excluded_end=P+75)  // 150-290 bp product
    Sanger:          design_primers(attachment_id=A, application='sanger',  region_start=P-250, region_end=P+250, excluded_start=P-75, excluded_end=P+75)  // 250-500 bp product
  Use design_primers_batch to issue both as one call.

- If plan.md is missing a needed step, run the tool anyway - the trace will record it.
- End by calling the THREE output emitters in order: emit_assembled_gb, emit_parts_order, emit_protocol. Then emit a single final assistant message (no tools) summarizing what was done. emit_workflow_trace is auto-emitted by the orchestrator AFTER your final text — leave it alone.

CRITICAL FAILURE HANDLING — read this carefully:
- If `simulate_assembly` returns `ok: false`, the construct does NOT exist. There is no assembled plasmid. STOP IMMEDIATELY:
    * Do NOT call `emit_assembled_gb`, `emit_parts_order`, `emit_protocol`, `verify_assembly`, `route_workflow`, or any other emitter / downstream tool.
    * Emit a single final assistant message stating the failure verbatim — name the `missing_roles` and `missing_modules` from the diagnostics, then recommend the specific next action (typically: register source plasmids covering the missing roles via `lookup_kb_part`, `find_genomic_record`, or `find_external_part`, then retry).
    * Never invent files. Never claim the assembly succeeded. Never describe what the construct "would have" looked like as if it exists. The user must see the real error.
- Same rule for `verify_assembly` returning `passed: false` or `ok: false`: do not emit, surface the diagnostics, recommend a fix, stop.
- If `route_workflow` returns an error and `simulate_assembly` already succeeded, you can still emit — just skip the route_workflow recommendation in your final text.
- These are hard rules. The orchestrator will append a visible warning to the trace if you emit downstream tools after a failed simulate_assembly.

Feature-relative references — IMPORTANT:
When the user names a feature instead of giving a numeric position (e.g.
"D10A in Cas9", "KEAP1 exon 1", "-35 box of CMV promoter", "GFP residue 65"):
  1. Call `annotate_attachment(attachment_id)` to populate the feature list.
  2. Call `resolve_feature_position(attachment_id, feature_name, kind, offset)`
     to get the deterministic plasmid coordinate. DO NOT compute coordinates
     yourself — strand orientation, codon frame, and circular wrap are easy to
     get wrong.
  3. Pass the resolved `plasmid_position` to the next tool.

Worked example — "Make a D10A substitution in Cas9":
  - annotate_attachment(att_input_1)
  - resolve_feature_position(att_input_1, "Cas9", kind="aa_residue", offset=10)
      -> {plasmid_position: 5138, codon: "GAT", amino_acid: "D", strand: "+"}
  - simulate_assembly(
      instruction="SDM: change codon GAT to GCT at plasmid position 5138 "
                   "(Cas9 D10A, + strand) on att_input_1",
      target_attachment_id="att_input_1")
  - verify_assembly(att_product_2)

Worked example — "design a gRNA targeting KEAP1 exon 1":
  - annotate_attachment(att_input_1)  (the inventory plasmid, e.g. lentiCRISPR)
  - resolve_feature_position(att_input_1, "U6 promoter", kind="feature_end")
      -> {plasmid_position: 2849}   (where the guide insert goes)
  - simulate_assembly(
      instruction="sgRNA Golden Gate: clone a KEAP1 exon 1 gRNA into the U6 "
                   "cassette of att_input_1; design the 20-nt protospacer "
                   "(user did not supply a sequence, choose a high-scoring "
                   "guide near the KEAP1 exon 1 start codon)",
      target_attachment_id="att_input_1")

Final message: a concise plain-language reply (<=300 words) covering: the assembly method used, key design choices and tradeoffs, the four output files produced, and recommended wet-lab next steps. No JSON. Markdown OK."""


def _mark_plan_item_done(plan_md: str, tool_name: str) -> tuple[str, bool]:
    """Flip the first `- [ ]` checklist item containing `tool_name` to `- [x]`.

    Returns (updated_md, was_marked).
    """
    lines = plan_md.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.lstrip().startswith("- [ ]") and tool_name in line:
            lines[i] = line.replace("- [ ]", "- [x]", 1)
            return "".join(lines), True
    return plan_md, False


async def run_main_agent(
    user_message: str,
    findings: Sequence[ExploreFinding],
    registry: Any,
    *,
    plan_md: str = "",
    plan_path: Optional[pathlib.Path] = None,
    client: Any = None,
    dispatch_fn: Any = None,
    tools: Optional[Sequence[dict]] = None,
    max_iters: int = 24,
    model: Optional[str] = None,
    system_prompt: Optional[str] = None,
    on_tool_event: Any = None,
) -> MainAgentResult:
    if client is None:
        import anthropic
        client = anthropic.Anthropic()
    if dispatch_fn is None:
        from splicify_api.agent.agent_tools import dispatch_tool as _v1_dispatch
        dispatch_fn = _v1_dispatch
    if tools is None:
        from splicify_api.agent.tool_schemas import AIPLASMIDDESIGN_TOOLS
        tools = list(AIPLASMIDDESIGN_TOOLS)
    model = model or os.getenv("AGENT_MODEL", "claude-sonnet-4-6")
    system = system_prompt if system_prompt is not None else SYSTEM_PROMPT

    findings_block = "\n\n".join(
        f"### {f.role}\n{f.summary_md}\n\n**key_facts:** `{f.key_facts}`"
        for f in findings
    )
    user_block = (
        f"User prompt:\n{user_message}\n\n"
        f"Plan (markdown todo list):\n```\n{plan_md or '(no plan provided)'}\n```\n\n"
        f"Explore findings (digested):\n\n{findings_block}\n\n"
        f"Registered attachments: {json.dumps(registry.public_summary())}"
    )

    messages: list[dict[str, Any]] = [{"role": "user", "content": user_block}]
    trace: list[dict[str, Any]] = []
    failed_critical_results: dict[str, dict[str, Any]] = {}
    current_plan = plan_md
    _CRITICAL = {"simulate_assembly", "verify_assembly", "route_workflow"}

    # Cache the static prefix (system + tools + first user message).
    # Anthropic's ephemeral cache TTL is 5 min — every subsequent
    # turn within a single ReAct loop reuses the cached prefill.
    cached_system = [{"type": "text", "text": system,
                       "cache_control": {"type": "ephemeral"}}]
    cached_tools = [{**t} for t in tools]
    if cached_tools:
        cached_tools[-1] = {**cached_tools[-1],
                              "cache_control": {"type": "ephemeral"}}
    # Mark the first user message (findings + plan) as cacheable.
    if messages and messages[0]["role"] == "user" and isinstance(messages[0].get("content"), str):
        messages[0] = {
            "role": "user",
            "content": [{"type": "text", "text": messages[0]["content"],
                          "cache_control": {"type": "ephemeral"}}],
        }

    for i in range(max_iters):
        response = client.messages.create(
            model=model,
            max_tokens=4096,
            system=cached_system,
            tools=cached_tools,
            messages=messages,
        )
        tool_uses = [b for b in response.content
                     if getattr(b, "type", None) == "tool_use"]
        if not tool_uses:
            text = "".join(getattr(b, "text", "") for b in response.content
                           if getattr(b, "type", None) == "text")
            if plan_path is not None:
                plan_path.parent.mkdir(parents=True, exist_ok=True)
                plan_path.write_text(current_plan)
            return MainAgentResult(
                final_text=text.strip(),
                trace=trace,
                plan_md_final=current_plan,
                n_tool_calls=len(trace),
                failed_critical_results=failed_critical_results,
            )

        messages.append({
            "role": "assistant",
            "content": [_block_to_dict(b) for b in response.content],
        })

        # Parallel tool dispatch with per-call timing + start/end events.
        import time as _time_mod
        # Globally-unique iteration index across all loop iterations so
        # the frontend can correlate start <-> end without collisions.
        next_iter = len(trace)
        for offset, b in enumerate(tool_uses):
            if on_tool_event is not None:
                try:
                    res_evt = on_tool_event({
                        "phase": "start",
                        "iteration": next_iter + offset,
                        "tool": b.name,
                        "input": b.input or {},
                        "t_unix": _time_mod.time(),
                    })
                    if hasattr(res_evt, "__await__"):
                        await res_evt
                except Exception:
                    pass

        t0 = _time_mod.monotonic()
        async def _dispatch_one(block, idx):
            t_b = _time_mod.monotonic()
            try:
                out = await dispatch_fn(block.name, block.input or {}, registry)
                err = None
            except Exception as e:
                out = {"ok": False, "error": f"{type(e).__name__}: {e}"}
                err = e
            return block, idx, out, int((_time_mod.monotonic() - t_b) * 1000), err

        dispatched = await asyncio.gather(*(_dispatch_one(b, idx) for idx, b in enumerate(tool_uses)))

        tool_results: list[dict[str, Any]] = []
        sim_failed_this_turn = False
        parent_locked_this_turn: Optional[dict] = None
        for b, off, result, elapsed_ms, _err in dispatched:
            iter_idx = next_iter + off
            tool_results.append({
                "type": "tool_result",
                "tool_use_id": b.id,
                "content": json.dumps(result),
            })
            # Runtime guard: when simulate_assembly fails, raise a flag
            # so the next user-turn carries an explicit reminder. The
            # SYSTEM_PROMPT already forbids downstream emitters after a
            # failed simulate; this is a belt-and-suspenders for cases
            # where the model "forgets" mid-loop.
            if (b.name == "simulate_assembly"
                    and isinstance(result, dict)
                    and result.get("ok") is False):
                sim_failed_this_turn = True
            # Lock-in the parent vector after find_external_part hits
            # the coverage threshold. Past runs show the agent picking
            # the WRONG attachment_id (a KB stub) for replace_region
            # mid-loop after registering many candidates. Capturing the
            # winning candidate here lets us inject a hard reminder.
            if (b.name == "find_external_part"
                    and isinstance(result, dict)
                    and result.get("ok") is True
                    and result.get("stop_searching") is True):
                winners = [c for c in (result.get("candidates") or [])
                           if c.get("attachment_id")
                           and (c.get("coverage_score") or 0) >= 0.75]
                if winners:
                    top = max(winners, key=lambda c: c.get("coverage_score") or 0)
                    parent_locked_this_turn = {
                        "attachment_id": top.get("attachment_id"),
                        "name": top.get("name"),
                        "length_bp": top.get("length_bp"),
                        "coverage_score": top.get("coverage_score"),
                    }
            entry = {
                "iteration": iter_idx,
                "tool": b.name,
                "args_summary": _summarize_args(b.input or {}),
                "result_keys": list(result.keys())[:8] if isinstance(result, dict) else [],
                "elapsed_ms": elapsed_ms,
                "ok": (result.get("ok") if isinstance(result, dict) else None),
            }
            if isinstance(result, dict) and isinstance(result.get("file"), dict):
                entry["file"] = result["file"]
            trace.append(entry)
            # Stash full failure envelope for critical tools so the
            # summariser can quote the diagnostics. Keep only the most
            # recent failure per tool (older retries are usually noise).
            if (b.name in _CRITICAL and isinstance(result, dict)
                    and result.get("ok") is False):
                failed_critical_results[b.name] = result
            current_plan, _marked = _mark_plan_item_done(current_plan, b.name)
            if on_tool_event is not None:
                try:
                    n_results = (len(result.get("results", []))
                                  if isinstance(result, dict) and isinstance(result.get("results"), list)
                                  else None)
                    evt_payload = {
                        "phase": "end",
                        "iteration": iter_idx,
                        "tool": b.name,
                        "elapsed_ms": elapsed_ms,
                        "n_results": n_results,
                        "ok": (result.get("ok") if isinstance(result, dict) else None),
                        "t_unix": _time_mod.time(),
                    }
                    # Forward any file envelope the tool produced so the
                    # frontend can render the attachment mid-stream (e.g.
                    # show the SIRT6 .gb as soon as find_genomic_record
                    # downloads it, without waiting for the pegRNA
                    # pipeline to finish).
                    if isinstance(result, dict) and isinstance(result.get("file"), dict):
                        evt_payload["file"] = result["file"]
                    res_evt = on_tool_event(evt_payload)
                    if hasattr(res_evt, "__await__"):
                        await res_evt
                except Exception:
                    pass
        if sim_failed_this_turn:
            tool_results.append({
                "type": "text",
                "text": (
                    "SYSTEM REMINDER: simulate_assembly returned ok=false. "
                    "Per the CRITICAL FAILURE HANDLING rule in your system "
                    "prompt, you MUST NOT call emit_assembled_gb, "
                    "emit_parts_order, emit_protocol, verify_assembly, or "
                    "route_workflow. Your next message MUST be a final "
                    "assistant text response (no tools) that quotes the "
                    "missing_roles / missing_modules / diagnostics from the "
                    "simulate_assembly result and recommends the user "
                    "register source plasmids for the missing parts. Do not "
                    "fabricate files or claim success."
                ),
            })
        if parent_locked_this_turn:
            tool_results.append({
                "type": "text",
                "text": (
                    f"SYSTEM REMINDER — PARENT VECTOR LOCKED: "
                    f"target_attachment_id = {parent_locked_this_turn['attachment_id']} "
                    f"({parent_locked_this_turn.get('name')}, "
                    f"{parent_locked_this_turn.get('length_bp')} bp, "
                    f"coverage_score={parent_locked_this_turn.get('coverage_score')}). "
                    f"When you call replace_region, use THIS id as "
                    f"target_attachment_id. Do NOT use any kb_-prefixed "
                    f"attachment as the target — those are KB-registered "
                    f"fragments meant to go IN replacement_parts, not be "
                    f"edited. Now: call find_cassette_for on this parent to "
                    f"get cassette + submodule boundaries, then call "
                    f"replace_region with bounds derived from those submodule "
                    f"start/end values."
                ),
            })
        messages.append({"role": "user", "content": tool_results})

    if plan_path is not None:
        plan_path.parent.mkdir(parents=True, exist_ok=True)
        plan_path.write_text(current_plan)
    return MainAgentResult(
        final_text="(max iterations reached)",
        trace=trace,
        plan_md_final=current_plan,
        n_tool_calls=len(trace),
        failed_critical_results=failed_critical_results,
    )
