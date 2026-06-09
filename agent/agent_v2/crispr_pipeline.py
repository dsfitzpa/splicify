"""CRISPR-design pipeline — replaces the stub_crispr.respond branch.

Mirrors the plasmid-cloning pipeline shape (Explore -> Plan -> Main ->
Summarizer) with a slimmer Explore phase: just TargetLocator (where to
edit) and GuideStrategist (how to edit). Plan and Main agents reuse the
existing implementations with CRISPR-specific system prompts swapped in.

Inputs: user_message + AttachmentRegistry (one or more .gb files registered).
Outputs: the same /agent_v2/chat envelope shape the plasmid pipeline returns.
"""
from __future__ import annotations

import asyncio
import pathlib
import time
from typing import Any, Optional

import agent_v2  # noqa: F401 - triggers path shim
from agent_v2 import memory, tools as v2_tools, crispr_clarification
from agent_v2.explore import ExploreFinding
from agent_v2.main_agent import run_main_agent
from agent_v2.subagents.guide_strategist import run_guide_strategist
from agent_v2.subagents.plan_agent import plan_path_for, run_plan_agent
from agent_v2.subagents.summarizer import run_summarizer
from agent_v2.subagents.target_locator import run_target_locator


CRISPR_PLAN_SYSTEM_PROMPT = """You are the Plan agent for the CRISPR-design workflow of an AI molecular-biology agent. You receive:
1. The user's original prompt.
2. Two digested ExploreFinding summaries (TargetLocator, GuideStrategist).

You DO NOT have raw DNA. You DO NOT have tools.

Your job: emit a SHORT markdown todo list that the Main Agent will execute. Each item must be a concrete tool call. The Main Agent has these tools:
- resolve_feature_position (if more residue->coord mapping is needed)
- design_guides (Cas9 / Cas12a sgRNA design — args: attachment_id, region_start, region_end, pam, guide_length, max_guides, score_method)
- design_pegrnas (prime editing — args: attachment_id, edit_start, edit_end, alt, edit_type, n_results, use_pe3)
- design_primers (NGS / Sanger primers — args: attachment_id, region_start, region_end, excluded_start, excluded_end, application='illumina'|'sanger', product_size_*)
- annotate_attachment, lookup_kb_part, analyze_design_intent, verify_assembly
- Output emitters: emit_guides_csv (sgRNA/pegRNA/primer/oligo table), emit_guides_gb (target .gb + appended sgRNA + primer features), emit_parts_order (sgRNA oligos + primer oligos), emit_protocol (emit_assembled_gb only when the user asked for an sgRNA-cloned plasmid). emit_workflow_trace is auto-emitted by the orchestrator AFTER your final text — DO NOT call it.

The plan should:
- Lead with the strategy from GuideStrategist (Cas9 sgRNA vs pegRNA, n_targets, PAM, scoring method).
- **Batch when n_targets > 1**: write ONE plan item that batches design_guides across all targets in parallel, then ONE plan item that batches design_primers (illumina) across all targets in parallel, then ONE plan item for design_primers (sanger). The Main Agent's harness will fan these out via asyncio.gather. Do NOT write a separate item per target.
- For each target the TargetLocator resolved: design the guide (design_guides for Cas9, design_pegrnas for prime edit). For multi-target prompts, ALL design_guides calls land in ONE plan item.
- After all guide designs: design NGS primers (application='illumina', region = guide +/-250 bp, excluded = guide +/-75 bp, product_size 150..290), then Sanger primers (application='sanger', region = guide +/-250 bp, excluded = guide +/-75 bp, product_size 250..500). One plan item per primer application, regardless of n_targets.
- End with: emit_guides_csv (every guide + primer + cloning oligo), emit_guides_gb (target_attachment_id with guides + primers appended as misc_RNA / primer_bind features), emit_parts_order (sgRNA oligos + primer oligos), emit_protocol (assembly_method='sgrna_gg' with NGS/Sanger sequencing steps appended via custom_steps). emit_workflow_trace is auto-emitted by the orchestrator — don't call it.

Format: one `## Plan` heading, then a Markdown checklist (`- [ ] 1. ...`). <=15 items. No commentary outside the heading + list."""


CRISPR_MAIN_SYSTEM_PROMPT = """You are the Main Agent for the CRISPR-design workflow of an AI molecular-biology agent. You orchestrate guide design, primer design for the editing readout, and output-file emission.

You have a full tool roster:
- annotate_attachment, analyze_design_intent, verify_assembly, route_workflow (lookup_kb_part is available but the CRISPR pipeline doesn't need it for standard nucleases / pegRNA scaffolds — skip it)
- resolve_feature_position
- design_guides (Cas9/Cas12a sgRNA), design_pegrnas (prime editing)
- design_primers (NGS / Sanger primers around the cut site)
- Output emitters: emit_guides_csv, emit_guides_gb, emit_assembled_gb (rarely used here), emit_parts_order, emit_protocol
  (emit_workflow_trace is auto-emitted by the orchestrator AFTER your
  final text response — do not call it; do not compose its args.)

You receive:
- The user prompt
- A plan.md todo list curated by the Plan agent (the harness will cross items off as you call matching tools)
- Two digested ExploreFinding summaries (TargetLocator, GuideStrategist)
- The registered attachments (no raw DNA)

Discipline:
- Follow plan.md top-to-bottom. Each tool call should match a `- [ ]` item.

- **PARALLEL TOOL CALLS ARE MANDATORY WHEN n_targets > 1. Read this carefully.**

  When the GuideStrategist's `n_targets` > 1 (e.g. 2 pegRNAs for R15C + CCR13-15RTC, or 3 sgRNAs across 3 residues), you MUST emit ALL the tool_use blocks for a given design step IN ONE RESPONSE. The harness dispatches every tool_use block in a single Claude response via `asyncio.gather` — sequential responses serialise the work and roughly multiply wall-clock by n_targets.

  Concrete batching for a 2-target prompt:

  - Response k:   2x resolve_feature_position (one per target)         <- IN ONE RESPONSE
  - Response k+1: 2x design_pegrnas (or design_guides)                  <- IN ONE RESPONSE
  - Response k+2: 2x design_primers(application='illumina')             <- IN ONE RESPONSE
  - Response k+3: 2x design_primers(application='sanger')               <- IN ONE RESPONSE
  - Response k+4: emit_guides_csv (ONE call carrying ALL targets)
  - Response k+5: emit_guides_gb  (ONE call carrying ALL targets)
  - Response k+6: emit_parts_order, emit_protocol  <- IN ONE RESPONSE
  - emit_workflow_trace is auto-emitted by the orchestrator after the loop

  Doing it sequentially (one tool_use per response) is wrong and ~2x slower for 2 targets. If you find yourself about to send a response with only a single design_* tool_use block for a multi-target prompt, STOP and add the sibling tool_use block(s) for the other targets into the same response.

  Plan.md crossing-off only flips one item per response, but the wall-clock saving comes from parallel dispatch, not from item-crossing pacing.

  The emitters (emit_guides_csv / emit_guides_gb / emit_parts_order / emit_protocol) are ALWAYS called ONCE each with aggregated args carrying every target's guides + primers + cloning oligos — never one emitter call per target. emit_workflow_trace is auto-emitted server-side after the loop, not by you.
- One design_pegrnas (or design_guides) call PER ENTRY in TargetLocator's `resolved_targets`. If the locator returned N targets, you MUST emit N design calls — never drop a target, never merge separate single-codon targets into a multi-codon span. Sanity-check yourself before emit_guides_csv: the descriptor + pegRNA names must mention every residue from the user's prompt.
- **STRICT RETRY BUDGET — 3 ATTEMPTS PER TARGET**. If a single target's design_pegrnas call fails 3 times in a row (length mismatch, no valid spacer in window, span crosses intron, etc.), ABANDON that target. Do NOT make a 4th attempt. Record the abandonment in your final summary with the failure reason. Common reasons to surface verbatim to the user: "spans an intron — pegRNA editing cannot read across exon boundaries"; "no NGG PAM within the spacer search window"; "residue numbers fall outside the picked isoform — confirm the protein accession". Proceed with whatever targets did design successfully.
- **ALL-TARGETS-FAILED EXIT**. If every target hit the 3-failure budget and there are zero successful designs, SKIP emit_guides_csv / emit_guides_gb / emit_parts_order / emit_protocol entirely. The pipeline should NOT emit empty output files. End your turn with a final assistant message stating which targets failed and why; the orchestrator will surface that as the chat reply.
- For each target:
  1. resolve_feature_position (skip if TargetLocator already resolved it).
     MULTI-RESIDUE / SPANNING edits (e.g. "E260D / C263T / R264Q" or any
     RTC-style combined edit): call resolve_feature_position for BOTH the
     LOWEST and HIGHEST residue, then combine into one design_pegrnas call
     where edit_start_1based = the lower bound and edit_end_1based = the
     upper bound from the resolver. Build `alt` as the + strand sequence:
     concatenate (in + strand order) the + strand bases for each codon,
     using the resolver's plus_strand_ref as the per-codon template and
     replacing the codon for each substituted residue with revcomp(new
     sense codon) on - strand genes (or just the new codon on + strand
     genes). NEVER guess edit_start/edit_end/alt for multi-residue spans
     — that is the #1 cause of wasted failed design_pegrnas calls.
  2. design_guides (Cas9) or design_pegrnas (prime editing). Region for design_guides = resolved plasmid_position +/-30 bp; max_guides=5.
  3. Pick the highest-scoring guide.
  4. design_primers application='illumina' for NGS amplicons. region = guide_start - 250 to guide_end + 250 (clamp to plasmid bounds). excluded_start = guide_start - 75, excluded_end = guide_end + 75.
  5. design_primers application='sanger' with the same region + excluded.
- After all targets, call the four emitters in order:
  1. emit_guides_csv — every sgRNA / pegRNA / ngRNA / primer / cloning oligo as one wide row each. Cloning oligos (oligo_top / oligo_bottom) are pre-computed by you for the user's destination vector (pX330=BbsI: top='CACCG'+spacer, bottom='AAAC'+revcomp(spacer)+'C'; lentiCRISPR=BsmBI).
  2. emit_guides_gb — target_attachment_id + descriptor ONLY. The dispatcher reuses the pegrnas/guides/primers/cloning_oligos you just passed to emit_guides_csv from a process-local cache, so RE-EMITTING them wastes tokens and slows the run by 30-60 s per repetition. Just pass target_attachment_id + descriptor.
  3. emit_parts_order — pass `descriptor` ONLY. The emitter auto-derives the parts list (every spacer cloning oligo top/bottom + every NGS / Sanger primer forward/reverse) from the data you just passed to emit_guides_csv. Re-emitting `parts=[...]` from scratch wastes 30-60 s of token generation per turn.
  4. emit_protocol — pass `assembly_method='sgrna_gg'` + `descriptor` ONLY. The template covers anneal -> ligate -> transform -> miniprep -> Sanger; ONLY pass `custom_steps` if the user explicitly asked for a non-default protocol step. Composing custom_steps verbatim on every run stalls the pipeline.
  After step 4: emit your final summary text. The orchestrator emits workflow_trace.txt automatically.
- Final assistant message: <=300 words summarising the chosen guides, primers, and the next wet-lab steps. No raw DNA in the reply. List each guide's spacer / score / PAM, and each primer pair's name + Tm + amplicon size."""


async def _run_crispr_pipeline(
    *,
    user_message: str,
    registry: Any,
    session_id: str,
    state: Optional[memory.SessionState],
    deps: Any,
    on_tool_event: Any = None,
) -> dict:
    # Per-turn output dir + plan path (mirrors _run_plasmid_pipeline).
    turn_id = f"turn_{(len(state.messages) // 2 + 1) if state else 1}_{int(time.time())}"
    turn_dir: Optional[str] = None
    plan_path: Optional[pathlib.Path] = None
    if deps.output_dir is not None:
        turn_dir_path = pathlib.Path(deps.output_dir) / session_id / turn_id
        turn_dir_path.mkdir(parents=True, exist_ok=True)
        turn_dir = str(turn_dir_path)
        plan_path = plan_path_for(session_id, turn_id, base=deps.output_dir)

    # Dispatch_fn binds output_dir for emit_* calls.
    if deps.dispatch_fn is not None:
        dispatch_fn = deps.dispatch_fn
    else:
        if turn_dir is None:
            dispatch_fn = v2_tools.dispatch_with_emitters
        else:
            async def _bound(name, args, registry, _td=turn_dir):
                return await v2_tools.dispatch_with_emitters(name, args, registry, output_dir=_td)
            dispatch_fn = _bound

    # Default tools: full roster (v1 + emitters + resolver + crispr tools).
    tools = deps.tools if deps.tools is not None else v2_tools.make_full_tool_roster()

    # Slim Explore phase: 2 subagents in parallel.
    findings: tuple[ExploreFinding, ExploreFinding] = await asyncio.gather(
        run_target_locator(
            user_message, registry,
            client=getattr(deps, "target_locator_client", None),
            dispatch_fn=dispatch_fn,
            on_tool_event=on_tool_event,
        ),
        run_guide_strategist(
            user_message, registry,
            client=getattr(deps, "guide_strategist_client", None),
            dispatch_fn=dispatch_fn,
            on_tool_event=on_tool_event,
        ),
    )

    # Missing-information short-circuit. Each Explore subagent may emit a
    # `missing_info` list of strings in its key_facts. If any are present,
    # skip Plan + Main + Summarizer entirely and return a clarification
    # envelope asking the user for the missing piece(s). Saves ~20-40 k
    # Anthropic tokens per under-specified request.
    missing: list[str] = []
    for f in findings:
        kf_missing = (f.key_facts or {}).get("missing_info") or []
        for entry in kf_missing:
            if entry and entry not in missing:
                missing.append(entry)
    if missing:
        envelope = crispr_clarification.respond(missing, user_message=user_message)
        envelope["turn_id"] = turn_id
        envelope["findings"] = [
            {"role": f.role, "summary_md": f.summary_md, "key_facts": f.key_facts}
            for f in findings
        ]
        return envelope

    plan_result = await run_plan_agent(
        user_message, findings,
        plan_path=plan_path,
        client=deps.plan_client,
        system_prompt=CRISPR_PLAN_SYSTEM_PROMPT,
    )

    main_result = await run_main_agent(
        user_message, findings, registry,
        plan_md=plan_result.plan_md, plan_path=plan_path,
        client=deps.main_client, dispatch_fn=dispatch_fn, tools=tools,
        system_prompt=CRISPR_MAIN_SYSTEM_PROMPT,
        on_tool_event=on_tool_event,
    )

    files = [entry["file"] for entry in main_result.trace if entry.get("file")]

    from agent_v2.outputs.workflow_trace import auto_emit_workflow_trace
    try:
        wt_result = await auto_emit_workflow_trace(
            session_id=session_id, turn_id=turn_id,
            main_result=main_result, findings=findings,
            plan_md=main_result.plan_md_final or "",
            output_dir=turn_dir, registry=registry,
        )
        if wt_result.get("ok") and wt_result.get("file"):
            files.append(wt_result["file"])
    except Exception:
        pass

    files = files or None

    summary = await run_summarizer(
        user_message,
        main_reply=main_result.final_text,
        files=files,
        findings=findings,
        decisions=(state.decisions if state else None),
        client=deps.summarizer_client,
    )

    return {
        "ok": True,
        "reply": summary.reply,
        "main_agent_draft": main_result.final_text,
        "files": files,
        "viz": None,
        "agent_trace": main_result.trace,
        "n_tool_calls": main_result.n_tool_calls,
        "error": None,
        "plan_md": main_result.plan_md_final,
        "findings": [
            {"role": f.role, "summary_md": f.summary_md, "key_facts": f.key_facts}
            for f in findings
        ],
        "turn_id": turn_id,
        "workflow": "crispr",
    }
