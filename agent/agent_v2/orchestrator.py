"""Orchestrator — the entry point that ties the agent_v2 pipeline together.

Flow:
  triage -> branch -> {rejection | crispr_stub | plasmid_pipeline} -> memory_save -> envelope

The plasmid pipeline runs the three Explore subagents in parallel via
asyncio.gather, then the Plan agent (digested findings only), then the Main
agent (full tool roster + plan.md crossing-off). When `deps.output_dir` is
set, plan.md and all emitter outputs land at `<output_dir>/<sid>/<tid>/`,
and `envelope["files"]` is populated from the main agent's trace.
"""
from __future__ import annotations

import asyncio
import pathlib
import time
from dataclasses import dataclass
from typing import Any, Awaitable, Callable, Optional

from agent_v2 import memory, rejection, stub_crispr, tools as v2_tools
from agent_v2.crispr_pipeline import _run_crispr_pipeline
from agent_v2.explore import ExploreFinding
from agent_v2.main_agent import run_main_agent
from agent_v2.subagents.method_router import run_method_router
from agent_v2.subagents.part_scout import run_part_scout
from agent_v2.subagents.plan_agent import plan_path_for, run_plan_agent
from agent_v2.subagents.summarizer import run_summarizer
from agent_v2.subagents.target_builder import run_target_builder
from agent_v2.triage import triage as _triage


@dataclass
class OrchestratorDeps:
    triage_client: Any = None
    part_scout_client: Any = None
    target_builder_client: Any = None
    method_router_client: Any = None
    target_locator_client: Any = None
    guide_strategist_client: Any = None
    plan_client: Any = None
    main_client: Any = None
    summarizer_client: Any = None
    dispatch_fn: Any = None
    tools: Any = None
    skip_memory: bool = False
    output_dir: Optional[str] = None  # base; emitters + plan.md land at <output_dir>/<sid>/<tid>/


def _make_dispatch_with_output_dir(output_dir: Optional[str]):
    """Wrap dispatch_with_emitters so it binds output_dir for every call."""
    if output_dir is None:
        return v2_tools.dispatch_with_emitters

    async def _bound(name, args, registry):
        return await v2_tools.dispatch_with_emitters(name, args, registry, output_dir=output_dir)
    return _bound


async def run_orchestrator(
    user_message: str,
    registry: Any,
    *,
    session_id: Optional[str] = None,
    deps: Optional[OrchestratorDeps] = None,
    on_triage: Optional[Callable[[Any], Awaitable[None]]] = None,
    on_tool_event: Optional[Callable[[dict], Any]] = None,
) -> dict:
    """Run the agent_v2 pipeline. Returns the /agent_v2/chat-shaped envelope.

    `on_triage` (optional) is awaited as soon as the triage classifier finishes
    so callers (e.g. the SSE chat-stream endpoint) can emit the shorthand
    summary to the client before the heavy work begins.
    """
    deps = deps or OrchestratorDeps()

    if not session_id:
        session_id = memory.mint_session_id()

    state = None if deps.skip_memory else memory.load(session_id)
    has_attachments = bool(registry.public_summary())

    triage_result = _triage(
        user_message,
        has_attachments=has_attachments,
        client=deps.triage_client,
    )
    if on_triage is not None:
        await on_triage(triage_result)

    if triage_result.intent == "REJECT":
        envelope = rejection.respond(
            shorthand=triage_result.shorthand,
            reason=triage_result.rejection_reason,
        )
    elif triage_result.intent == "CRISPR_GUIDE":
        envelope = await _run_crispr_pipeline(
            user_message=user_message,
            registry=registry,
            session_id=session_id,
            state=state,
            deps=deps,
            on_tool_event=on_tool_event,
        )
    elif triage_result.intent == "PLASMID_QA":
        envelope = await _run_qa_pipeline(
            user_message=user_message,
            registry=registry,
            deps=deps,
            on_tool_event=on_tool_event,
        )
        # QA-pipeline redirect: re-route into the requested workflow.
        from agent_v2.redirect import is_redirect_envelope
        if is_redirect_envelope(envelope):
            redirected_from = triage_result.intent
            redirect_reason = envelope.get("redirect_reason", "")
            new_intent = envelope.get("redirect_to")
            if on_triage is not None:
                try:
                    await on_triage(type("T", (), {
                        "intent": new_intent,
                        "shorthand": f"Redirected: {redirect_reason}"[:80],
                        "is_new_topic": triage_result.is_new_topic,
                        "rejection_reason": None,
                    })())
                except Exception:
                    pass
            if new_intent == "PLASMID_CLONING":
                envelope = await _run_plasmid_pipeline(
                    user_message=user_message, registry=registry,
                    session_id=session_id, state=state, deps=deps,
                )
            elif new_intent == "CRISPR_GUIDE":
                envelope = await _run_crispr_pipeline(
                    user_message=user_message, registry=registry,
                    session_id=session_id, state=state, deps=deps,
                )
            elif new_intent == "REJECT":
                envelope = rejection.respond(
                    shorthand=triage_result.shorthand,
                    reason=redirect_reason or triage_result.rejection_reason,
                )
            envelope["redirected_from"] = redirected_from
            envelope["redirect_reason"] = redirect_reason
            envelope["intent"] = new_intent
    else:
        envelope = await _run_plasmid_pipeline(
            user_message=user_message,
            registry=registry,
            session_id=session_id,
            state=state,
            deps=deps,
            on_tool_event=on_tool_event,
        )

    envelope["session_id"] = session_id
    envelope["intent"] = envelope.get("intent") or triage_result.intent
    envelope["shorthand"] = triage_result.shorthand
    envelope["is_new_topic"] = triage_result.is_new_topic

    if envelope.get("redirected_from"):
        from agent_v2.redirect import merge_redirect_into_reply
        envelope = merge_redirect_into_reply(envelope)

    if not deps.skip_memory:
        prev_messages = state.messages if state else []
        prev_decisions = state.decisions if state else []
        new_state = memory.SessionState(
            session_id=session_id,
            messages=prev_messages + [
                {"role": "user", "content": user_message},
                {"role": "assistant", "content": envelope.get("reply", "")},
            ],
            registry_summary=registry.public_summary(),
            decisions=prev_decisions,
            last_user_message=user_message,
        )
        memory.save(new_state)

    return envelope


async def _run_qa_pipeline(
    *,
    user_message: str,
    registry: Any,
    deps: OrchestratorDeps,
    on_tool_event: Optional[Callable[[dict[str, Any]], Any]] = None,
) -> dict:
    """Plasmid Q&A — routes to the deterministic interpreter (PlasmidIndex
    + Sonnet tool-use loop) instead of the cloning subagent pipeline.

    Walks the AttachmentRegistry, annotates each plasmid through
    annotate_llm_cached (cache hits when the user already uploaded the
    file via /annotate-on-upload), builds a PlasmidRegistry, and asks
    the interpreter. Wraps the answer in the same chat-envelope shape
    that the cloning + crispr + rejection branches produce so the
    frontend doesn't need a separate render path."""
    from splicify_api.annotation_cache import annotate_llm_cached
    from agent_v2.interpreter.plasmid_registry import PlasmidRegistry
    from agent_v2.interpreter.agent import run_interpreter

    qa_reg = PlasmidRegistry()
    skipped: list[str] = []
    for att in registry.items.values():
        try:
            env = await annotate_llm_cached(att.sequence, circular=att.circular)
            env.setdefault("sequence", att.sequence)
            qa_reg.register(att.attachment_id, env, name=att.name)
        except Exception as e:
            skipped.append(f"{att.attachment_id}: {type(e).__name__}: {e}")

    # Empty registry is intentionally NOT short-circuited — the
    # interpreter's find_external_part tool can still answer "do you know
    # about lentiCRISPR v2?" without any local plasmids. The system
    # prompt handles the empty case directly.

    try:
        result = await run_interpreter(user_message, qa_reg,
                                       client=deps.main_client,
                                       on_tool_event=on_tool_event)
    except Exception as e:
        return {
            "ok": False,
            "reply": f"Interpreter failed: {type(e).__name__}: {e}",
            "files": None, "viz": None,
            "agent_trace": [],
            "n_tool_calls": 0,
            "error": str(e),
        }

    # If the interpreter asked to redirect, surface it as a redirect-
    # shaped envelope. The orchestrator picks this up and reruns the
    # request through the requested workflow.
    if getattr(result, "redirect", None) is not None:
        return {
            "ok": True,
            "redirect_to": result.redirect.target_intent,
            "redirect_reason": result.redirect.reason,
            "redirect_notes": getattr(result.redirect, "notes", ""),
            "redirect_findings": getattr(result.redirect, "findings", {}) or {},
            "files": None, "viz": None,
            "agent_trace": [
                {"iteration": i, "tool": t.get("tool"),
                 "args_summary": str(t.get("input", "")),
                 "result_keys": [],
                 "n_results": t.get("n_results")}
                for i, t in enumerate(result.trace)
            ],
            "n_tool_calls": result.n_tool_calls,
            "error": None,
        }

    return {
        "ok": True,
        "reply": result.answer,
        "files": None,
        "viz": None,
        "agent_trace": [
            {"iteration": i, "tool": t.get("tool"),
             "args_summary": str(t.get("input", "")),
             "result_keys": [],
             "n_results": t.get("n_results")}
            for i, t in enumerate(result.trace)
        ],
        "n_tool_calls": result.n_tool_calls,
        "error": None,
        "citations": result.citations,
        "skipped_attachments": skipped,
    }





async def _emit_phase(on_tool_event, phase_name, iteration, phase="start", elapsed_ms=0):
    if on_tool_event is None:
        return
    import time as _t
    try:
        evt = on_tool_event({"phase": phase, "iteration": iteration,
                              "tool": phase_name, "input": {},
                              "elapsed_ms": elapsed_ms, "n_results": None,
                              "ok": True, "t_unix": _t.time()})
        if hasattr(evt, "__await__"):
            await evt
    except Exception:
        pass


async def _run_phase(on_tool_event, name, iteration, coro):
    import time as _t
    await _emit_phase(on_tool_event, name, iteration, phase="start")
    t0 = _t.monotonic()
    try:
        result = await coro
    finally:
        elapsed_ms = int((_t.monotonic() - t0) * 1000)
    await _emit_phase(on_tool_event, name, iteration, phase="end",
                      elapsed_ms=elapsed_ms)
    return result


async def _run_plasmid_pipeline(
    *,
    user_message: str,
    registry: Any,
    session_id: str,
    state: Optional[memory.SessionState],
    deps: OrchestratorDeps,
    on_tool_event: Any = None,
) -> dict:
    # Resolve per-turn output dir + plan path
    turn_id = f"turn_{(len(state.messages) // 2 + 1) if state else 1}_{int(time.time())}"
    turn_dir: Optional[str] = None
    plan_path: Optional[pathlib.Path] = None
    if deps.output_dir is not None:
        turn_dir_path = pathlib.Path(deps.output_dir) / session_id / turn_id
        turn_dir_path.mkdir(parents=True, exist_ok=True)
        turn_dir = str(turn_dir_path)
        plan_path = plan_path_for(session_id, turn_id, base=deps.output_dir)

    # Build the dispatch_fn that binds output_dir for emit_* calls.
    if deps.dispatch_fn is not None:
        dispatch_fn = deps.dispatch_fn
    else:
        dispatch_fn = _make_dispatch_with_output_dir(turn_dir)

    # Default tools to the full roster (v1 + 4 emitters) when not injected.
    tools = deps.tools if deps.tools is not None else v2_tools.make_full_tool_roster()

    # Hard 120 s cap on the Explore phase. Parallel sessions share the
    # Anthropic API key + NCBI/PubMed rate gates; without the cap, a
    # runaway Explore in one session can stall every other request.
    # On timeout we surface placeholder findings and let the Main
    # agent proceed with local KB knowledge.
    EXPLORE_TIMEOUT_S = 120.0
    explore_coro = asyncio.gather(
        run_part_scout(user_message, registry,
                        client=deps.part_scout_client, dispatch_fn=dispatch_fn),
        run_target_builder(user_message, registry,
                            client=deps.target_builder_client, dispatch_fn=dispatch_fn),
        run_method_router(user_message, registry,
                           client=deps.method_router_client, dispatch_fn=dispatch_fn),
        return_exceptions=True,
    )
    await _emit_phase(on_tool_event, "_phase_explore", -100, phase="start")
    import time as _t_explore
    _t_explore_start = _t_explore.monotonic()
    try:
        raw_findings = await asyncio.wait_for(explore_coro, timeout=EXPLORE_TIMEOUT_S)
    except asyncio.TimeoutError:
        raw_findings = []
    elapsed_ms = int((_t_explore.monotonic() - _t_explore_start) * 1000)
    await _emit_phase(on_tool_event, "_phase_explore", -100, phase="end",
                      elapsed_ms=elapsed_ms)
    # Wrap any per-subagent exception into a placeholder so the rest
    # of the pipeline doesn't see None or raise.
    findings_list = []
    for role, item in zip(("part_scout", "target_builder", "method_router"),
                           raw_findings if raw_findings else []):
        if isinstance(item, ExploreFinding):
            findings_list.append(item)
        else:
            err = (type(item).__name__ if isinstance(item, BaseException)
                   else "timeout or empty result")
            findings_list.append(ExploreFinding(
                role=role,
                summary_md=f"(Explore subagent {role} aborted: {err})",
                key_facts={"aborted": True},
                references=[a["attachment_id"] for a in registry.public_summary()],
                trace=[],
            ))
    while len(findings_list) < 3:
        findings_list.append(ExploreFinding(
            role="missing", summary_md="(missing)", key_facts={},
            references=[], trace=[],
        ))
    findings = tuple(findings_list)

    plan_result = await _run_phase(
        on_tool_event, "_phase_plan", -99,
        run_plan_agent(
            user_message, findings,
            plan_path=plan_path, client=deps.plan_client,
        ),
    )

    await _emit_phase(on_tool_event, "_phase_main", -98, phase="start")
    main_result = await run_main_agent(
        user_message, findings, registry,
        plan_md=plan_result.plan_md, plan_path=plan_path,
        client=deps.main_client, dispatch_fn=dispatch_fn, tools=tools,
        on_tool_event=on_tool_event,
    )

    # Collect any file envelopes the main agent's emitters produced.
    files = [entry["file"] for entry in main_result.trace if entry.get("file")]

    # Safety-net: when the Main agent skipped emit_assembled_gb but the
    # registry has attachments, synthesize an assembled.gb server-side
    # from the registered products. Guarantees the user always gets the
    # headline deliverable + a visualization, even when the LLM gave up
    # threading attachment_ids through simulate_assembly.
    emitted_assembled = any(
        (e.get("tool") == "emit_assembled_gb") for e in (main_result.trace or [])
    )
    if not emitted_assembled and hasattr(registry, "items") and registry.items:
        try:
            from agent_v2.outputs.assembled_gb import emit_assembled_gb
            from agent_v2.tools import simulate_assembly_fast
            # Run the fast deterministic assembler with empty args so it
            # falls through to "concatenate everything registered."
            sim_result = await simulate_assembly_fast(
                {"target_attachment_id": None,
                 "inventory_attachment_ids": [],
                 "instruction": "auto-assembled from registry on Main-agent skip"},
                registry, output_dir=turn_dir,
            )
            new_aid = sim_result.get("product_attachment_id")
            if new_aid:
                a_result = await emit_assembled_gb(
                    {"attachment_id": new_aid}, registry, output_dir=turn_dir,
                )
                if isinstance(a_result, dict) and a_result.get("file"):
                    files.append(a_result["file"])
                    main_result.trace.append({
                        "iteration": 999,
                        "tool": "emit_assembled_gb",
                        "args_summary": f"auto-emitted from registry (product={new_aid})",
                        "result_keys": list(a_result.keys())[:8],
                        "file": a_result["file"],
                    })
        except Exception as _e:
            # Best-effort fallback; never break the pipeline on this.
            main_result.trace.append({
                "iteration": 999,
                "tool": "emit_assembled_gb",
                "args_summary": f"auto-emit failed: {type(_e).__name__}: {_e}",
                "result_keys": [],
            })

    # Server-side workflow_trace emission. The LLM stopped being allowed
    # to call emit_workflow_trace (it kept hanging the run while composing
    # the args dict + final text in one turn, pushing past Vercel's
    # 5-minute SSE-proxy timeout). Emit it here instead — same dispatcher
    # path so output_dir + base64 envelope contract is identical.
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
    except Exception as e:
        # Trace emission is non-critical — surface the failure on the
        # envelope but let the pipeline keep going.
        envelope_workflow_trace_error = f"{type(e).__name__}: {e}"
    else:
        envelope_workflow_trace_error = None

    files = files or None

    # Polish the final reply via the summarizer (falls back to main draft on error).
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
        "viz": None,    # populated by router after pipeline (re-annotates products)
        "agent_trace": main_result.trace,
        "n_tool_calls": main_result.n_tool_calls,
        "error": None,
        "plan_md": main_result.plan_md_final,
        "findings": [
            {"role": f.role, "summary_md": f.summary_md, "key_facts": f.key_facts}
            for f in findings
        ],
        "turn_id": turn_id,
    }
