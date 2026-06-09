"""
AI-powered explanation generation using OpenAI.
Generates expert, human-readable explanations of bioinformatics design results.
"""
from __future__ import annotations

import os
from typing import Any, Dict, Optional


_GIBSON_SYSTEM = """You are an expert molecular biologist with 15+ years of cloning experience.
You have just completed a Gibson Assembly primer design. Write a concise expert explanation
(175–275 words, 3 paragraphs) covering:
1. What was designed — fragment count, junction strategy, key overlap/annealing parameters. If pre-design routing info is provided, briefly explain why Gibson Assembly was selected (e.g., fragment count, cost, feasibility vs. alternatives).
2. Quality flags — flag any overlap scores <60, uniqueness <70, strong hairpins (dG < -3 kcal/mol), or warnings; confirm good quality if all pass
3. Next steps — PCR amplification, assembly reaction conditions, verification strategy

Use plain scientific language. Reference specific numbers from the data. Do not repeat generic protocol boilerplate."""

_PCR_SYSTEM = """You are a PCR primer design expert with 15+ years of experience.
You have just designed primers for a PCR reaction. Write a concise expert explanation
(150–200 words, 3 paragraphs) covering:
1. Primer quality — Tm values, length, any secondary structure concerns (hairpin, dimer Tm)
2. Expected reaction performance — product size, any mispriming risk or excluded region
3. Optimization tips — recommended polymerase, annealing temperature, extension time

Use plain scientific language. Reference specific numbers from the data."""

_BATCH_PCR_SYSTEM = """You are a PCR primer design expert with 15+ years of experience.
You have designed primers for multiple templates. Write a concise expert explanation
(150–200 words, 3 paragraphs) covering:
1. Overview — how many primer pairs, Tm range, any outliers
2. Quality highlights — note any pairs with concerning secondary structures or Tm spread
3. Multiplexing / practical tips — pooling, annealing temperature compromise if needed

Use plain scientific language. Reference specific numbers from the data."""

_INV_GIB_SYSTEM = """You are an expert molecular biologist specializing in Gibson Assembly from existing DNA inventory.
Write a concise expert explanation (150–200 words, 3 paragraphs) covering:
1. Coverage — what fraction of the target is covered by inventory fragments vs. synthesis gaps. If pre-design routing info is provided, briefly explain why this workflow was selected.
2. Fragment strategy — orientation of fragments, any wrap-around hits, source plasmids
3. Next steps — ordering synthesis fragments (if needed), Gibson assembly conditions, verification

Use plain scientific language. Reference specific numbers from the data (target length, coverage %, number of fragments)."""


async def generate_explanation(
    intent: str,
    result: Dict[str, Any],
    user_message: str = "",
    predesign_context: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """
    Generate an expert explanation of the design result using OpenAI.
    Returns None if OpenAI is unavailable or the call fails.

    Args:
        intent: The workflow intent (gibson_design, pcr_design, etc.)
        result: The design result dictionary
        user_message: The original user message
        predesign_context: Optional pre-design routing context including:
            - selected_workflow: The chosen workflow method
            - target_length: Target plasmid length
            - num_parts: Number of parts/fragments
            - num_compatible_workflows: Number of compatible alternatives
            - selected_cost, selected_time_days, selected_risk: Metrics for chosen workflow
            - alternative_workflows: List of alternative workflows evaluated
            - parts_summary: Summary of resolved parts
    """
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        return None

    try:
        from openai import AsyncOpenAI

        client = AsyncOpenAI(api_key=api_key)
        system_prompt, user_content = _build_prompt(intent, result, user_message, predesign_context)

        response = await client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_content},
            ],
            max_tokens=600,
        )

        return response.choices[0].message.content

    except Exception:
        return None


def _render_routing_audit(predesign_context: Optional[Dict[str, Any]]) -> str:
    """Render the target_from_inventory routing audit as a text block the
    LLM prompt can consume. Safe to call with None — returns empty string."""
    if not predesign_context:
        return ""
    ra = predesign_context.get("routing_audit") or {}
    if not ra:
        return ""
    chosen_wf = ra.get("chosen_workflow") or ""
    chosen_rat = ra.get("chosen_rationale") or ""
    score = ra.get("chosen_score")
    work = ra.get("chosen_work_estimate")
    succ = ra.get("chosen_success_estimate")
    val_mode = ra.get("chosen_validation_mode") or ""
    score_str = f"{score:.3f}" if isinstance(score, (int, float)) else "?"
    succ_str = f"{succ:.2f}" if isinstance(succ, (int, float)) else "?"
    out = [
        "",
        "METHOD SELECTION CONTEXT (from target_from_inventory router):",
        f"- Router pick: {chosen_wf} (score={score_str}, success_estimate={succ_str}, "
        f"work_estimate={work}, validation={val_mode})",
        f"- Router rationale: {chosen_rat}",
    ]
    alts = ra.get("alternatives") or []
    feas = [a for a in alts if a.get("feasible")]
    if feas:
        out.append("- Alternatives considered (feasible):")
        for a in feas[:8]:
            _s = a.get("score")
            _ss = f"{_s:.3f}" if isinstance(_s, (int, float)) else "?"
            out.append(
                f"  * {a.get('workflow')} (score={_ss}): {a.get('rationale', '')}"
            )
    rejected = [a for a in alts if not a.get("feasible")]
    if rejected:
        out.append(
            "- Rejected workflows: "
            + ", ".join(a.get("workflow", "?") for a in rejected[:10])
        )
    return "\n".join(out)


def _build_prompt(
    intent: str,
    result: Dict[str, Any],
    user_message: str,
    predesign_context: Optional[Dict[str, Any]] = None,
) -> tuple[str, str]:
    """Build (system_prompt, user_content) tuned to each intent."""

    if intent == "gibson_design":
        return _GIBSON_SYSTEM, _gibson_user(result, user_message, predesign_context)

    if intent == "pcr_design":
        return _PCR_SYSTEM, _pcr_user(result, user_message)

    if intent == "multi_pcr_design":
        if "results" in result:
            return _BATCH_PCR_SYSTEM, _batch_pcr_user(result, user_message)
        return _PCR_SYSTEM, _pcr_user(result, user_message)

    if intent == "inv_gib":
        return _INV_GIB_SYSTEM, _inv_gib_user(result, user_message, predesign_context)

    if intent == "sdm_design":
        return _GENERIC_CLONING_SYSTEM, _generic_cloning_user(
            intent="sdm_design", result=result, user_message=user_message,
            predesign_context=predesign_context,
        )

    if intent == "golden_gate_primer_design":
        return _GENERIC_CLONING_SYSTEM, _generic_cloning_user(
            intent="golden_gate_primer_design", result=result, user_message=user_message,
            predesign_context=predesign_context,
        )

    # Generic fallback — still pass predesign_context via the generic builder
    # so any intent with routing data gets the method-selection section.
    return _GENERIC_CLONING_SYSTEM, _generic_cloning_user(
        intent=intent, result=result, user_message=user_message,
        predesign_context=predesign_context,
    )


def _gibson_user(result: Dict[str, Any], user_message: str, predesign_context: Optional[Dict[str, Any]] = None) -> str:
    primers = result.get("primers_by_fragment", [])
    junctions = result.get("junctions", [])
    viz = result.get("viz", {}) or {}
    assembled_len = len(viz.get("sequence", ""))

    # Collect primer details
    primer_lines = []
    for p in primers:
        frag = p.get("fragment", "?")
        if p.get("needs_primers"):
            fwd_tm = p.get("forward_anneal_tm")
            rev_tm = p.get("reverse_anneal_tm")
            fwd_len = len(p.get("forward_primer", ""))
            rev_len = len(p.get("reverse_primer", ""))
            fwd_score = p.get("forward_extension_total_score")
            rev_score = p.get("reverse_extension_total_score")
            fwd_warn = p.get("forward_extension_warnings", "")
            rev_warn = p.get("reverse_extension_warnings", "")
            fwd_hp = p.get("forward_extension_hairpin_dg")
            rev_hp = p.get("reverse_extension_hairpin_dg")
            line = (
                f"  {frag}: Fwd {fwd_len} bp (anneal Tm {fwd_tm:.1f}°C, score {fwd_score:.0f}"
                if fwd_tm and fwd_score else f"  {frag}: Fwd {fwd_len} bp"
            )
            if fwd_hp is not None and fwd_hp < -3:
                line += f", hairpin dG={fwd_hp:.1f} kcal/mol ⚠"
            if fwd_warn:
                line += f" [{fwd_warn}]"
            line += f" | Rev {rev_len} bp"
            if rev_tm:
                line += f" (anneal Tm {rev_tm:.1f}°C"
                if rev_score:
                    line += f", score {rev_score:.0f}"
                line += ")"
            if rev_hp is not None and rev_hp < -3:
                line += f", hairpin dG={rev_hp:.1f} kcal/mol ⚠"
            if rev_warn:
                line += f" [{rev_warn}]"
            primer_lines.append(line)
        else:
            primer_lines.append(f"  {frag}: no new primers needed ({p.get('reason', '')})")

    # Collect junction details
    junction_lines = []
    for j in junctions:
        ov_seq = j.get("overlap_sequence", "")
        ov_len = j.get("overlap_length", "?")
        ov_tm = j.get("overlap_tm")
        score = j.get("overlap_score")
        uniq = j.get("overlap_uniqueness_score")
        warn = j.get("overlap_warnings", "")
        src = j.get("source", "?")
        line = f"  {j.get('from')}→{j.get('to')}: {ov_len} bp overlap ({src})"
        if ov_tm:
            line += f", Tm {ov_tm:.1f}°C"
        if score is not None:
            line += f", score {score:.0f}"
            if score < 60:
                line += " ⚠ LOW"
        if uniq is not None and uniq < 70:
            line += f", uniqueness {uniq:.0f} ⚠ LOW"
        if warn:
            line += f" [{warn}]"
        junction_lines.append(line)

    assembly = result.get("assembly", "circular")

    # Add pre-design routing context if available
    routing_info = ""
    if predesign_context:
        target_len = predesign_context.get("target_length", "?")
        num_parts = predesign_context.get("num_parts", len(primers))
        num_compatible = predesign_context.get("num_compatible_workflows", "?")
        cost = predesign_context.get("selected_cost", "?")
        alternatives = predesign_context.get("alternative_workflows", [])

        routing_info = (
            f"\nPRE-DESIGN ROUTING:\n"
            f"- Target: {target_len} bp, {predesign_context.get('target_topology', 'circular')}\n"
            f"- Parts resolved: {num_parts}\n"
            f"- Compatible workflows evaluated: {num_compatible}\n"
            f"- Selected: Gibson Assembly (estimated cost ${cost:.0f})\n"
        )

        if alternatives:
            alt_names = [a["method"] for a in alternatives[:2] if a.get("compatible")]
            if alt_names:
                routing_info += f"- Alternative workflows: {', '.join(alt_names)}\n"

    return (
        f"User request: {user_message or '(not provided)'}\n\n"
        f"DESIGN SUMMARY:\n"
        f"- Assembly: {assembly}, {len(primers)} fragment(s), {len(junctions)} junction(s)\n"
        f"- Assembled construct: {assembled_len} bp{routing_info}\n\n"
        f"PRIMERS:\n" + "\n".join(primer_lines) + "\n\n"
        f"JUNCTIONS:\n" + "\n".join(junction_lines)
    )


def _pcr_user(result: Dict[str, Any], user_message: str) -> str:
    left_tm = result.get("left_tm")
    right_tm = result.get("right_tm")
    product = result.get("product_size", "?")
    left_primer = result.get("left_primer", "")
    right_primer = result.get("right_primer", "")
    excluded = result.get("excluded_region") or {}
    l_scores = result.get("left_scores") or {}
    r_scores = result.get("right_scores") or {}
    left_mis = result.get("left_mispriming_sites", [])
    right_mis = result.get("right_mispriming_sites", [])

    lines = [
        f"User request: {user_message or '(not provided)'}",
        "",
        "DESIGN RESULT:",
        f"- Left primer ({len(left_primer)} bp): {left_primer}",
        f"  Tm: {left_tm:.1f}°C" if left_tm else "  Tm: N/A",
        f"  Hairpin Tm: {l_scores.get('hairpin_th'):.1f}°C" if l_scores.get('hairpin_th') else "",
        f"  Self-dimer Tm: {l_scores.get('any_th'):.1f}°C" if l_scores.get('any_th') else "",
        f"  Mispriming sites: {len(left_mis)}" if left_mis else "",
        f"- Right primer ({len(right_primer)} bp): {right_primer}",
        f"  Tm: {right_tm:.1f}°C" if right_tm else "  Tm: N/A",
        f"  Hairpin Tm: {r_scores.get('hairpin_th'):.1f}°C" if r_scores.get('hairpin_th') else "",
        f"  Self-dimer Tm: {r_scores.get('any_th'):.1f}°C" if r_scores.get('any_th') else "",
        f"  Mispriming sites: {len(right_mis)}" if right_mis else "",
        f"- Product size: {product} bp",
    ]
    if excluded.get("length"):
        lines.append(f"- Excluded region: positions {excluded['start']}–{excluded['start'] + excluded['length']} ({excluded['length']} bp)")

    return "\n".join(l for l in lines if l)


def _batch_pcr_user(result: Dict[str, Any], user_message: str) -> str:
    results = result.get("results", [])
    tms = [(r.get("left_tm"), r.get("right_tm")) for r in results]
    all_tms = [t for pair in tms for t in pair if t is not None]
    tm_min = min(all_tms) if all_tms else None
    tm_max = max(all_tms) if all_tms else None

    lines = [
        f"User request: {user_message or '(not provided)'}",
        "",
        f"BATCH PCR — {len(results)} template(s):",
    ]
    for r in results:
        name = r.get("template_name", f"Template_{r.get('template_index', 0) + 1}")
        lt = r.get("left_tm")
        rt = r.get("right_tm")
        ps = r.get("product_size", "?")
        lt_str = f"{lt:.1f}°C" if lt else "N/A"
        rt_str = f"{rt:.1f}°C" if rt else "N/A"
        lines.append(f"  {name}: left Tm {lt_str}, right Tm {rt_str}, product {ps} bp")

    if tm_min is not None:
        lines.append(f"\nTm range across all primers: {tm_min:.1f}°C – {tm_max:.1f}°C")

    return "\n".join(lines)


def _inv_gib_user(result: Dict[str, Any], user_message: str, predesign_context: Optional[Dict[str, Any]] = None) -> str:
    summary = result.get("inv_gib_summary", {})
    fragments = result.get("fragments_in", [])
    target_len = summary.get("target_len", 0)
    covered = summary.get("covered_bp", 0)
    pct = covered / target_len * 100 if target_len else 0
    n_inv = summary.get("emitted_inventory_fragments", 0)
    n_synth = summary.get("synth_gap_count", 0)

    frag_lines = []
    for frag in fragments:
        if not isinstance(frag, dict):
            continue
        name = frag.get("name", "?")
        start = frag.get("target_start", frag.get("start", "?"))
        end = frag.get("target_end", frag.get("end", "?"))
        length = frag.get("length_bp", "?")
        orient = frag.get("source_orientation", "+")
        src = frag.get("source_inventory", "")
        frag_lines.append(f"  {name}: {start}–{end} ({length} bp, {orient}){', from ' + src if src else ''}")

    # Add pre-design routing context if available
    routing_info = ""
    if predesign_context:
        num_compatible = predesign_context.get("num_compatible_workflows", "?")
        cost = predesign_context.get("selected_cost", "?")
        ra = predesign_context.get("routing_audit") or {}

        routing_info = (
            f"\nPRE-DESIGN ROUTING:\n"
            f"- Compatible workflows evaluated: {num_compatible}\n"
            f"- Selected: Inventory-based Gibson Assembly (estimated cost ${cost:.0f})\n"
        )

        # Full target_from_inventory routing audit, when present, so the
        # reply can cite the actual workflow decision + rationale + score
        # table rather than guessing at design intent.
        if ra:
            chosen_rat = ra.get("chosen_rationale") or ""
            score = ra.get("chosen_score")
            work = ra.get("chosen_work_estimate")
            succ = ra.get("chosen_success_estimate")
            val_mode = ra.get("chosen_validation_mode") or ""
            chosen_wf = ra.get("chosen_workflow") or ""
            score_str = f"{score:.3f}" if isinstance(score, (int, float)) else "?"
            succ_str = f"{succ:.2f}" if isinstance(succ, (int, float)) else "?"
            routing_info += (
                f"- Router decision: {chosen_wf} "
                f"(score={score_str}, success_estimate={succ_str}, "
                f"work_estimate={work}, validation={val_mode})\n"
                f"- Router rationale: {chosen_rat}\n"
            )
            alts = ra.get("alternatives") or []
            feas_alts = [a for a in alts if a.get("feasible")]
            if feas_alts:
                routing_info += "- Alternatives considered (feasible):\n"
                for a in feas_alts[:6]:
                    _s = a.get("score")
                    _ss = f"{_s:.3f}" if isinstance(_s, (int, float)) else "?"
                    routing_info += (
                        f"  * {a.get('workflow')} (score={_ss}): "
                        f"{a.get('rationale', '')}\n"
                    )
            rejected = [a for a in alts if not a.get("feasible")]
            if rejected:
                routing_info += (
                    "- Rejected workflows: "
                    + ", ".join(a.get("workflow", "?") for a in rejected[:8])
                    + "\n"
                )

    return (
        f"User request: {user_message or '(not provided)'}\n\n"
        f"TARGET: {target_len} bp\n"
        f"COVERAGE: {covered} / {target_len} bp ({pct:.0f}%)\n"
        f"INVENTORY FRAGMENTS: {n_inv}\n"
        f"SYNTHESIS GAPS: {n_synth}{routing_info}\n\n"
        f"FRAGMENTS:\n" + "\n".join(frag_lines)
    )


def _summarize_result(intent: str, result: Dict[str, Any]) -> str:
    """Kept for backward compatibility — returns brief summary string."""
    if intent == "gibson_design":
        junctions = result.get("junctions", [])
        primers = result.get("primers_by_fragment", [])
        n_primers = sum(1 for p in primers if p.get("needs_primers"))
        viz = result.get("viz", {}) or {}
        seq_len = len(viz.get("sequence", ""))
        return (
            f"Gibson assembly: {len(junctions)} junction(s), "
            f"{n_primers} fragment(s) need primers, "
            f"assembled construct {seq_len} bp"
        )

    if intent in ("pcr_design", "multi_pcr_design") and "left_tm" in result:
        left_tm = result.get("left_tm")
        right_tm = result.get("right_tm")
        product = result.get("product_size")
        left_str = f"{left_tm:.1f}°C" if left_tm is not None else "N/A"
        right_str = f"{right_tm:.1f}°C" if right_tm is not None else "N/A"
        return f"PCR: left Tm={left_str}, right Tm={right_str}, product={product} bp"

    if intent == "multi_pcr_design" and "results" in result:
        count = result.get("count", len(result.get("results", [])))
        tms = [r.get("left_tm") for r in result.get("results", []) if r.get("left_tm") is not None]
        avg_tm = sum(tms) / len(tms) if tms else None
        avg_str = f", avg left Tm={avg_tm:.1f}°C" if avg_tm else ""
        return f"Batch PCR: {count} primer pair(s) designed{avg_str}"

    if intent == "inv_gib":
        summary = result.get("inv_gib_summary", {})
        covered = summary.get("covered_bp", 0)
        total = summary.get("target_len", 1)
        pct = covered / total * 100 if total else 0
        frags = summary.get("emitted_inventory_fragments", 0)
        synth = summary.get("synth_gap_count", 0)
        return (
            f"Inventory Gibson: {covered} / {total} bp covered ({pct:.0f}%), "
            f"{frags} inventory fragment(s), {synth} synthesis gap(s)"
        )

    if intent == "annotate_gb":
        count = result.get("annotation_count", 0)
        return f"Plasmid annotation: {count} feature(s) identified"

    return str(result)[:300]

_GENERIC_CLONING_SYSTEM = (
    "You are an expert molecular biologist writing a concise design rationale "
    "(2-3 short paragraphs). When a METHOD SELECTION CONTEXT is provided, cite "
    "which alternative workflows the router evaluated and why the chosen one "
    "was preferred. Never claim the inventory is empty if fragments or coverage "
    "numbers are given. Be specific: name the cloning method, cite the bp length, "
    "and reference actual features where the context supplies them."
)


def _generic_cloning_user(
    intent: str,
    result: Dict[str, Any],
    user_message: str,
    predesign_context: Optional[Dict[str, Any]] = None,
) -> str:
    """Cloning-agnostic LLM prompt builder that tacks the routing audit onto
    a compact summary of whatever result was returned. Used for sdm_design,
    golden_gate_primer_design, and as the generic fallback."""
    summary = _summarize_result(intent, result)
    routing_info = _render_routing_audit(predesign_context)
    return (
        f"User request: {user_message or '(not provided)'}\n\n"
        f"INTENT: {intent}\n"
        f"RESULT SUMMARY: {summary}\n"
        f"{routing_info}\n\n"
        "Explain the design in 2-3 short paragraphs for a scientist, "
        "leading with method selection if routing context is present."
    )

