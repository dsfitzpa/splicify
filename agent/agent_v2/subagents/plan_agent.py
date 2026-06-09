"""Plan subagent.

Receives ONLY the three ExploreFinding digests + the user prompt. Has no
tools, no parent context, no raw DNA, no full annotation lists. Emits a
markdown todo list (`plan.md`) that the Main Agent will read + edit as it
crosses items off during the execute phase.

Mirrors the Claude Code blog-post pattern:
> "The Plan Agent did not carry all the context from the main agent nor the
> Explore subagents ... it only contains the summarization of the Explore
> subagents' findings."
"""
from __future__ import annotations

import os
import pathlib
from dataclasses import dataclass
from typing import Any, Optional, Sequence

from agent_v2.explore import ExploreFinding


@dataclass
class PlanResult:
    plan_md: str
    plan_path: Optional[pathlib.Path]
    n_steps: int


SYSTEM_PROMPT = """You are the Plan agent for an AI molecular-biology agent. You receive:
1. The user's original prompt.
2. Three digested ExploreFinding summaries (PartScout, TargetBuilder, MethodRouter).

You DO NOT have raw DNA. You DO NOT have full annotation lists. You DO NOT have tools.

Your job: emit a SHORT markdown todo list that the Main Agent will execute. Each item must be a concrete tool call or output emission. The Main Agent has these tools available:
- find_external_part, annotate_attachment, find_features, find_cassette_for, ask_plasmid, lookup_kb_part, golden_gate_assemble, digest_plasmid, find_primer_binding_sites, score_sanger_primer, route_workflow, verify_assembly, compare_to_choice, web_search
- cassette_swap (PRIMARY assembly tool — deterministic end-to-end cassette swap: auto-picks parent, runs find_cassette_for, computes bounds, calls replace_region. Use this for any named-feature swap like "replace PuroR with mCherry inside the Cas9 cassette".)
- replace_region (lower-level in-place plasmid edit — use only when you need numeric bounds explicitly, e.g. additive insertion at a known MCS position)
- emit_assembled_gb, emit_parts_order, emit_protocol, emit_workflow_trace (output emitters)
- `simulate_assembly` and `graft_parts` are ARCHIVED — never list them in the plan.

ROUTING RULES (apply these every time):

1. **Always start from a backbone plasmid.** Even when the user doesn't name one, the plan's first execute step is `find_external_part(description=..., required_features=[...])` to locate a parent. Every design uses `replace_region` against a backbone — there is no de-novo concatenation path anymore.

2. **For ANY modification to an existing plasmid, plan `replace_region`.** Cassette swaps, CDS replacements, marker substitutions, tag additions/removals, sgRNA stuffer cloning, insertions at MCS positions — all of these are `replace_region`. The tool keeps the parent intact and only swaps one region. Cassette signatures the Main agent should recognize:
   - User specifies `[A]-[X]-[B]` and parent has `[A]-[Y]-[B]` → swap `[X]` for `[Y]` (anchor shared at ≥80% nt identity).
   - User specifies a functional element (selection marker, fluorescent reporter, Pol II promoter, Pol III promoter, polyA, 2A peptide, NLS, origin) and parent has a different member of the same group in the same module → functional substitution.
   - User wants to add/remove a small tag (FLAG, HA, Myc, His, V5) on a CDS in the parent → tag edit.
   - User provides an sgRNA spacer and parent has a guide_expression_cassette stuffer → stuffer replacement.
   - User wants to add a new cassette and parent has an MCS / lacZα / inter-module gap → zero-bp insertion at the MCS via replace_region.

3. **Interchangeability groups** (use these to recognize substitutions):
   - Bacterial selection: AmpR, KanR, CmR, SpcR, TetR, GentR
   - Mammalian selection: PuroR, NeoR, HygroR, BlastR, BleoR(Zeocin), BsdR
   - Fluorescent reporter (group spans all colors AND can substitute for a mammalian selection marker when used as a selection readout): EGFP, mCherry, mOrange, BFP, YFP, mRuby, mScarlet, iRFP, etc.
   - Mammalian Pol II promoter: CMV, EF-1α, CAG, PGK, UBC, SV40 early
   - Pol III promoter: U6, H1, 7SK
   - PolyA: bGH, SV40, rbglob, hGH, SPA
   - 2A peptide: P2A, T2A, E2A, F2A
   - NLS: SV40 NLS, nucleoplasmin NLS, c-myc NLS
   - Bacterial origin: pUC, ColE1, pMB1, p15A
   The KB's `feature_class` field is the authoritative grouping when available.

4. **Insertion points are picked automatically** by replace_region when not specified — prefer an MCS (especially in lacZα for blue-white screening backbones), otherwise the largest inter-module gap.

5. **Refusals**: when no replacement rule fits cleanly, the Main agent will return `ok=false` from replace_region with diagnostics; the plan should anticipate this and include a fallback "ask user to clarify" step rather than chaining downstream emitters past a failure.

Plan format — REQUIRED step ordering for ANY modification of an existing plasmid:

1. **`find_external_part`** — locate the parent vector (skip when a backbone is already uploaded).
2. **`find_cassette_for(target_attachment_id, query=<anchor feature>)`** — REQUIRED before any replace_region. Returns the functional cassette boundaries (cassette_start, cassette_end) plus the ordered submodule list. The Main agent uses these to compute exact region_start / region_end inputs to step N. Do NOT list `find_features` here as a substitute — find_cassette_for already returns the matching features and adds the cassette context the agent needs to pick boundaries correctly.
3. **`find_features`** — only when you need a SPECIFIC sub-feature boundary that find_cassette_for didn't surface (e.g. a flanking restriction site).
4. **`lookup_kb_part`** for any NEW parts to insert. Read each response's `parent_match` field: if populated, the KB hit is redundant with the parent — do NOT pass that attachment_id to replace_region; slice the parent instead.
5. **(optional) `ask_plasmid`** — for ambiguous design decisions where rules 2-3 aren't enough (e.g. "where should I cut to add a second sgRNA cassette without disrupting Cas9 expression?"). Slow; reach for it only when find_cassette_for can't resolve.
6. **EXACTLY ONE `cassette_swap`** call. Args: `cassette_anchor` (the feature that identifies which cassette to edit, usually the CDS to keep, e.g. 'Cas9'), `keep_through` (last feature to preserve — excision starts at this feature's end + 1), `excise_through` (last feature to delete — excision ends at this feature's end), `replacement_parts` (new content only — literal sequences and KB hits without parent_match). The tool computes the exact region_start/region_end internally from find_cassette_for + the named features. Only fall back to bare `replace_region` when you need explicit numeric bounds for an additive insertion (e.g. zero-bp excision at an MCS position).
7. **`verify_assembly`**, then **`emit_assembled_gb`**, **`emit_parts_order`**, **`emit_protocol`**.
- `emit_workflow_trace` is auto-emitted by the orchestrator — do not list it.
- <=15 items. No commentary outside the heading + list.

Output structure: `## Plan` heading then a numbered Markdown checklist (`- [ ] 1. tool_name — what this call does`)."""


def plan_path_for(session_id: str, turn_id: str,
                   base: str = "/var/data/agent_v2_runs") -> pathlib.Path:
    return pathlib.Path(base) / session_id / turn_id / "plan.md"


def _format_finding(f: ExploreFinding) -> str:
    return (
        f"### {f.role}\n"
        f"{f.summary_md}\n\n"
        f"**key_facts:** `{f.key_facts}`\n"
    )


async def run_plan_agent(
    user_message: str,
    findings: Sequence[ExploreFinding],
    *,
    plan_path: Optional[pathlib.Path] = None,
    client: Any = None,
    model: Optional[str] = None,
    system_prompt: Optional[str] = None,
) -> PlanResult:
    if client is None:
        import anthropic
        client = anthropic.Anthropic()
    model = model or os.getenv("AGENT_MODEL", "claude-sonnet-4-6")
    system = system_prompt if system_prompt is not None else SYSTEM_PROMPT

    user_block = (
        f"User prompt:\n{user_message}\n\n"
        "Explore findings (digested):\n\n"
        + "\n\n".join(_format_finding(f) for f in findings)
    )

    response = client.messages.create(
        model=model,
        max_tokens=2048,
        system=system,
        messages=[{"role": "user", "content": user_block}],
    )
    text = "".join(getattr(b, "text", "") for b in response.content
                   if getattr(b, "type", None) == "text")
    plan_md = text.strip() or "## Plan\n- [ ] 1. (Plan agent emitted no content)"

    if plan_path is not None:
        plan_path.parent.mkdir(parents=True, exist_ok=True)
        plan_path.write_text(plan_md)
    n_steps = sum(1 for line in plan_md.splitlines() if line.lstrip().startswith("- ["))

    return PlanResult(plan_md=plan_md, plan_path=plan_path, n_steps=n_steps)
