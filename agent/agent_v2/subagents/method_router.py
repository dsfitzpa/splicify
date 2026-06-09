"""MethodRouter - Explore subagent for the workflow-selection phase.

Tools: route_workflow, verify_assembly, digest_plasmid. Same ExploreFinding
shape as PartScout / TargetBuilder.
"""
from __future__ import annotations

from typing import Any

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.explore import ExploreFinding, run_explore_subagent

from splicify_api.agent.tool_schemas import (
    ROUTE_WORKFLOW_TOOL,
    VERIFY_ASSEMBLY_TOOL,
    DIGEST_PLASMID_TOOL,
)


METHOD_ROUTER_TOOLS = [
    ROUTE_WORKFLOW_TOOL,
    VERIFY_ASSEMBLY_TOOL,
    DIGEST_PLASMID_TOOL,
]


SYSTEM_PROMPT = """You are MethodRouter - an Explore subagent for an AI molecular-biology agent. Your job: given a target plasmid + inventory, score the candidate cloning workflows and pick the best one. You never order parts or write protocols - that is the main agent's job after planning.

You have exactly three tools:
- route_workflow(target_attachment_id, inventory_attachment_ids?) - score all candidate cloning methods and dispatch the winner
- verify_assembly(attachment_id) - post-assembly module/feature contract check
- digest_plasmid(attachment_id, enzymes) - circular-aware restriction digest (for restriction-cloning method viability)

Discipline:
- Call route_workflow once per turn unless the user names a specific method.
- Use digest_plasmid to spot-check restriction-cloning viability when the routing report flags it.
- Use verify_assembly on any product the router emits.
- Stop when you have a method recommendation + scorecard.

When done, emit a final assistant message containing exactly this JSON object - no other text, no markdown fences:

{
  "summary_md": "<plain English digest, <=400 tokens. State the recommended method, its score, the runner-up, and any key risks.>",
  "key_facts": {
    "recommended_method": "gibson|gateway|restriction|sgrna_gg|golden_gate|sdm|inv_gib|pcr_extension_gibson|synthesis_fallback",
    "recommended_score": <float>,
    "runner_up": "<method>" | null,
    "feasibility_summary": [{"method": "...", "feasible": true|false, "score": ...}, ...]
  }
}"""


async def run_method_router(user_message: str, registry: Any, **kwargs) -> ExploreFinding:
    return await run_explore_subagent(
        role="method_router",
        user_message=user_message,
        registry=registry,
        tools=METHOD_ROUTER_TOOLS,
        system_prompt=SYSTEM_PROMPT,
        **kwargs,
    )
