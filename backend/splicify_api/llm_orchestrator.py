"""
LLMOrchestrator — interface for the (future) LLM that reasons over a
PlasmidSpec, the assembled target, and the annotation diff to suggest
edit operations a deterministic step couldn't resolve on its own.

v1 leaves this as a no-op pass-through. The deterministic edit_ops path in
describe_plasmid_handler.py handles every case where the missing module is
KB-resolvable; the orchestrator is reserved for cases where:
  - a module is named but has no KB entry
  - a module is described in prose ("a strong inducible mammalian promoter")
    rather than named, and the deterministic role-keyword path under-resolves
  - the assembled target satisfies the spec by role but not by intent
    (e.g. correct module type, wrong host compatibility)

Two public call sites are wired with a guard:
  - describe_plasmid_handler.describe_plasmid (after deterministic edit_ops)
  - chat.py _execute_unified_predesign (after spec_diff)
Both check `is_enabled()` and skip the call otherwise.
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Protocol

from .plasmid_spec import PlasmidSpec

logger = logging.getLogger("orchestrator")


@dataclass
class OrchestratorEdit:
    """An edit the orchestrator proposes. Mirrors describe_plasmid_handler.EditOp
    but is decoupled so the orchestrator package does not import from a
    specific handler.
    """
    op:        str  # "add" | "remove" | "replace"
    target:    str
    role:      Optional[str] = None
    sequence:  Optional[str] = None
    strategy:  str = "synthesis"  # "primer_tail" | "synthesis"
    rationale: str = ""
    confidence: float = 0.0


class LLMOrchestrator(Protocol):
    """Protocol the eventual LLM-backed orchestrator will satisfy."""

    async def review(
        self,
        *,
        spec: PlasmidSpec,
        target_modules: List[Dict[str, Any]],
        spec_diff: Dict[str, Any],
        unresolved_edits: List[Dict[str, Any]],
    ) -> List[OrchestratorEdit]:
        ...


# ---------------------------------------------------------------------------
# v1 implementation: no-op
# ---------------------------------------------------------------------------
class NoOpOrchestrator:
    async def review(
        self,
        *,
        spec: PlasmidSpec,
        target_modules: List[Dict[str, Any]],
        spec_diff: Dict[str, Any],
        unresolved_edits: List[Dict[str, Any]],
    ) -> List[OrchestratorEdit]:
        if unresolved_edits:
            logger.info(
                "orchestrator disabled — leaving %d unresolved edits for the user",
                len(unresolved_edits),
            )
        return []


# ---------------------------------------------------------------------------
# Public entry points (used by call sites — keep stable across phases)
# ---------------------------------------------------------------------------
def is_enabled() -> bool:
    """The orchestrator is opt-in via PLASMID_LLM_ORCHESTRATOR=1.
    v1 ships disabled; flip this when the LLM backend is implemented."""
    return os.getenv("PLASMID_LLM_ORCHESTRATOR", "").strip() in ("1", "true", "yes")


_ORCHESTRATOR: Optional[LLMOrchestrator] = None


def get_orchestrator() -> LLMOrchestrator:
    global _ORCHESTRATOR
    if _ORCHESTRATOR is None:
        _ORCHESTRATOR = NoOpOrchestrator()
    return _ORCHESTRATOR


async def review_design(
    *,
    spec: PlasmidSpec,
    target_modules: List[Dict[str, Any]],
    spec_diff: Dict[str, Any],
    unresolved_edits: Optional[List[Dict[str, Any]]] = None,
) -> List[OrchestratorEdit]:
    """Optional gap-filler. Returns an empty list when disabled."""
    if not is_enabled():
        return []
    try:
        return await get_orchestrator().review(
            spec=spec,
            target_modules=target_modules,
            spec_diff=spec_diff,
            unresolved_edits=unresolved_edits or [],
        )
    except Exception as exc:  # pragma: no cover — best effort
        logger.warning("orchestrator review failed: %s", exc)
        return []
