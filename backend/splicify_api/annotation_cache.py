"""
Sequence-hash annotation cache.

Annotation is deterministic for a given sequence + topology, so callers that
annotate the same sequence multiple times in one request (PartResolver,
TargetPlasmidBuilder, target_from_inventory_router, describe_plasmid handler)
can hit memory instead of re-running the full pipeline.

Two annotation depths:
  - "full"            → /annotate_sequence_with_hierarchy: features + modules
                         + interactions + cloning features. Use on parts and
                         on the final response target.
  - "modules_only"    → hierarchy with cloning-feature pass disabled. Use on
                         intermediate candidate target plasmids during
                         workflow assessment.
"""
from __future__ import annotations

import hashlib
import logging
from typing import Any, Dict, Optional

logger = logging.getLogger("annotation_cache")

_CACHE: Dict[str, Dict[str, Any]] = {}
_MAX_ENTRIES = 256


def _key(sequence: str, *, circular: bool, depth: str) -> str:
    h = hashlib.sha1(sequence.upper().encode("ascii", errors="ignore")).hexdigest()
    return f"{depth}:{int(bool(circular))}:{h}"


def _evict_if_needed() -> None:
    if len(_CACHE) <= _MAX_ENTRIES:
        return
    # Drop the oldest ~25 % of entries (insertion order on dicts is preserved)
    drop = max(1, len(_CACHE) // 4)
    for k in list(_CACHE.keys())[:drop]:
        _CACHE.pop(k, None)


def get_cached(sequence: str, *, circular: bool, depth: str) -> Optional[Dict[str, Any]]:
    return _CACHE.get(_key(sequence, circular=circular, depth=depth))


def put_cached(
    sequence: str, *, circular: bool, depth: str, result: Dict[str, Any]
) -> None:
    _evict_if_needed()
    _CACHE[_key(sequence, circular=circular, depth=depth)] = result


async def annotate_cached(
    sequence: str,
    *,
    circular: bool = True,
    depth: str = "full",
) -> Dict[str, Any]:
    """Annotate a sequence with caching, via the canonical
    /annotate_sequence_llm pipeline (rule-based modules, Pol2 cassettes,
    interaction_builder, cloning features, ORF filter + per-AA
    translation strips).

    `depth`:
      - "full"          everything, including cloning_features.
      - "modules_only"  cloning_features stripped for cheaper intermediate
                        target-assessment calls (the rule-based modules
                        and translations still come along).
    """
    cached = get_cached(sequence, circular=circular, depth=depth)
    if cached is not None:
        logger.debug("annotation cache hit (%s, %d bp)", depth, len(sequence))
        return cached

    from .plannotate_router import (
        AnnotateSequenceRequest,
        annotate_sequence_llm_endpoint,
    )

    req = AnnotateSequenceRequest(
        sequence=sequence, circular=circular, detailed=True, hierarchical=True,
    )
    result = await annotate_sequence_llm_endpoint(req)

    if depth == "modules_only" and isinstance(result, dict):
        # Shallow copy so we don't mutate the LLM-endpoint cache below it.
        result = {**result}
        result.pop("cloning_features", None)
        result.pop("cut_profile", None)

    put_cached(sequence, circular=circular, depth=depth, result=result)
    logger.debug("annotation cache miss → stored (%s, %d bp)", depth, len(sequence))
    return result


async def annotate_llm_cached(
    sequence: str,
    *,
    circular: bool = True,
) -> Dict[str, Any]:
    """Backwards-compatible alias for `annotate_cached(depth="full")`.

    Both functions now route through the same canonical
    /annotate_sequence_llm pipeline; the legacy "annotate_cached uses
    annotate_sequence_with_hierarchy / annotate_llm_cached uses
    annotate_sequence_llm" split is gone.
    """
    return await annotate_cached(sequence, circular=circular, depth="full")


def clear_cache() -> None:
    _CACHE.clear()
