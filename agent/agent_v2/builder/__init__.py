"""Builder package: parts → virtual construct → verify → revise loop.

Designed as a parallel pipeline branch from the current cloning
pipeline. Inputs: an `IntentSpec` (required modules / interactions /
function) + an `AttachmentRegistry` of source plasmids. Outputs: a
verified `VirtualConstruct` ready to hand to a cloning-method
workflow.

Constraints:
  - Every Part carries 50 bp of upstream + downstream junction
    sequence from its source plasmid so the builder can detect
    compatibility for Gibson / Golden Gate / restriction methods.
  - Builder does NOT add new sequence. If a verifier check fails it
    can only reorder, reorient, swap parts, or fail.
  - Verifier annotates the materialized construct through
    annotate_llm_cached and runs deterministic checks against the
    IntentSpec.
"""
