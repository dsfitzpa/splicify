"""Splicify — bioinformatics tools for molecular biology cloning and CRISPR
guide design.

Public Python API. Anything not re-exported from this module is internal
to ``splicify_api`` and may change without notice between releases.

Import style::

    from splicify_api import AttachmentRegistry, annotate_cached, design_guides
    # or, for the agent runtime:
    from splicify_api import (
        AttachmentRegistry,
        AIPLASMIDDESIGN_TOOLS,
        dispatch_tool,
        annotate_llm_cached,
    )
"""

# Attachment registry + agent dispatcher (used by every agent integration).
from splicify_api.agent.agent_tools import (
    AttachmentRegistry,
    dispatch_tool,
    extract_seq_from_genbank,
    tool_lookup_kb_part,
    _gb_for_attachment,
)

# Anthropic tool schemas.
from splicify_api.agent.tool_schemas import (
    AIPLASMIDDESIGN_TOOLS,
    ANNOTATE_ATTACHMENT_TOOL,
)

# Annotation cache (SHA1-keyed; hierarchy + LLM endpoints).
from splicify_api.annotation_cache import (
    annotate_cached,
    annotate_llm_cached,
)

# External lookups (NCBI fetch + rate gate live inside).
from splicify_api.external_search import (
    search_ncbi_gene,
    fetch_ncbi_genbank,
)

# CRISPR guide design (Doench 2014 + heuristics).
from splicify_api.guide_designer import design_guides

# Primer design (Primer3 + application-aware Sanger / Illumina paths).
from splicify_api.pcr import (
    PrimerRequest,
    design_primers,
)

# Prime editing (PE3 XGBoost port of easy_prime).
from splicify_api.pegrna_designer import design_pegrnas

__version__ = "0.1.0"

__all__ = [
    "AttachmentRegistry",
    "dispatch_tool",
    "extract_seq_from_genbank",
    "tool_lookup_kb_part",
    "_gb_for_attachment",
    "AIPLASMIDDESIGN_TOOLS",
    "ANNOTATE_ATTACHMENT_TOOL",
    "annotate_cached",
    "annotate_llm_cached",
    "search_ncbi_gene",
    "fetch_ncbi_genbank",
    "design_guides",
    "PrimerRequest",
    "design_primers",
    "design_pegrnas",
    "__version__",
]
