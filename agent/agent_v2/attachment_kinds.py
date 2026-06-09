"""Per-request store mapping attachment_id -> classification + cached annotation.

AttachmentRegistry (v1) is shared with the v1 production agent and we don't
want to change its public contract. Instead, the v2 router runs the file_kind
classifier at upload time, stashes the result here, and (for genomic-kind
files) keeps the raw GenBank text so the annotator can be invoked lazily
when a subagent or tool needs it.

Process-local. No Redis dependency. Cleared per test via clear_all_kinds().
"""
from __future__ import annotations

import threading
from dataclasses import dataclass, field
from typing import Optional

from agent_v2.file_kind import FileKind


@dataclass
class CachedAttachment:
    kind: FileKind
    gb_text: Optional[str] = None
    # Lazily built GenomicAnnotation. Untyped here to avoid circular imports;
    # consumers cast via agent_v2.genomic_annotator.GenomicAnnotation.
    _genomic_annotation: object = field(default=None, repr=False)


_CACHE: dict[str, CachedAttachment] = {}
_LOCK = threading.Lock()


def stash_kind(attachment_id: str, kind: FileKind, *, gb_text: Optional[str] = None) -> None:
    with _LOCK:
        _CACHE[attachment_id] = CachedAttachment(kind=kind, gb_text=gb_text)


def get_kind(attachment_id: str) -> Optional[FileKind]:
    with _LOCK:
        c = _CACHE.get(attachment_id)
    return c.kind if c is not None else None


def kind_str(attachment_id: str, default: str = "plasmid") -> str:
    """Return the kind name ('plasmid' / 'genomic' / 'unknown'), or the default."""
    fk = get_kind(attachment_id)
    return fk.kind if fk is not None else default


def get_genomic_annotation(attachment_id: str):
    """Lazy build of the GenomicAnnotation for a genomic-kind attachment.

    Returns None if the attachment isn't genomic kind or has no stashed gb_text.
    """
    with _LOCK:
        c = _CACHE.get(attachment_id)
    if c is None or c.kind.kind != "genomic" or not c.gb_text:
        return None
    if c._genomic_annotation is None:
        from agent_v2.genomic_annotator import annotate_genomic_gb
        c._genomic_annotation = annotate_genomic_gb(c.gb_text)
    return c._genomic_annotation


def clear_all_kinds() -> None:
    """Test helper. Resets the process-local map."""
    with _LOCK:
        _CACHE.clear()

def set_annotation(attachment_id: str, annotation) -> None:
    """Eagerly stash a pre-built GenomicAnnotation in the cache so the next
    call to get_genomic_annotation returns immediately. Called from
    find_genomic_record_tool after the NCBI .gb is registered and the
    annotation has already been built once for the digest."""
    with _LOCK:
        c = _CACHE.get(attachment_id)
        if c is not None:
            c._genomic_annotation = annotation
