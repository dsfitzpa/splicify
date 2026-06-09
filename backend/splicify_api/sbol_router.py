"""
FastAPI router for SBOL3 export / import.

Exposes two endpoints under /sbol:

- POST /sbol/export   – in: annotation payload (same shape as the
  `annotate_sequence_llm` response). Out: SBOL3 document as a string +
  MIME-typed file download.

- POST /sbol/import   – in: a raw SBOL3 document (string or uploaded file).
  Out: flat annotation payload (features + interactions + sequence) that the
  frontend viewer / backend workflows can consume directly.

Design goals:
- Don't duplicate annotation logic. The export endpoint trusts its caller to
  provide enriched annotations + interactions, matching what the main pipeline
  already returns.
- Keep the import endpoint structural — just convert SBOL3 → feature list;
  re-running the pipeline is the caller's job if they want richer output.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException, UploadFile, File, Form
from pydantic import BaseModel, Field

from .sbol_io import (
    DEFAULT_NAMESPACE,
    export_annotation_to_sbol3,
    import_sbol3,
)

router = APIRouter(prefix="/sbol", tags=["sbol"])

# --------------------------------------------------------------------------- #
# Request / response models
# --------------------------------------------------------------------------- #


class SBOLExportRequest(BaseModel):
    sequence: str = Field(..., description="Raw DNA sequence (ACGT)")
    annotations: List[Dict[str, Any]] = Field(
        default_factory=list,
        description="hierarchical_annotations from annotate_sequence_llm",
    )
    interactions: List[Dict[str, Any]] = Field(
        default_factory=list,
        description="interactions payload from annotate_sequence_llm",
    )
    plasmid_name: str = Field("plasmid", description="Display name for the Component")
    circular: bool = True
    namespace: str = DEFAULT_NAMESPACE
    file_format: str = Field(
        "sorted nt",
        description="One of: 'sorted nt', 'ttl', 'xml', 'json-ld', 'nt11'",
    )


class SBOLExportResponse(BaseModel):
    ok: bool
    document: str
    file_format: str
    mime_type: str
    namespace: str
    component_count: int
    interaction_count: int


class SBOLImportResponse(BaseModel):
    ok: bool
    plasmid_name: str
    sequence: str
    annotations: List[Dict[str, Any]]
    interactions: List[Dict[str, Any]]
    annotation_count: int
    interaction_count: int


_MIME_FOR_FORMAT = {
    "sorted nt": "application/n-triples",
    "nt11": "application/n-triples",
    "ttl": "text/turtle",
    "xml": "application/rdf+xml",
    "json-ld": "application/ld+json",
}

# --------------------------------------------------------------------------- #
# Endpoints
# --------------------------------------------------------------------------- #


@router.post("/export", response_model=SBOLExportResponse)
async def export_sbol(request: SBOLExportRequest) -> SBOLExportResponse:
    """Serialize an annotation payload to an SBOL3 document."""
    try:
        doc = export_annotation_to_sbol3(
            sequence=request.sequence,
            annotations=request.annotations,
            interactions=request.interactions,
            plasmid_name=request.plasmid_name,
            circular=request.circular,
            namespace=request.namespace,
        )
        serialized = doc.write_string(request.file_format)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=400, detail=f"SBOL3 export failed: {e}") from e

    return SBOLExportResponse(
        ok=True,
        document=serialized,
        file_format=request.file_format,
        mime_type=_MIME_FOR_FORMAT.get(request.file_format, "text/plain"),
        namespace=request.namespace,
        component_count=len(request.annotations),
        interaction_count=len(request.interactions),
    )


@router.post("/import", response_model=SBOLImportResponse)
async def import_sbol(
    document: Optional[str] = Form(None),
    file_format: Optional[str] = Form(None),
    file: Optional[UploadFile] = File(None),
) -> SBOLImportResponse:
    """Parse an SBOL3 document (inline or uploaded) into a flat annotation payload."""
    if document is None and file is None:
        raise HTTPException(status_code=400, detail="Provide either `document` or `file`.")

    if file is not None:
        raw = (await file.read()).decode("utf-8", errors="replace")
    else:
        raw = document or ""

    try:
        parsed = import_sbol3(raw, file_format=file_format)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=400, detail=f"SBOL3 import failed: {e}") from e

    return SBOLImportResponse(
        ok=True,
        plasmid_name=parsed.get("plasmid_name", ""),
        sequence=parsed.get("sequence", ""),
        annotations=parsed.get("annotations", []),
        interactions=parsed.get("interactions", []),
        annotation_count=len(parsed.get("annotations", [])),
        interaction_count=len(parsed.get("interactions", [])),
    )


__all__ = ["router"]
