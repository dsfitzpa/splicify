"""
FastAPI router for the AIPlasmidDesign agent (v3).

Endpoints:
  POST /agent/chat   — accepts EITHER multipart/form-data (with file uploads)
                        OR application/json. Same logical shape as /api/chat.
                        Optional `choices_json` for MCQ benchmark scoring.
  GET  /agent/health — anthropic SDK installed? API key set?

Response shape (always JSON, even on errors):
  {
    reply:         str,
    files:         [ {fileName, dataBase64} ]   # GenBank files emitted by tool calls
    viz:           {sequence, annotations, ...} # primary product viz (for CircularPlasmidViewer)
    agent_trace:   [...],
    n_tool_calls:  int,
    error:         str | null,
    ok:            bool,
  }
"""
from __future__ import annotations

import base64
import json
import logging
import os
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, File, Form, Request, UploadFile
from fastapi.responses import JSONResponse

from .agent_loop import run_agent

logger = logging.getLogger("agent.router")

router = APIRouter(prefix="/agent", tags=["agent"])


@router.get("/health")
async def health():
    try:
        import anthropic  # noqa: F401
        sdk_ok = True
    except ImportError:
        sdk_ok = False
    return {
        "anthropic_sdk_installed": sdk_ok,
        "api_key_set":             bool(os.environ.get("ANTHROPIC_API_KEY")),
        "model":                   os.environ.get("AGENT_MODEL", "claude-opus-4-7"),
        "max_tokens":              int(os.environ.get("AGENT_MAX_TOKENS", "16384")),
    }


def _gb_for_attachment_payload(name: str, sequence: str, circular: bool) -> str:
    """Render a minimal GenBank file for an attachment (matches the format the
    agent_tools helper uses but kept inline so we don't import private state)."""
    seq = (sequence or "").lower()
    topo = "circular" if circular else "linear"
    safe_name = "".join(c for c in (name or "product") if c.isalnum() or c in "_-")[:16] or "product"
    lines = [
        f"LOCUS       {safe_name:<16} {len(seq):>8} bp    DNA     {topo:<8} SYN 01-JAN-2026",
        f"DEFINITION  {name}",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(seq)}",
        '                     /organism="synthetic DNA construct"',
        '                     /mol_type="other DNA"',
        "ORIGIN",
    ]
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60]
        groups = [chunk[j:j + 10] for j in range(0, len(chunk), 10)]
        lines.append(f"{i+1:>9} " + " ".join(groups))
    lines.append("//")
    return "\n".join(lines) + "\n"


async def _build_response_payload(result: Dict[str, Any]) -> Dict[str, Any]:
    """Convert agent run result into the frontend-friendly response shape.
    Surfaces every product attachment registered during the run as a
    downloadable GenBank file plus a primary viz block."""
    files: List[Dict[str, str]] = []
    viz: Optional[Dict[str, Any]] = None

    registry = result.get("_registry")
    if registry is not None:
        try:
            from .annotation_cache_passthrough import maybe_annotate
        except Exception:
            maybe_annotate = None
        for att in registry.items.values():
            if att.role == "product":
                gb = _gb_for_attachment_payload(att.name, att.sequence, att.circular)
                files.append({
                    "fileName":   f"{att.name}.gb",
                    "dataBase64": base64.b64encode(gb.encode("utf-8")).decode("ascii"),
                })
                if viz is None:
                    viz = {
                        "type":     "plasmid",
                        "title":    att.name,
                        "sequence": att.sequence,
                        "circular": att.circular,
                    }
    return {
        "ok":           result.get("error") is None,
        "reply":        result.get("reply", ""),
        "files":        files or None,
        "viz":          viz,
        "agent_trace":  result.get("trace", []),
        "n_tool_calls": result.get("n_tool_calls", 0),
        "error":        result.get("error"),
    }


@router.post("/chat")
async def agent_chat(request: Request):
    """Dual-input: multipart/form-data OR application/json.
    Returns JSON envelope with reply + emitted files + viz.
    """
    content_type = (request.headers.get("content-type") or "").lower()
    message:        str = ""
    session_id:     str = ""
    choices_json:   str = ""
    target_gb:      Optional[str] = None
    inventory_gbs:  List[str] = []

    try:
        if "multipart/form-data" in content_type:
            form = await request.form()
            message = str(form.get("message") or "")
            session_id = str(form.get("session_id") or "")
            choices_json = str(form.get("choices_json") or "")
            f = form.get("file")
            if hasattr(f, "read"):
                raw = await f.read()
                target_gb = raw.decode("utf-8", errors="ignore")
            invs = form.getlist("inventory_files") if hasattr(form, "getlist") else []
            for inv in invs:
                if hasattr(inv, "read"):
                    raw = await inv.read()
                    inventory_gbs.append(raw.decode("utf-8", errors="ignore"))
        elif "application/json" in content_type:
            body = await request.json()
            if not isinstance(body, dict):
                return JSONResponse(
                    {"ok": False, "reply": "JSON body must be an object",
                     "error": "bad_request", "files": None, "viz": None,
                     "agent_trace": [], "n_tool_calls": 0},
                    status_code=400,
                )
            message = str(body.get("message") or "")
            session_id = str(body.get("session_id") or "")
            choices = body.get("choices") if isinstance(body.get("choices"), list) else None
            if choices:
                choices_json = json.dumps(choices)
            # Allow inline target / inventory as GenBank-text strings
            tgb = body.get("target_genbank")
            if isinstance(tgb, str) and tgb.strip():
                target_gb = tgb
            invs = body.get("inventory_genbank") or []
            if isinstance(invs, list):
                for v in invs:
                    if isinstance(v, str) and v.strip():
                        inventory_gbs.append(v)
        else:
            # Fallback: try to parse as form first, then JSON
            try:
                form = await request.form()
                message = str(form.get("message") or "")
                session_id = str(form.get("session_id") or "")
            except Exception:
                try:
                    body = await request.json()
                    message = str((body or {}).get("message") or "")
                except Exception:
                    return JSONResponse(
                        {"ok": False, "reply": (
                            f"Unsupported Content-Type: {content_type!r}. "
                            "Send multipart/form-data or application/json."),
                         "error": "unsupported_content_type", "files": None,
                         "viz": None, "agent_trace": [], "n_tool_calls": 0},
                        status_code=415,
                    )
    except Exception as e:
        logger.exception("agent/chat input parse failed")
        return JSONResponse(
            {"ok": False, "reply": f"Bad request: {e}",
             "error": "bad_request", "files": None, "viz": None,
             "agent_trace": [], "n_tool_calls": 0},
            status_code=400,
        )

    choices = None
    if choices_json:
        try:
            parsed = json.loads(choices_json)
            if isinstance(parsed, list):
                choices = parsed
        except Exception as e:
            logger.warning("agent/chat: bad choices_json: %s", e)

    logger.info(
        "agent/chat: ct=%s msg_len=%d target=%s inv_count=%d choices=%d",
        content_type.split(";")[0] or "?",
        len(message or ""), bool(target_gb), len(inventory_gbs),
        len(choices) if choices else 0,
    )

    result = await run_agent(
        user_message=message,
        target_gb=target_gb,
        inventory_gbs=inventory_gbs,
        choices=choices,
    )
    payload = await _build_response_payload(result)
    return JSONResponse(payload, status_code=200)
