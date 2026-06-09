"""emit_assembled_gb — first of four output emitters.

Wraps v1's `_gb_for_attachment` GenBank serialiser. Returns a chat-envelope
file dict (`{fileName, dataBase64}`) and optionally writes the GenBank text
to `<output_dir>/assembled.gb` for the workflow trace to reference.
"""
from __future__ import annotations

import base64
import pathlib
from typing import Any, Optional

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.outputs import prefixed_filename, derive_descriptor


async def emit_assembled_gb(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    aid = args.get("attachment_id")
    att = registry.get(aid) if aid else None
    if att is None:
        return {"ok": False, "error": f"unknown attachment_id: {aid!r}"}

    from splicify_api.agent.agent_tools import _gb_for_attachment
    gb_text = _gb_for_attachment(att)

    _descriptor = derive_descriptor(args)
    file_envelope = {
        "fileName": prefixed_filename("assembled.gb", _descriptor),
        "dataBase64": base64.b64encode(gb_text.encode("utf-8")).decode("ascii"),
    }

    written_path: Optional[str] = None
    if output_dir is not None:
        out = pathlib.Path(output_dir) / prefixed_filename("assembled.gb", _descriptor)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(gb_text)
        written_path = str(out)

    return {
        "ok": True,
        "file": file_envelope,
        "attachment_id": aid,
        "length_bp": len(att.sequence),
        "topology": "circular" if att.circular else "linear",
        "written_path": written_path,
    }
