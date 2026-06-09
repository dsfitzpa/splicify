"""
FastAPI router for the plasmid LM.

Endpoints:
  POST /plasmid_lm/generate   — generate a plasmid design from a description
  GET  /plasmid_lm/health     — is the model loaded + what's its training state

Mounted in main.py. Level-1 scope: token stream + structural parse only; no
sequence materialization, no validator. Materializer/validator slot in later.
"""
from __future__ import annotations
import logging
from typing import Any

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

logger = logging.getLogger("plasmid_lm.router")

router = APIRouter(prefix="/plasmid_lm", tags=["plasmid_lm"])


class GenerateRequest(BaseModel):
    description: str = Field(..., min_length=1, max_length=1000)
    max_new_tokens: int = Field(256, ge=1, le=1024)
    temperature: float = Field(0.8, gt=0.0, le=2.0)
    top_k: int = Field(40, ge=1, le=500)
    seed: int | None = None


@router.get("/health")
def health() -> dict[str, Any]:
    try:
        from .lm_inference import get_lm_service
        svc = get_lm_service()
        return {
            "ok": True,
            "model_loaded": True,
            "training_step": svc.training_step,
            "vocab_size": svc.vocab.size,
            "params_M": round(svc.model.num_params() / 1e6, 2),
            "device": svc.device,
            "note": ("experimental CPU smoke-training checkpoint; "
                     "structural grammar only, description conditioning weak"),
        }
    except Exception as exc:
        logger.exception("LM health check failed")
        return {"ok": False, "model_loaded": False, "error": str(exc)}


@router.post("/generate")
def generate(req: GenerateRequest) -> dict[str, Any]:
    try:
        from .lm_inference import get_lm_service, format_lm_reply
        svc = get_lm_service()
        result = svc.generate(
            description=req.description,
            max_new_tokens=req.max_new_tokens,
            temperature=req.temperature,
            top_k=req.top_k,
            seed=req.seed,
        )
        result["reply_markdown"] = format_lm_reply(result)
        return {"ok": True, **result}
    except Exception as exc:
        logger.exception("LM generate failed")
        raise HTTPException(status_code=500, detail=str(exc))
