"""Webhook-driven deploy router for agent_v2_api.

Mirrors `aiplasmiddesign_api/backend/splicify_api/deploy_router.py` —
same HMAC-SHA256 GitHub-webhook contract, same out-of-band `systemd-run`
dispatch pattern. The only differences are paths + secret location +
the systemd unit name so this can sit alongside the v1 deploy without
collisions.

POST /agent_v2/admin/deploy
    Body: GitHub push-event JSON.
    Auth: X-Hub-Signature-256 HMAC of the raw body against the shared
          secret at /etc/agent_v2/webhook-secret (or the
          AGENT_V2_WEBHOOK_SECRET env var; env wins).

On a verified push to refs/heads/main the handler queues an
out-of-band deploy via systemd-run:

    systemd-run --unit=agent_v2-deploy --on-active=1 \
        /usr/local/bin/agent_v2-deploy

The deploy script (in this repo at deploy/agent_v2-deploy.sh) stashes
any local edits, fast-forwards main, and restarts the agent_v2_api
service. systemd-run detaches the work so the HTTP response goes out
before the API restarts itself.
"""
from __future__ import annotations

import hashlib
import hmac
import json
import logging
import os
import subprocess
from pathlib import Path

from fastapi import APIRouter, Request, Response, status


router = APIRouter()
logger = logging.getLogger(__name__)

SECRET_FILE = Path("/etc/agent_v2/webhook-secret")
DEPLOY_SCRIPT = "/usr/local/bin/agent_v2-deploy"
DEPLOY_UNIT = "agent_v2-deploy"


def _load_secret() -> bytes | None:
    """Resolve the webhook secret. Env wins over file so rotations can
    happen via systemd EnvironmentFile without touching disk."""
    env_secret = os.environ.get("AGENT_V2_WEBHOOK_SECRET")
    if env_secret:
        return env_secret.encode("utf-8")
    try:
        return SECRET_FILE.read_bytes().strip()
    except FileNotFoundError:
        return None


def _verify_signature(secret: bytes, body: bytes, header: str | None) -> bool:
    if not header or not header.startswith("sha256="):
        return False
    expected = "sha256=" + hmac.new(secret, body, hashlib.sha256).hexdigest()
    return hmac.compare_digest(expected, header)


@router.post("/admin/deploy")
async def github_webhook_deploy(request: Request) -> Response:
    """Receive GitHub push webhooks and queue a deploy on main pushes."""
    secret = _load_secret()
    if secret is None:
        logger.error("deploy webhook called but no secret is configured")
        return Response(status_code=status.HTTP_503_SERVICE_UNAVAILABLE)

    body = await request.body()
    sig = request.headers.get("X-Hub-Signature-256")
    if not _verify_signature(secret, body, sig):
        logger.warning("deploy webhook signature mismatch")
        return Response(status_code=status.HTTP_403_FORBIDDEN)

    event = request.headers.get("X-GitHub-Event", "")
    if event == "ping":
        return {"pong": True}
    if event != "push":
        return Response(status_code=status.HTTP_204_NO_CONTENT)

    try:
        payload = json.loads(body or b"{}")
    except json.JSONDecodeError:
        return Response(status_code=status.HTTP_400_BAD_REQUEST)

    ref = payload.get("ref")
    if ref != "refs/heads/main":
        # Only main pushes deploy; everything else is a no-op acknowledgement.
        return {"queued": False, "ref": ref, "reason": "non-main ref"}

    head_sha = (payload.get("after") or "")[:12]

    try:
        subprocess.run(
            [
                "systemd-run",
                f"--unit={DEPLOY_UNIT}",
                "--on-active=1",
                "--description=agent_v2_api webhook-triggered deploy",
                DEPLOY_SCRIPT,
            ],
            check=True,
            capture_output=True,
            timeout=10,
        )
    except FileNotFoundError:
        logger.error("systemd-run not found; cannot dispatch deploy")
        return Response(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR)
    except subprocess.CalledProcessError as e:
        logger.error("systemd-run failed: %s", e.stderr.decode("utf-8", "replace"))
        return Response(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR)

    return {"queued": True, "ref": "main", "head": head_sha}
