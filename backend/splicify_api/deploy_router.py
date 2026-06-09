"""
Webhook-driven deploy router.

POST /api/admin/deploy
    Body: GitHub push-event JSON.
    Auth: X-Hub-Signature-256 HMAC of the raw body against a shared
          secret kept at /etc/aiplasmiddesign/webhook-secret (or the
          AIPLASMIDDESIGN_WEBHOOK_SECRET env var; env wins if set).

On a verified push to refs/heads/main, the handler queues an out-of-band
deploy via systemd-run:

    systemd-run --unit=aiplasmiddesign-deploy --on-active=1 \\
        /usr/local/bin/aiplasmiddesign-deploy

The deploy script (in this repo at deploy/aiplasmiddesign-deploy.sh)
stashes any local changes, fast-forwards main, and restarts the API
service. systemd-run lets the HTTP response go out before the API
process is restarted by itself.
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

SECRET_FILE = Path("/etc/aiplasmiddesign/webhook-secret")
DEPLOY_SCRIPT = "/usr/local/bin/aiplasmiddesign-deploy"
DEPLOY_UNIT = "aiplasmiddesign-deploy"


def _load_secret() -> bytes | None:
    """Resolve the webhook secret.

    Env var wins over the on-disk file so secret rotation can happen
    without touching the filesystem (set in the systemd unit's
    EnvironmentFile). Returns None when neither is set — the endpoint
    fails closed in that case.
    """
    env_secret = os.environ.get("AIPLASMIDDESIGN_WEBHOOK_SECRET")
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
    # Constant-time compare to avoid leaking the prefix on timing channels.
    return hmac.compare_digest(expected, header)


@router.post("/api/admin/deploy")
async def github_webhook_deploy(request: Request) -> Response:
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
        # Non-main pushes are no-ops; we only deploy main.
        return {"queued": False, "ref": ref, "reason": "non-main ref"}

    head_sha = (payload.get("after") or "")[:12]

    # Detach the deploy so the API can finish this request before it
    # restarts itself. --on-active=1 schedules a one-shot timer that
    # fires ~1 s later; ExecStart runs as root (same user as the API).
    try:
        subprocess.run(
            [
                "systemd-run",
                f"--unit={DEPLOY_UNIT}",
                "--on-active=1",
                "--description=AI Plasmid Design webhook-triggered deploy",
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
