"""agent_v2_api — FastAPI entry."""
from fastapi import FastAPI

from agent_v2.deploy_router import router as agent_v2_deploy_router
from agent_v2.router import router as agent_v2_router

app = FastAPI(title="agent_v2_api", version="0.1.0")
app.include_router(agent_v2_router, prefix="/agent_v2")
# Webhook-driven auto-deploy. Mirrors v1 aiplasmiddesign_api's
# /api/admin/deploy contract: POST /agent_v2/admin/deploy with a
# GitHub push payload HMAC-signed against the shared secret at
# /etc/agent_v2/webhook-secret (or AGENT_V2_WEBHOOK_SECRET env var).
app.include_router(agent_v2_deploy_router, prefix="/agent_v2")
