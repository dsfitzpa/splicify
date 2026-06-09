"""Redis-backed conversation memory for agent_v2.

One key per session at `agent_v2:session:<id>`. JSON-serialised SessionState
holds the rolling Anthropic message list, an AttachmentRegistry public summary
(no raw DNA), a decisions ledger, and the last user message. TTL 30 days.

DB 1 by default to isolate from v1's implicit DB 0; override with REDIS_URL.
"""
from __future__ import annotations

import json
import os
import time
import uuid
from dataclasses import asdict, dataclass, field
from typing import Any, Optional

import redis

_KEY_PREFIX = "agent_v2:session:"
_TTL_SECONDS = 30 * 24 * 60 * 60  # 30 days


@dataclass
class SessionState:
    session_id: str
    messages: list[dict[str, Any]] = field(default_factory=list)
    registry_summary: list[dict[str, Any]] = field(default_factory=list)
    decisions: list[dict[str, Any]] = field(default_factory=list)
    last_user_message: str = ""
    updated_at: float = 0.0


def _client() -> "redis.Redis":
    url = os.getenv("REDIS_URL", "redis://localhost:6379/1")
    return redis.from_url(url, decode_responses=True)


def _key(session_id: str) -> str:
    return f"{_KEY_PREFIX}{session_id}"


def mint_session_id() -> str:
    return f"sess_{uuid.uuid4().hex[:12]}"


def load(session_id: str) -> Optional[SessionState]:
    raw = _client().get(_key(session_id))
    if not raw:
        return None
    return SessionState(**json.loads(raw))


def save(state: SessionState) -> None:
    state.updated_at = time.time()
    _client().setex(_key(state.session_id), _TTL_SECONDS, json.dumps(asdict(state)))


def reset(session_id: str) -> None:
    _client().delete(_key(session_id))


def is_new_topic(session_id: Optional[str], user_message: str) -> bool:
    """Heuristic: new topic if no session id, or stored session is empty.

    The LLM-based new-topic agent (subagents/new_topic.py) will refine this
    later; this baseline lets the orchestrator route correctly today.
    """
    if not session_id:
        return True
    state = load(session_id)
    if state is None or not state.messages:
        return True
    return False
