"""Tests for POST /agent_v2/chat-stream — orchestrator monkeypatched."""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest
from fastapi.testclient import TestClient

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.triage import TriageResult
from main import app


def _parse_sse(body: str) -> list[dict]:
    """Tiny SSE parser. Returns a list of {event, data} dicts.

    Skips comment-only chunks (lines starting with :) such as the
    initial :start padding the endpoint emits to flush CDN buffers.
    """
    out = []
    for chunk in body.strip().split("\n\n"):
        stripped = chunk.strip()
        if not stripped:
            continue
        non_comment_lines = [ln for ln in chunk.splitlines() if not ln.startswith(":") and ln.strip()]
        if not non_comment_lines:
            continue
        evt = {"event": "message", "data": ""}
        for line in chunk.splitlines():
            if line.startswith("event: "):
                evt["event"] = line[len("event: "):]
            elif line.startswith("data: "):
                evt["data"] = line[len("data: "):]
        out.append(evt)
    return out


@pytest.fixture
def patched_run_orchestrator(monkeypatch):
    """Replace the orchestrator with a controllable fake.

    The fake calls the supplied `on_triage` callback then returns a fixed
    envelope, so we can assert on the streamed event order + payloads.
    """
    captured = {}

    async def fake(message, registry, *, session_id=None, deps=None, on_triage=None, on_tool_event=None):
        captured["message"] = message
        captured["registry_size"] = len(registry.public_summary())
        if on_triage is not None:
            await on_triage(TriageResult(
                intent="PLASMID_CLONING",
                shorthand="Gibson into pUC19",
                is_new_topic=True,
                rejection_reason=None,
            ))
        return {
            "ok": True,
            "reply": "Done.",
            "files": None,
            "viz": None,
            "session_id": session_id or "sess_test",
            "intent": "PLASMID_CLONING",
            "shorthand": "Gibson into pUC19",
            "is_new_topic": True,
            "agent_trace": [],
            "n_tool_calls": 0,
            "error": None,
        }

    import agent_v2.router as r
    monkeypatch.setattr(r, "run_orchestrator", fake)
    return captured


def test_chat_stream_emits_shorthand_then_envelope(patched_run_orchestrator):
    client = TestClient(app)
    r = client.post("/agent_v2/chat-stream", json={"message": "Build a Gibson"})
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/event-stream")

    events = _parse_sse(r.text)
    assert len(events) == 2
    assert events[0]["event"] == "shorthand"
    assert events[1]["event"] == "envelope"

    import json as _json
    sh = _json.loads(events[0]["data"])
    assert sh["shorthand"] == "Gibson into pUC19"
    assert sh["intent"] == "PLASMID_CLONING"
    assert sh["is_new_topic"] is True

    env = _json.loads(events[1]["data"])
    assert env["ok"] is True
    assert env["reply"] == "Done."


def test_chat_stream_passes_target_genbank_into_registry(patched_run_orchestrator):
    gb = ("LOCUS test 800 bp DNA circular\nORIGIN\n"
          + "        1 " + "acgtacgtac " * 8 + "\n//\n")
    client = TestClient(app)
    r = client.post(
        "/agent_v2/chat-stream",
        json={"message": "annotate", "target_genbank": gb},
    )
    assert r.status_code == 200
    # The fake captured the registry size:
    assert patched_run_orchestrator["registry_size"] == 1


def test_chat_stream_inventory_list(patched_run_orchestrator):
    gb1 = "ACGTACGTAC" * 30
    gb2 = "GGCCAATTCC" * 30
    client = TestClient(app)
    r = client.post(
        "/agent_v2/chat-stream",
        json={"message": "...", "inventory_genbank": [gb1, gb2]},
    )
    assert r.status_code == 200
    assert patched_run_orchestrator["registry_size"] == 2


def test_chat_stream_rejects_non_json_gracefully(patched_run_orchestrator):
    client = TestClient(app)
    # POST with no body / wrong content-type — orchestrator still gets called
    # with empty message; we just need to verify it doesn't 500.
    r = client.post("/agent_v2/chat-stream",
                    content=b"not json", headers={"content-type": "text/plain"})
    assert r.status_code == 200
    events = _parse_sse(r.text)
    assert any(e["event"] == "envelope" for e in events)


def test_chat_stream_orchestrator_exception_emits_error_event(monkeypatch):
    async def boom(message, registry, *, session_id=None, deps=None, on_triage=None, on_tool_event=None):
        if on_triage is not None:
            await on_triage(TriageResult(intent="REJECT", shorthand="x",
                                          is_new_topic=True, rejection_reason="r"))
        raise RuntimeError("orchestrator crashed mid-run")

    import agent_v2.router as rmod
    monkeypatch.setattr(rmod, "run_orchestrator", boom)

    client = TestClient(app)
    r = client.post("/agent_v2/chat-stream", json={"message": "x"})
    events = _parse_sse(r.text)
    event_names = [e["event"] for e in events]
    assert "shorthand" in event_names
    assert "error" in event_names
    err_evt = next(e for e in events if e["event"] == "error")
    import json as _json
    assert "RuntimeError" in _json.loads(err_evt["data"])["error"]
