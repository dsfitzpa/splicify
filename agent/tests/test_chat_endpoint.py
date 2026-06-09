"""Tests for POST /agent_v2/chat (non-streaming, multipart-or-JSON)."""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest
from fastapi.testclient import TestClient

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2.triage import TriageResult
from main import app


_GB = ("LOCUS test 800 bp DNA circular SYN 01-JAN-2026\nORIGIN\n"
       "        1 " + "acgtacgtac " * 8 + "\n//\n")


@pytest.fixture
def patched_run_orchestrator(monkeypatch):
    """Replace run_orchestrator with a fake that captures inputs and returns a fixed envelope."""
    captured: dict = {}

    async def fake(message, registry, *, session_id=None, deps=None, on_triage=None, on_tool_event=None):
        captured["message"] = message
        captured["session_id"] = session_id
        captured["registry"] = [a["attachment_id"] for a in registry.public_summary()]
        if on_triage is not None:
            await on_triage(TriageResult(intent="PLASMID_CLONING",
                                          shorthand="test", is_new_topic=True,
                                          rejection_reason=None))
        return {
            "ok": True,
            "reply": "echo: " + message,
            "files": None,
            "viz": None,
            "session_id": session_id or "sess_x",
            "intent": "PLASMID_CLONING",
            "shorthand": "test",
            "is_new_topic": True,
            "agent_trace": [],
            "n_tool_calls": 0,
            "error": None,
        }

    import agent_v2.router as r
    monkeypatch.setattr(r, "run_orchestrator", fake)
    return captured


def test_chat_json_path(patched_run_orchestrator):
    client = TestClient(app)
    r = client.post(
        "/agent_v2/chat",
        json={"message": "Build a Gibson", "target_genbank": _GB},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True
    assert body["reply"] == "echo: Build a Gibson"
    assert patched_run_orchestrator["message"] == "Build a Gibson"
    assert patched_run_orchestrator["registry"] == ["att_product_1"]


def test_chat_multipart_with_target_file(patched_run_orchestrator):
    client = TestClient(app)
    r = client.post(
        "/agent_v2/chat",
        data={"message": "annotate this", "session_id": "sess_user_supplied"},
        files={"file": ("my_plasmid.gb", _GB, "application/octet-stream")},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True
    assert patched_run_orchestrator["message"] == "annotate this"
    assert patched_run_orchestrator["session_id"] == "sess_user_supplied"
    assert patched_run_orchestrator["registry"] == ["att_product_1"]


def test_chat_multipart_with_inventory_files(patched_run_orchestrator):
    client = TestClient(app)
    r = client.post(
        "/agent_v2/chat",
        data={"message": "best clone path?"},
        files=[
            ("file", ("target.gb", _GB, "application/octet-stream")),
            ("inventory_files", ("inv1.gb", _GB, "application/octet-stream")),
            ("inventory_files", ("inv2.gb", _GB, "application/octet-stream")),
        ],
    )
    assert r.status_code == 200
    # Target + 2 inventory = 3 registered attachments
    assert patched_run_orchestrator["registry"] == [
        "att_product_1", "att_product_2", "att_product_3",
    ]


def test_chat_orchestrator_exception_returns_error_envelope(monkeypatch):
    async def boom(message, registry, *, session_id=None, deps=None, on_triage=None, on_tool_event=None):
        raise RuntimeError("orchestrator crashed mid-run")
    import agent_v2.router as rmod
    monkeypatch.setattr(rmod, "run_orchestrator", boom)

    client = TestClient(app)
    r = client.post("/agent_v2/chat", json={"message": "x"})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is False
    assert "RuntimeError" in body["error"]
    assert "crashed" in body["error"]


def test_chat_empty_body_runs_with_empty_message(patched_run_orchestrator):
    """No body / wrong content-type: orchestrator still runs with empty message."""
    client = TestClient(app)
    r = client.post("/agent_v2/chat", content=b"", headers={"content-type": "text/plain"})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True
    assert patched_run_orchestrator["message"] == ""
    assert patched_run_orchestrator["registry"] == []


def test_chat_stream_now_supports_multipart(patched_run_orchestrator):
    """Bonus: the /chat-stream endpoint also benefits from the unified parser
    so it can accept multipart bodies just like /chat does."""
    client = TestClient(app)
    r = client.post(
        "/agent_v2/chat-stream",
        data={"message": "build"},
        files={"file": ("plasmid.gb", _GB, "application/octet-stream")},
    )
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/event-stream")
    assert patched_run_orchestrator["registry"] == ["att_product_1"]
