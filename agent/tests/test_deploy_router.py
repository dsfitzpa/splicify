"""Tests for the agent_v2_api auto-deploy webhook.

Covers HMAC-signature contract, GitHub event routing (ping / push /
others), main-branch gating, and the systemd-run dispatch path. Real
systemd-run is never invoked — it's patched out via monkeypatch so the
test suite can run anywhere.
"""
import hashlib
import hmac
import json
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest
from fastapi.testclient import TestClient

import agent_v2  # noqa: F401


SECRET = b"test-secret-for-the-deploy-webhook"


def _sig(body: bytes) -> str:
    return "sha256=" + hmac.new(SECRET, body, hashlib.sha256).hexdigest()


@pytest.fixture
def client(monkeypatch, tmp_path):
    """TestClient with secret resolution monkeypatched + systemd-run mocked."""
    monkeypatch.setenv("AGENT_V2_WEBHOOK_SECRET", SECRET.decode())
    from main import app
    return TestClient(app)


def _capture_subprocess_run(monkeypatch):
    """Patch subprocess.run inside deploy_router to capture the dispatch
    args without actually invoking systemd-run."""
    import agent_v2.deploy_router as dr
    calls = []

    class _DummyCompleted:
        returncode = 0
        stdout = b""
        stderr = b""

    def fake_run(*args, **kwargs):
        calls.append({"args": args, "kwargs": kwargs})
        return _DummyCompleted()

    monkeypatch.setattr(dr.subprocess, "run", fake_run)
    return calls


# ---------------------------------------------------------------------------
# HMAC signature contract
# ---------------------------------------------------------------------------
def test_missing_signature_returns_403(client):
    body = json.dumps({"ref": "refs/heads/main"}).encode()
    r = client.post(
        "/agent_v2/admin/deploy",
        content=body,
        headers={"X-GitHub-Event": "push", "Content-Type": "application/json"},
    )
    assert r.status_code == 403


def test_bad_signature_returns_403(client):
    body = json.dumps({"ref": "refs/heads/main"}).encode()
    r = client.post(
        "/agent_v2/admin/deploy",
        content=body,
        headers={
            "X-GitHub-Event": "push",
            "X-Hub-Signature-256": "sha256=" + "0" * 64,
            "Content-Type": "application/json",
        },
    )
    assert r.status_code == 403


def test_no_secret_configured_returns_503(monkeypatch, tmp_path):
    """When the env secret is unset AND /etc/agent_v2/webhook-secret is
    missing, the endpoint fails closed at 503."""
    monkeypatch.delenv("AGENT_V2_WEBHOOK_SECRET", raising=False)
    import agent_v2.deploy_router as dr
    monkeypatch.setattr(dr, "SECRET_FILE", tmp_path / "does-not-exist")
    from main import app
    c = TestClient(app)
    body = b"{}"
    r = c.post(
        "/agent_v2/admin/deploy",
        content=body,
        headers={"X-GitHub-Event": "push", "X-Hub-Signature-256": _sig(body)},
    )
    assert r.status_code == 503


# ---------------------------------------------------------------------------
# GitHub event routing
# ---------------------------------------------------------------------------
def test_ping_event_returns_pong(client, monkeypatch):
    _capture_subprocess_run(monkeypatch)
    body = b"{}"
    r = client.post(
        "/agent_v2/admin/deploy",
        content=body,
        headers={"X-GitHub-Event": "ping", "X-Hub-Signature-256": _sig(body)},
    )
    assert r.status_code == 200
    assert r.json() == {"pong": True}


def test_non_push_non_ping_event_returns_204(client, monkeypatch):
    _capture_subprocess_run(monkeypatch)
    body = b"{}"
    r = client.post(
        "/agent_v2/admin/deploy",
        content=body,
        headers={"X-GitHub-Event": "pull_request", "X-Hub-Signature-256": _sig(body)},
    )
    assert r.status_code == 204


def test_non_main_branch_push_does_not_dispatch(client, monkeypatch):
    calls = _capture_subprocess_run(monkeypatch)
    body = json.dumps({
        "ref": "refs/heads/feat/something", "after": "deadbeef" * 5,
    }).encode()
    r = client.post(
        "/agent_v2/admin/deploy",
        content=body,
        headers={"X-GitHub-Event": "push", "X-Hub-Signature-256": _sig(body)},
    )
    assert r.status_code == 200
    out = r.json()
    assert out["queued"] is False
    assert out["ref"] == "refs/heads/feat/something"
    assert calls == [], "systemd-run must NOT fire on non-main pushes"


def test_main_push_dispatches_via_systemd_run(client, monkeypatch):
    calls = _capture_subprocess_run(monkeypatch)
    body = json.dumps({
        "ref": "refs/heads/main",
        "after": "abc123def456" + "00" * 14,
    }).encode()
    r = client.post(
        "/agent_v2/admin/deploy",
        content=body,
        headers={"X-GitHub-Event": "push", "X-Hub-Signature-256": _sig(body)},
    )
    assert r.status_code == 200
    out = r.json()
    assert out == {"queued": True, "ref": "main", "head": "abc123def456"}
    # systemd-run invoked exactly once with the expected unit + script.
    assert len(calls) == 1
    args = calls[0]["args"][0]
    assert args[0] == "systemd-run"
    assert "--unit=agent_v2-deploy" in args
    assert "--on-active=1" in args
    assert args[-1] == "/usr/local/bin/agent_v2-deploy"


def test_invalid_json_body_returns_400(client, monkeypatch):
    _capture_subprocess_run(monkeypatch)
    body = b"this is not json"
    r = client.post(
        "/agent_v2/admin/deploy",
        content=body,
        headers={"X-GitHub-Event": "push", "X-Hub-Signature-256": _sig(body)},
    )
    assert r.status_code == 400


# ---------------------------------------------------------------------------
# Secret loading precedence
# ---------------------------------------------------------------------------
def test_env_secret_wins_over_file(monkeypatch, tmp_path):
    """If both AGENT_V2_WEBHOOK_SECRET and the file are set, env wins
    so secret rotation can happen without touching disk."""
    file_secret = b"FILE_SECRET"
    env_secret = b"ENV_SECRET"
    secret_file = tmp_path / "webhook-secret"
    secret_file.write_bytes(file_secret)

    monkeypatch.setenv("AGENT_V2_WEBHOOK_SECRET", env_secret.decode())
    import agent_v2.deploy_router as dr
    monkeypatch.setattr(dr, "SECRET_FILE", secret_file)
    assert dr._load_secret() == env_secret


def test_file_secret_used_when_env_unset(monkeypatch, tmp_path):
    monkeypatch.delenv("AGENT_V2_WEBHOOK_SECRET", raising=False)
    file_secret = b"FILE_SECRET\n   "  # trailing whitespace stripped
    secret_file = tmp_path / "webhook-secret"
    secret_file.write_bytes(file_secret)
    import agent_v2.deploy_router as dr
    monkeypatch.setattr(dr, "SECRET_FILE", secret_file)
    assert dr._load_secret() == b"FILE_SECRET"
