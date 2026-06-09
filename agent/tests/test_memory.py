"""Round-trip and TTL tests for agent_v2.memory."""
import os
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
os.environ.setdefault("REDIS_URL", "redis://localhost:6379/1")

import pytest

from agent_v2 import memory


@pytest.fixture
def fresh_session():
    sid = memory.mint_session_id()
    yield sid
    memory.reset(sid)


def test_mint_session_id_format():
    sid = memory.mint_session_id()
    assert sid.startswith("sess_")
    assert len(sid) == len("sess_") + 12


def test_load_missing_returns_none():
    assert memory.load("sess_nonexistent_zzzz") is None


def test_save_and_load_roundtrip(fresh_session):
    state = memory.SessionState(
        session_id=fresh_session,
        messages=[{"role": "user", "content": "Hello"}],
        registry_summary=[{"attachment_id": "att_1", "name": "test"}],
        decisions=[{"choice": "Gibson", "reason": "PCR-friendly"}],
        last_user_message="Hello",
    )
    memory.save(state)
    loaded = memory.load(fresh_session)
    assert loaded is not None
    assert loaded.session_id == fresh_session
    assert loaded.messages == [{"role": "user", "content": "Hello"}]
    assert loaded.registry_summary[0]["attachment_id"] == "att_1"
    assert loaded.decisions[0]["choice"] == "Gibson"
    assert loaded.updated_at > 0


def test_reset_removes_state(fresh_session):
    memory.save(memory.SessionState(session_id=fresh_session, last_user_message="x"))
    assert memory.load(fresh_session) is not None
    memory.reset(fresh_session)
    assert memory.load(fresh_session) is None


def test_is_new_topic_no_session():
    assert memory.is_new_topic(None, "anything") is True
    assert memory.is_new_topic("", "anything") is True


def test_is_new_topic_existing_session(fresh_session):
    memory.save(memory.SessionState(
        session_id=fresh_session,
        messages=[{"role": "user", "content": "Hi"}],
    ))
    assert memory.is_new_topic(fresh_session, "follow-up") is False


def test_ttl_set(fresh_session):
    import redis as redis_mod

    memory.save(memory.SessionState(session_id=fresh_session, last_user_message="x"))
    client = redis_mod.from_url(os.environ["REDIS_URL"], decode_responses=True)
    ttl = client.ttl(memory._key(fresh_session))
    assert 29 * 24 * 60 * 60 < ttl <= 30 * 24 * 60 * 60
