import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from fastapi.testclient import TestClient
from main import app


def test_health():
    client = TestClient(app)
    r = client.get("/agent_v2/health")
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True
    assert body["service"] == "agent_v2_api"
