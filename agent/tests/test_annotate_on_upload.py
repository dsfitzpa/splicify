"""Tests for POST /agent_v2/annotate-on-upload — annotate_llm_cached mocked."""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest
from fastapi.testclient import TestClient

import agent_v2  # noqa: F401 — triggers path shim
from main import app


_GB_TEXT = """LOCUS       test_plasmid          800 bp    DNA     circular SYN 01-JAN-2026
DEFINITION  test
FEATURES             Location/Qualifiers
     source          1..800
                     /organism="synthetic DNA construct"
                     /mol_type="other DNA"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
       61 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
      121 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""


@pytest.fixture
def patched_annotate_llm_cached(monkeypatch):
    async def fake(seq, *, circular=True, depth="full"):
        return {
            "annotations": [{"name": "GFP", "start": 100, "end": 800}],
            "modules": [{"name": "expression cassette"}],
            "cloning_features": [{"name": "EcoRI", "start": 50}],
            "interactions": [],
            "hierarchical_annotations": [{"name": "root"}],
        }
    import splicify_api.annotation_cache as ac
    monkeypatch.setattr(ac, "annotate_llm_cached", fake)
    return fake


def test_annotate_on_upload_happy_path(patched_annotate_llm_cached):
    client = TestClient(app)
    r = client.post(
        "/agent_v2/annotate-on-upload",
        files={"file": ("test_plasmid.gb", _GB_TEXT, "application/octet-stream")},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True
    assert body["title"] == "test_plasmid"
    assert body["length_bp"] > 0
    viz = body["viz"]
    assert viz["type"] == "plasmid"
    assert viz["title"] == "test_plasmid"
    assert viz["circular"] is True
    assert viz["annotations"][0]["name"] == "GFP"
    assert viz["modules"][0]["name"] == "expression cassette"
    assert viz["cloning_features"][0]["name"] == "EcoRI"
    assert viz["hierarchical_annotations"][0]["name"] == "root"


def test_annotate_on_upload_empty_file(patched_annotate_llm_cached):
    client = TestClient(app)
    r = client.post(
        "/agent_v2/annotate-on-upload",
        files={"file": ("empty.gb", "", "application/octet-stream")},
    )
    body = r.json()
    assert body["ok"] is False
    assert "empty" in body["error"].lower()


def test_annotate_on_upload_no_sequence_in_file(patched_annotate_llm_cached):
    """A file with non-DNA content (no ORIGIN block, no DNA letters) should
    surface a clear error rather than running the pipeline on garbage.
    """
    client = TestClient(app)
    r = client.post(
        "/agent_v2/annotate-on-upload",
        files={"file": ("readme.md", "Just some markdown notes!!!\n# Header\n",
                        "text/markdown")},
    )
    body = r.json()
    assert body["ok"] is False
    assert "no plasmid sequence" in body["error"].lower()


def test_annotate_on_upload_annotation_failure(monkeypatch):
    async def boom(seq, *, circular=True, depth="full"):
        raise RuntimeError("annotation pipeline crashed")
    import splicify_api.annotation_cache as ac
    monkeypatch.setattr(ac, "annotate_llm_cached", boom)

    client = TestClient(app)
    r = client.post(
        "/agent_v2/annotate-on-upload",
        files={"file": ("test.gb", _GB_TEXT, "application/octet-stream")},
    )
    body = r.json()
    assert body["ok"] is False
    assert "annotation failed" in body["error"]
    assert "RuntimeError" in body["error"]


def test_annotate_on_upload_strips_extension_from_title(patched_annotate_llm_cached):
    client = TestClient(app)
    r = client.post(
        "/agent_v2/annotate-on-upload",
        files={"file": ("my-plasmid.genbank", _GB_TEXT, "application/octet-stream")},
    )
    body = r.json()
    assert body["title"] == "my-plasmid"


def test_annotate_on_upload_accepts_raw_dna_no_locus(patched_annotate_llm_cached):
    """`extract_seq_from_genbank` falls back to raw DNA letters when no
    ORIGIN block is present. Verify the endpoint accepts that shape too.
    """
    client = TestClient(app)
    r = client.post(
        "/agent_v2/annotate-on-upload",
        files={"file": ("raw.txt", "ACGTACGTACGTACGT" * 50,
                        "application/octet-stream")},
    )
    body = r.json()
    assert body["ok"] is True
    assert body["length_bp"] == 16 * 50
