"""Smoke-test that v1 splicify_api modules are reachable via the path shim."""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401 — triggers path shim
from splicify_api.agent.agent_tools import (
    AttachmentRegistry,
    _strip_sequences,
)


def test_attachment_registry_creates_and_registers():
    reg = AttachmentRegistry()
    aid = reg.register_product("smoke_test", "ACGTACGTACGT", circular=True)
    assert aid == "att_product_1"
    summary = reg.public_summary()
    assert len(summary) == 1
    assert summary[0]["attachment_id"] == aid
    assert "sequence" not in summary[0]


def test_strip_sequences_redacts_seq_keys():
    payload = {
        "name": "ok",
        "sequence": "ACGTACGT",
        "nested": {"genbank": "LOCUS ..."},
    }
    out = _strip_sequences(payload)
    assert out["name"] == "ok"
    assert out["sequence"] == "[redacted]"
    assert out["nested"]["genbank"] == "[redacted]"


def test_strip_sequences_redacts_long_dna_strings():
    long_dna = "ACGT" * 16  # 64 chars, matches [ACGTNacgtn\\s]+
    out = _strip_sequences({"label": long_dna})
    assert out["label"].startswith("[redacted DNA,")
