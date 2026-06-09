"""Tests for emit_guides_csv + emit_guides_gb."""
import asyncio
import base64
import io
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import agent_v2  # noqa: F401
from agent_v2 import tools as v2_tools, attachment_kinds
from agent_v2.file_kind import classify_genbank
from agent_v2.outputs.guides_csv import emit_guides_csv
from agent_v2.outputs.guides_gb import emit_guides_gb
from splicify_api.agent.agent_tools import (
    AttachmentRegistry, extract_seq_from_genbank,
)


FIXTURES = pathlib.Path(__file__).parent / "fixtures"


def _sample_guide(name="g1", start=100, end=120, score=78.0):
    return {
        "name": name, "spacer": "ATGCATGCATGCATGCATGC", "pam": "TGG",
        "start": start, "end": end, "direction": 1,
        "score": score, "score_method": "doench2014",
        "gc_fraction": 0.5, "n_offtargets": 0,
        "target_attachment_id": "att_product_1",
    }


def _sample_pegrna(name="p1", with_ngrna=True):
    p = {
        "name": name, "rank": 1, "predicted_efficiency": 0.42,
        "spacer": "GCGATGCATGCATGCATGCA", "pam": "AGG",
        "spacer_start": 1000, "spacer_end": 1020, "direction": 1,
        "cas9_score": 65.0,
        "rtt": "ATGCATGCATGC", "rtt_length": 12,
        "pbs": "ATGCATGCATGCA", "pbs_length": 13,
        "scaffold": "GTTTT...",
        "full_pegrna": "ACGT" * 30, "full_pegrna_length": 120,
        "is_dpam": False, "is_pe3b": True,
        "edit_type": "substitution", "edit_ref": "A", "edit_alt": "G",
        "target_attachment_id": "att_product_1",
    }
    if with_ngrna:
        p["ngrna"] = {
            "spacer": "TTACGTACGTACGTACGTAC", "pam": "GGG",
            "start": 1100, "end": 1120, "strand": "-",
            "cas9_score": 55.0, "nick_to_pegRNA": 80,
            "is_pe3b": True,
        }
    return p


def _sample_primer_pair(application="sanger"):
    return {
        "application": application,
        "left_primer": "ACGTACGTACGTACGTACGT",
        "left_annealing": "ACGTACGTACGTACGTACGT",
        "left_adapter": "",
        "left_tm": 60.0,
        "right_primer": "TGCATGCATGCATGCATGCA",
        "right_annealing": "TGCATGCATGCATGCATGCA",
        "right_adapter": "",
        "right_tm": 59.5,
        "product_size": 280,
        "pair_label": "KEAP1_amp",
        "name_fwd": "KEAP1_amp_SF",
        "name_rev": "KEAP1_amp_SR",
        "left_pos_plasmid": 850,
        "right_pos_plasmid": 1130,
        "region_start_plasmid": 600,
        "region_end_plasmid": 1400,
        "target_attachment_id": "att_product_1",
    }


# ---------------------------------------------------------------------------
# emit_guides_csv
# ---------------------------------------------------------------------------
def test_guides_csv_writes_one_row_per_entry():
    reg = AttachmentRegistry()
    reg.register_product("test", "ACGT" * 500, circular=True)
    result = asyncio.run(emit_guides_csv(
        {
            "guides": [_sample_guide(), _sample_guide(name="g2", start=200, end=220)],
            "pegrnas": [_sample_pegrna()],            # produces 1 pegRNA + 1 ngRNA row
            "primers": [_sample_primer_pair()],       # produces 2 rows (fwd + rev)
            "cloning_oligos": [{"name": "g1_oligo", "spacer": "ATGC",
                                  "oligo_top": "CACCGATGC", "oligo_bottom": "AAACGCATC"}],
        }, reg,
    ))
    assert result["ok"] is True
    assert result["file"]["fileName"] == "guides.csv"
    # 2 sgRNA + 1 pegRNA + 1 ngRNA + 2 primers + 1 oligo = 7 rows
    assert result["n_rows"] == 7
    assert result["n_sgRNAs"] == 2
    assert result["n_pegRNAs"] == 1
    assert result["n_primer_pairs"] == 1
    assert result["n_cloning_oligos"] == 1

    csv_text = base64.b64decode(result["file"]["dataBase64"]).decode("utf-8")
    lines = csv_text.strip().splitlines()
    assert lines[0].split(",")[0] == "type"
    assert sum(1 for l in lines if l.startswith("sgRNA,")) == 2
    assert sum(1 for l in lines if l.startswith("pegRNA,")) == 1
    assert sum(1 for l in lines if l.startswith("ngRNA,")) == 1
    assert sum(1 for l in lines if l.startswith("primer_sanger_fwd,")) == 1
    assert sum(1 for l in lines if l.startswith("primer_sanger_rev,")) == 1
    assert sum(1 for l in lines if l.startswith("cloning_oligo,")) == 1


def test_guides_csv_writes_to_disk(tmp_path):
    reg = AttachmentRegistry()
    reg.register_product("test", "ACGT" * 200, circular=True)
    result = asyncio.run(emit_guides_csv(
        {"guides": [_sample_guide()]}, reg, output_dir=str(tmp_path),
    ))
    assert result["written_path"] == str(tmp_path / "guides.csv")
    body = (tmp_path / "guides.csv").read_text()
    assert "sgRNA" in body
    assert "ATGCATGCATGCATGCATGC" in body


def test_guides_csv_empty_input_writes_header_only():
    reg = AttachmentRegistry()
    result = asyncio.run(emit_guides_csv({}, reg))
    assert result["ok"] is True
    assert result["n_rows"] == 0
    csv_text = base64.b64decode(result["file"]["dataBase64"]).decode("utf-8")
    assert csv_text.splitlines()[0].startswith("type,")


# ---------------------------------------------------------------------------
# emit_guides_gb
# ---------------------------------------------------------------------------
def test_guides_gb_unknown_attachment():
    reg = AttachmentRegistry()
    result = asyncio.run(emit_guides_gb({"target_attachment_id": "att_nope"}, reg))
    assert result["ok"] is False
    assert "unknown" in result["error"].lower()


def test_guides_gb_minimal_record_fallback_when_no_cached_gb():
    """No gb_text stashed -> the emitter builds a minimal record from the
    registry sequence + appends features."""
    attachment_kinds.clear_all_kinds()
    reg = AttachmentRegistry()
    aid = reg.register_product("test_plasmid", "ACGT" * 500, circular=True)
    result = asyncio.run(emit_guides_gb(
        {"target_attachment_id": aid,
         "guides": [_sample_guide(start=100, end=120)],
         "primers": [_sample_primer_pair()]},
        reg,
    ))
    assert result["ok"] is True
    assert result["n_features_added"] == 1 + 2  # 1 sgRNA + 2 primer_bind
    gb_text = base64.b64decode(result["file"]["dataBase64"]).decode("utf-8")
    assert gb_text.startswith("LOCUS")
    assert "misc_RNA" in gb_text
    assert "primer_bind" in gb_text
    # sgRNA name appears as /label.
    assert "g1" in gb_text


def test_guides_gb_uses_cached_gb_text_when_available():
    """If router cached the original .gb, the emitter reuses it so existing
    annotations (gene/mRNA/CDS) survive."""
    attachment_kinds.clear_all_kinds()
    gb = (FIXTURES / "KEAP1.gb").read_text()
    reg = AttachmentRegistry()
    aid = reg.register_product("KEAP1", extract_seq_from_genbank(gb), circular=False)
    attachment_kinds.stash_kind(aid, classify_genbank(gb), gb_text=gb)

    result = asyncio.run(emit_guides_gb(
        {"target_attachment_id": aid,
         "guides": [_sample_guide(start=2300, end=2320)]},
        reg,
    ))
    assert result["ok"] is True
    assert result["n_features_added"] == 1
    gb_text = base64.b64decode(result["file"]["dataBase64"]).decode("utf-8")
    # KEAP1 annotations survive.
    assert "/gene=\"KEAP1\"" in gb_text
    # New sgRNA feature appended.
    assert "misc_RNA" in gb_text


def test_guides_gb_writes_to_disk(tmp_path):
    attachment_kinds.clear_all_kinds()
    reg = AttachmentRegistry()
    aid = reg.register_product("test", "ACGT" * 500, circular=True)
    result = asyncio.run(emit_guides_gb(
        {"target_attachment_id": aid, "guides": [_sample_guide()]},
        reg, output_dir=str(tmp_path),
    ))
    assert result["written_path"] == str(tmp_path / "guides.gb")
    assert (tmp_path / "guides.gb").read_text().startswith("LOCUS")


def test_guides_gb_reverse_strand_primer_placement():
    attachment_kinds.clear_all_kinds()
    reg = AttachmentRegistry()
    aid = reg.register_product("test", "ACGT" * 500, circular=True)
    pp = _sample_primer_pair()
    pp["right_pos_plasmid"] = 1130
    pp["right_annealing"] = "TGCATGCATGCATGCATGCA"  # 20 nt
    result = asyncio.run(emit_guides_gb(
        {"target_attachment_id": aid, "primers": [pp]}, reg,
    ))
    gb_text = base64.b64decode(result["file"]["dataBase64"]).decode("utf-8")
    # Reverse primer should be on complement strand. SeqIO writes that as
    # `complement(...)` location.
    assert "complement" in gb_text


# ---------------------------------------------------------------------------
# Dispatch chain + roster
# ---------------------------------------------------------------------------
def test_dispatch_routes_emit_guides_csv():
    reg = AttachmentRegistry()
    reg.register_product("test", "ACGT" * 100, circular=True)
    result = asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_guides_csv", {"guides": [_sample_guide()]}, reg,
    ))
    assert result["ok"] is True
    assert result["file"]["fileName"] == "guides.csv"


def test_dispatch_routes_emit_guides_gb():
    attachment_kinds.clear_all_kinds()
    reg = AttachmentRegistry()
    aid = reg.register_product("test", "ACGT" * 100, circular=True)
    result = asyncio.run(v2_tools.dispatch_with_emitters(
        "emit_guides_gb", {"target_attachment_id": aid}, reg,
    ))
    assert result["ok"] is True
    assert result["file"]["fileName"] == "guides.gb"


def test_full_tool_roster_includes_new_emitters():
    roster = v2_tools.make_full_tool_roster()
    names = {t["name"] for t in roster}
    assert "emit_guides_csv" in names
    assert "emit_guides_gb" in names
    # Existing emitters still present.
    assert "emit_assembled_gb" in names
    assert "emit_parts_order" in names
