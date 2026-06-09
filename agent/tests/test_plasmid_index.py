"""Unit tests for PlasmidIndex — fed a synthetic annotation envelope
matching the shape produced by /plannotate/annotate_sequence_llm."""
from __future__ import annotations

import agent_v2  # noqa: F401 — path shim
from agent_v2.interpreter.plasmid_index import PlasmidIndex


def _make_envelope():
    """A 100-bp 'plasmid' with a single ORF that has a 10-aa N-terminal
    tag followed by a 30-aa Cas9-like region, plus a Pol II expression
    cassette interaction and a guide_expression_cassette module."""
    sequence = (
        "ATGGACTACAAAGATGATGACGATAAA"          # 27 nt: MDYKDDDDK + M (FLAG-ish)
        "GTGAAGGAGCTGCTGAAGGAGCTGCTG"          # 27 nt
        "AAGGAGCTGCTGAAGGAGCTGCTGAAG"          # 27 nt
        "GAGCTGCTGAAGGAGCTGCTGAAGTAA"          # 27 nt + stop
        "ATGCGTCTC"                            # 9 nt — start + BsmBI flank
    )
    return {
        "sequence": sequence,
        "annotations": [
            {"name": "FLAG", "type": "CDS", "start": 0, "end": 24, "direction": 1,
              "kb_data": {"gene_name": "FLAG", "feature_name": "FLAG tag"}},
            {"name": "Cas9", "type": "CDS", "start": 24, "end": 105, "direction": 1,
              "kb_data": {"gene_name": "Cas9", "protein_name": "CRISPR-associated endonuclease 9"}},
            {"name": "U6 promoter", "type": "promoter", "start": 0, "end": 10, "direction": 1,
              "kb_data": {"feature_class": "promoter"}},
            {"name": "gRNA scaffold", "type": "misc_RNA", "start": 90, "end": 100, "direction": 1,
              "kb_data": {}},
        ],
        "modules": [
            {
                "module_type": "guide_expression_cassette",
                "rule_id": "POL3-GG-01",
                "name": "Pol III guide cassette (U6)",
                "start": 0, "end": 100, "strand": 1,
                "submodules": [
                    {"module_type": "pol3_promoter", "start": 0, "end": 10, "name": "U6"},
                    {"module_type": "sgrna_scaffold", "start": 90, "end": 100, "name": "gRNA scaffold"},
                ],
                "golden_gate": {"enzyme": "BsmBI", "stuffer_name": "Stuffer-001"},
                "notes": "BsmBI Golden Gate",
            },
            {
                "module_type": "mammalian_pol2_expression_cassette",
                "name": "Pol II Cas9 cassette",
                "start": 0, "end": 110, "strand": 1,
                "submodules": [
                    {"module_type": "upstream_regulatory_module", "name": "EF-1α", "start": 0, "end": 10,
                      "metadata": {"promoter_name": "EF-1α"}},
                    {"module_type": "downstream_regulatory_module", "name": "bGH polyA", "start": 100, "end": 110,
                      "metadata": {"polya_signal": "bGH polyA"}},
                ],
            },
        ],
        "interactions": [],
        "hierarchical_annotations": [
            {
                "name": "Translation (40 aa)",
                "module_type": "translation",
                "layer": "translation",
                "start": 0, "end": 108, "direction": 1,
                "metadata": {
                    "aa_length": 35,
                    "aa_sequence": "MDYKDDDDKVKELLKELLKELLKELLKELLKELLKE",
                    "feature_regions": [
                        {"name": "FLAG", "aa_start": 1, "aa_end": 8, "feature_type": "CDS"},
                        {"name": "Cas9", "aa_start": 9, "aa_end": 35, "feature_type": "CDS"},
                    ],
                },
            },
        ],
        "cloning_features": {
            "features": [
                {"name": "BsmBI", "feature_family": "restriction_site_IIs",
                  "subtype": "BsmBI", "start": 100, "end": 106},
            ],
        },
    }


def test_summary_lists_module_types_and_features():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope(), name="test.gb")
    s = idx.summary()
    assert s["n_modules"] == 2
    assert "guide_expression_cassette" in s["module_types"]
    assert "Cas9" in s["feature_names"]
    assert s["length_bp"] == 117


def test_find_modules_by_type_substring():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    hits = idx.find_modules("guide_expression")
    assert len(hits) == 1
    assert hits[0]["rule_id"] == "POL3-GG-01"
    assert hits[0]["golden_gate"]["enzyme"] == "BsmBI"


def test_find_modules_by_rule_id():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    hits = idx.find_modules("POL3-GG")
    assert len(hits) == 1


def test_find_features_by_gene_name():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    hits = idx.find_features("Cas9")
    assert len(hits) == 1
    assert hits[0]["kb_data"]["gene_name"] == "Cas9"


def test_find_features_by_kb_protein_name():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    hits = idx.find_features("crispr-associated")
    assert len(hits) == 1


def test_lookup_amino_acid_within_cas9():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    # AA position 2 within Cas9 — feature_regions says Cas9 starts at aa 9
    # so this is the 10th aa of the ORF. The 10th codon of the fixture is GTG (Valine).
    res = idx.lookup_amino_acid("Cas9", 2)
    assert res is not None and res["ok"]
    assert res["letter"] == "V"
    assert res["aa_position_in_orf"] == 10
    assert res["aa_position_in_feature"] == 2
    assert res["amino_acid"] == "Valine"
    assert res["codon"] == "GTG"


def test_lookup_amino_acid_in_flag_tag():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    res = idx.lookup_amino_acid("FLAG", 1)
    assert res["letter"] == "M"
    assert res["aa_position_in_orf"] == 1


def test_lookup_amino_acid_out_of_range():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    res = idx.lookup_amino_acid("Cas9", 500)
    assert res is not None and res["ok"] is False
    assert "past the end" in res["reason"]


def test_expression_cassette_for_falls_back_to_module():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    res = idx.expression_cassette_for("Cas9")
    assert res is not None
    assert res["promoter"]["promoter_name"] == "EF-1α"
    assert res["polyA"]["polya_signal"] == "bGH polyA"


def test_infer_application_picks_crispr_mammalian():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    res = idx.infer_application()
    assert "CRISPR" in res["application"] or "guide" in res["application"].lower()
    assert res["confidence"] in {"high", "medium"}


def test_find_cloning_features_matches_bsmbi():
    idx = PlasmidIndex.from_envelope("p1", _make_envelope())
    hits = idx.find_cloning_features("BsmBI")
    assert len(hits) == 1
    assert hits[0]["start"] == 100
