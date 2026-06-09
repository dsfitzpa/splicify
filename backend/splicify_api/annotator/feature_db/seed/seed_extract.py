#!/usr/bin/env python3
"""
seed_extract.py — Extract seed anchor features from SNAPGENE_REBUILD_DB.

Queries the module library for known features and writes them as a FASTA file
that the pipeline uses to bootstrap gold-tier canonical features.
"""
from __future__ import annotations

import json
import logging
import os
import sys
from pathlib import Path
import statistics

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
log = logging.getLogger("seed_extract")

HERE = Path(__file__).parent  # seed/
DB_DIR = HERE.parent          # feature_db/

SEED_FASTA = HERE / "seed_anchor_features.fasta"
SEED_META = HERE / "seed_meta.json"

DATABASE_URL = os.environ.get(
    "DATABASE_URL",
    "postgresql://splicify:0livia-Chizz0@127.0.0.1:5433/SNAPGENE_REBUILD_DB",
)

# Each entry: (pipeline_canonical_id, category, db_search_terms, min_length, display_name, db_canonical_id_hint)
SEED_TARGETS: list[tuple[str, str, list[str], int, str, str]] = [
    # Promoters
    ("PROMOTER_CMV_ENH_PROM", "PROMOTER_POL2", ["PROMOTER_CMV", "ENHANCER_CMV", "CMV"], 100, "CMV enhancer/promoter", "ENHANCER_CMV"),
    ("PROMOTER_EF1A", "PROMOTER_POL2", ["PROMOTER_EF1A", "EF1A"], 100, "EF1a", ""),
    ("PROMOTER_CAG", "PROMOTER_POL2", ["PROMOTER_CAG", "CAG"], 100, "CAG", ""),
    ("PROMOTER_PGK", "PROMOTER_POL2", ["PROMOTER_PGK", "PGK"], 100, "PGK", ""),
    ("PROMOTER_SV40_EARLY", "PROMOTER_POL2", ["PROMOTER_SV40", "SV40"], 50, "SV40 early", ""),
    ("PROMOTER_RSV_LTR", "PROMOTER_POL2", ["PROMOTER_RSV", "RSV"], 100, "RSV LTR", ""),
    ("PROMOTER_CBA", "PROMOTER_POL2", ["PROMOTER_CBA", "CBA"], 100, "CBA", ""),
    ("PROMOTER_TRE3G", "PROMOTER_POL2", ["PROMOTER_TRE", "TRE3G", "TRE"], 100, "TRE3G", ""),
    ("PROMOTER_U6_HUMAN", "PROMOTER_POL3", ["PROMOTER_U6_HUMAN", "PROMOTER_U6", "U6"], 100, "U6 human", "PROMOTER_U6_HUMAN"),
    ("PROMOTER_H1", "PROMOTER_POL3", ["PROMOTER_H1", "H1"], 100, "H1", "PROMOTER_H1"),
    ("PROMOTER_BACT_T7", "PROMOTER_BACTERIAL", ["PROMOTER_T7", "T7"], 15, "T7 promoter", ""),
    ("PROMOTER_BACT_SP6", "PROMOTER_BACTERIAL", ["PROMOTER_SP6", "SP6"], 15, "SP6 promoter", ""),
    ("PROMOTER_BACT_LAC", "PROMOTER_BACTERIAL", ["PROMOTER_LAC", "LAC"], 15, "lac promoter", ""),
    # Resistance markers
    ("CDS_BLA_TEM1", "CDS_RESISTANCE", ["MARKER_AMP", "CDS_AMP", "CDS_BLA", "AMP"], 200, "AmpR (TEM-1)", "MARKER_AMP"),
    ("CDS_KANR_APH3", "CDS_RESISTANCE", ["MARKER_KAN", "CDS_KAN", "KAN"], 200, "KanR", "MARKER_KAN"),
    ("CDS_CMR_CAT", "CDS_RESISTANCE", ["CDS_CMR", "CDS_CAT", "CMR", "CAT"], 200, "CmR", ""),
    ("CDS_PUROR_PAC", "CDS_RESISTANCE", ["CDS_puror", "CDS_PURO", "PURO", "PAC"], 200, "PuroR", "CDS_puror"),
    ("CDS_HYGROR_HPH", "CDS_RESISTANCE", ["CDS_hphmx6", "CDS_HYGRO", "HYGRO", "HPH"], 200, "HygroR", "CDS_hphmx6"),
    ("CDS_NEOR_APH3", "CDS_RESISTANCE", ["CDS_neor_kanr", "CDS_NEO", "NEO", "G418"], 200, "NeoR", "CDS_neor_kanr"),
    ("CDS_BLASTR_BSD", "CDS_RESISTANCE", ["CDS_bsd", "CDS_BLAST", "BSD"], 200, "BlastR", "CDS_bsd"),
    ("CDS_ZEOR_SHBLE", "CDS_RESISTANCE", ["CDS_ZEO", "SH_BLE"], 200, "ZeoR", ""),
    ("CDS_SPCR_AADA", "CDS_RESISTANCE", ["CDS_SPEC", "AADA", "SPEC"], 200, "SpecR", ""),
    # Reporters
    ("CDS_REPORTER_EGFP", "CDS_REPORTER", ["CDS_gfpuv", "CDS_EGFP", "CDS_GFP", "EGFP"], 200, "EGFP", "CDS_gfpuv"),
    ("CDS_REPORTER_MCHERRY", "CDS_REPORTER", ["CDS_mcherry2", "CDS_MCHERRY", "MCHERRY"], 200, "mCherry", "CDS_mcherry2"),
    ("CDS_REPORTER_LUC2", "CDS_REPORTER", ["CDS_LUC", "LUC2", "FLUC"], 200, "Luc2", ""),
    ("CDS_REPORTER_LACZ", "CDS_REPORTER", ["CDS_LACZ", "LACZ"], 200, "LacZ", ""),
    # Nucleases
    ("CDS_NUCLEASE_SPCAS9", "CDS_NUCLEASE", ["CDS_CAS9", "CDS_SPCAS9", "CAS9"], 200, "SpCas9", "CDS_CAS9"),
    # PolyA signals
    ("POLYA_BGH", "POLYA_SIGNAL", ["POLYA_BGH", "BGH"], 50, "BGH polyA", ""),
    ("POLYA_SV40_LATE", "POLYA_SIGNAL", ["POLYA_SV40", "SV40_LATE"], 50, "SV40 late polyA", ""),
    ("POLYA_HGH", "POLYA_SIGNAL", ["POLYA_HGH", "HGH"], 50, "hGH polyA", ""),
    # Origins
    ("ORI_COLE1_PUC", "ORI_BACTERIAL", ["ORI_GENERIC", "ORI_COLE1", "ORI_PUC", "COLE1"], 100, "ColE1/pUC ori", "ORI_GENERIC"),
    ("ORI_P15A", "ORI_BACTERIAL", ["ORI_P15A", "P15A"], 100, "p15A ori", ""),
    ("ORI_F1_PHAGE", "ORI_PHAGE", ["ORI_F1", "F1_ORI"], 100, "f1 phage ori", "ORI_F1"),
    # Terminators
    ("TERMINATOR_BACT_RRNB_T1T2", "TERMINATOR_BACTERIAL", ["TERMINATOR_RRNB_T1", "TERMINATOR_RRNB", "RRNB"], 50, "rrnB T1/T2", "TERMINATOR_RRNB_T1"),
    ("TERMINATOR_BACT_T7", "TERMINATOR_BACTERIAL", ["TERMINATOR_T7"], 20, "T7 terminator", ""),
    # Viral elements
    ("VIRAL_ELEMENT_LENTI_WPRE", "VIRAL_ELEMENT_LENTI", ["LENTI_WPRE", "LENTI_ELEMENT_WPRE", "WPRE"], 200, "WPRE", "LENTI_WPRE"),
    ("VIRAL_ELEMENT_LENTI_HIV_RRE", "VIRAL_ELEMENT_LENTI", ["LENTI_RRE", "LENTI_ELEMENT_RRE", "RRE"], 100, "HIV RRE", "LENTI_RRE"),
    ("VIRAL_ELEMENT_LENTI_HIV_5LTR", "VIRAL_ELEMENT_LENTI", ["LENTI_LTR_5", "LENTI_LTR", "LENTI_ELEMENT_5_LTR", "5_LTR"], 100, "HIV 5' LTR", "LENTI_LTR_5"),
    ("VIRAL_ELEMENT_LENTI_HIV_3LTR_DU3", "VIRAL_ELEMENT_LENTI", ["LENTI_LTR_3", "LENTI_ELEMENT_3_LTR", "3_LTR"], 100, "HIV 3' LTR dU3", "LENTI_LTR_3"),
    ("VIRAL_ELEMENT_AAV_ITR_AAV2", "VIRAL_ELEMENT_AAV", ["AAV_ITR", "MISC_3_itr", "ITR"], 50, "AAV2 ITR", "MISC_3_itr"),
    # Introns
    ("INTRON_IME_CHIMERIC", "INTRON_IME", ["INTRON_CHIMERIC", "CHIMERIC_INTRON"], 50, "chimeric intron", ""),
    # gRNA scaffold
    ("GRNA_SCAFFOLD_SPCAS9", "GRNA_SCAFFOLD", ["MISC_GRNA_SCAFFOLD", "GRNA_SCAFFOLD"], 20, "SpCas9 sgRNA scaffold", ""),
    # Signals
    ("MISC_SIGNAL_KOZAK", "MISC_SIGNAL", ["MISC_KOZAK", "KOZAK"], 5, "Kozak", ""),
    ("MISC_SIGNAL_IRES", "MISC_SIGNAL", ["MISC_IRES", "IRES"], 100, "IRES", ""),
    ("MISC_SIGNAL_T2A", "MISC_SIGNAL", ["LINKER_2A", "CDS_T2A", "T2A"], 30, "T2A", "LINKER_2A"),
    ("MISC_SIGNAL_P2A", "MISC_SIGNAL", ["LINKER_2A", "CDS_P2A", "P2A"], 30, "P2A", "LINKER_2A"),
    # NLS
    ("CDS_TAG_SV40_NLS", "CDS_TAG", ["CDS_SV40_NLS", "SV40_NLS"], 10, "SV40 NLS", "CDS_SV40_NLS"),
    # Additional bacterial elements
    ("PROMOTER_BACT_TAC", "PROMOTER_BACTERIAL", ["PROMOTER_TAC", "TAC"], 15, "tac promoter", ""),
    ("PROMOTER_BACT_ARABAD", "PROMOTER_BACTERIAL", ["PROMOTER_ARABAD", "ARABAD"], 100, "araBAD", ""),
]

MAX_SEED_LENGTH_BY_CATEGORY: dict[str, int] = {
    "PROMOTER_POL2": 2500,
    "PROMOTER_POL3": 500,
    "PROMOTER_BACTERIAL": 250,
    "CDS_RESISTANCE": 4000,
    "CDS_REPORTER": 4000,
    "CDS_NUCLEASE": 5000,
    "POLYA_SIGNAL": 1000,
    "ORI_BACTERIAL": 2500,
    "ORI_PHAGE": 2500,
    "TERMINATOR_BACTERIAL": 500,
    "VIRAL_ELEMENT_LENTI": 2000,
    "VIRAL_ELEMENT_AAV": 500,
    "INTRON_IME": 1500,
    "GRNA_SCAFFOLD": 500,
    "MISC_SIGNAL": 1000,
    "CDS_TAG": 250,
}

MAX_SEED_LENGTH_OVERRIDES: dict[str, int] = {
    "PROMOTER_CAG": 2500,
    "PROMOTER_SV40_EARLY": 1000,
    "PROMOTER_CMV_ENH_PROM": 1200,
    "PROMOTER_RSV_LTR": 1000,
    "PROMOTER_TRE3G": 600,
    "ORI_COLE1_PUC": 1200,
    "ORI_P15A": 1500,
    "ORI_F1_PHAGE": 1200,
    "VIRAL_ELEMENT_LENTI_HIV_5LTR": 1200,
    "VIRAL_ELEMENT_LENTI_HIV_3LTR_DU3": 1200,
    "VIRAL_ELEMENT_AAV_ITR_AAV2": 250,
    "TERMINATOR_BACT_RRNB_T1T2": 300,
    "POLYA_SV40_LATE": 500,
}

EXPECTED_SEED_TYPES: dict[str, dict[str, tuple[str, ...]]] = {
    "PROMOTER_POL2": {
        "canonical_types": ("promoter", "enhancer"),
        "feature_types": ("promoter", "enhancer"),
    },
    "PROMOTER_POL3": {
        "canonical_types": ("promoter",),
        "feature_types": ("promoter",),
    },
    "PROMOTER_BACTERIAL": {
        "canonical_types": ("promoter",),
        "feature_types": ("promoter",),
    },
    "CDS_RESISTANCE": {
        "canonical_types": ("cds", "protein", "gene"),
        "feature_types": ("cds", "gene"),
    },
    "CDS_REPORTER": {
        "canonical_types": ("cds", "protein", "gene"),
        "feature_types": ("cds", "gene"),
    },
    "CDS_NUCLEASE": {
        "canonical_types": ("cds", "protein", "gene"),
        "feature_types": ("cds", "gene"),
    },
    "CDS_TAG": {
        "canonical_types": ("cds", "protein", "gene"),
        "feature_types": ("cds", "gene"),
    },
    "POLYA_SIGNAL": {
        "canonical_types": ("polya_signal", "terminator"),
        "feature_types": ("polyA_signal", "terminator"),
    },
    "ORI_BACTERIAL": {
        "canonical_types": ("origin",),
        "feature_types": ("rep_origin",),
    },
    "ORI_PHAGE": {
        "canonical_types": ("origin",),
        "feature_types": ("rep_origin",),
    },
    "TERMINATOR_BACTERIAL": {
        "canonical_types": ("terminator",),
        "feature_types": ("terminator",),
    },
    "VIRAL_ELEMENT_LENTI": {
        "canonical_types": ("lentiviral_element", "misc"),
        "feature_types": ("misc_feature", "LTR", "repeat_region"),
    },
    "VIRAL_ELEMENT_AAV": {
        "canonical_types": ("misc",),
        "feature_types": ("repeat_region", "misc_feature"),
    },
    "INTRON_IME": {
        "canonical_types": ("intron",),
        "feature_types": ("intron",),
    },
    "GRNA_SCAFFOLD": {
        "canonical_types": ("guide_scaffold", "misc_rna"),
        "feature_types": ("misc_RNA",),
    },
    "MISC_SIGNAL": {
        "canonical_types": ("misc", "kozak_element"),
        "feature_types": ("misc_feature", "RBS"),
    },
}


def connect_db():
    try:
        import psycopg
        return psycopg.connect(DATABASE_URL)
    except ImportError:
        try:
            import psycopg2
            return psycopg2.connect(DATABASE_URL)
        except ImportError:
            log.error("Neither psycopg3 nor psycopg2 available")
            sys.exit(1)


def _revcomp(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ATGCN", "TACGN"))[::-1]


def _extract_feature_seq(module_seq: str, start: int, end: int, strand: int) -> list[str]:
    candidates: list[str] = []
    if 0 <= start < end <= len(module_seq):
        candidates.append(module_seq[start:end])
        if end < len(module_seq):
            candidates.append(module_seq[start:end + 1])
    if 1 <= start <= end <= len(module_seq):
        candidates.append(module_seq[start - 1:end])
    deduped: list[str] = []
    seen: set[str] = set()
    for seq in candidates:
        if strand == -1:
            seq = _revcomp(seq)
        if seq and seq not in seen:
            seen.add(seq)
            deduped.append(seq)
    return deduped


def _seed_type_score(category: str, canonical_type: str | None, feature_type: str | None) -> int:
    hints = EXPECTED_SEED_TYPES.get(category, {})
    canonical_type = (canonical_type or "").lower()
    feature_type = (feature_type or "").lower()
    score = 0
    if canonical_type and canonical_type in {v.lower() for v in hints.get("canonical_types", ())}:
        score += 2
    if feature_type and feature_type in {v.lower() for v in hints.get("feature_types", ())}:
        score += 2
    return score


def _select_candidate(candidates: list[dict]) -> tuple[str, str] | None:
    if not candidates:
        return None
    best_type_score = max(c["type_score"] for c in candidates)
    typed = [c for c in candidates if c["type_score"] == best_type_score]
    median_len = statistics.median(c["length_bp"] for c in typed)
    typed.sort(
        key=lambda c: (
            abs(c["length_bp"] - median_len),
            -c["type_score"],
            c["length_bp"],
        )
    )
    best = typed[0]
    return best["sequence"], best["canonical_id"]


def max_seed_length(pipeline_canonical_id: str, category: str) -> int:
    return MAX_SEED_LENGTH_OVERRIDES.get(
        pipeline_canonical_id,
        MAX_SEED_LENGTH_BY_CATEGORY.get(category, 5000),
    )


def query_module_sequence(
    conn,
    *,
    pipeline_canonical_id: str,
    category: str,
    search_terms: list[str],
    min_length: int,
) -> tuple[str, str] | None:
    """
    Search for a feature matching any of the search terms and return the feature span.

    Seed extraction must use the annotated feature coordinates from module_features,
    not the full module sequence. Otherwise promoters/origins inherit whole plasmid-
    length modules and poison the seed set.
    """
    cur = conn.cursor()
    max_length = max_seed_length(pipeline_canonical_id, category)
    for term in search_terms:
        candidates: list[dict] = []
        cur.execute(
            """
            SELECT
                m.sequence,
                mf.canonical_id,
                mf.canonical_type,
                mf.feature_type,
                mf.start,
                mf."end",
                COALESCE(mf.strand, 1)
            FROM modules m
            JOIN module_features mf ON mf.module_id = m.id
            WHERE LOWER(mf.canonical_id) = LOWER(%s)
              AND mf.start IS NOT NULL
              AND mf."end" IS NOT NULL
            LIMIT 100
            """,
            (term,),
        )
        for row in cur.fetchall():
            module_seq, canonical_id, canonical_type, feature_type, start, end, strand = row
            if start is None or end is None or start < 0 or end <= start or end > len(module_seq):
                continue
            type_score = _seed_type_score(category, canonical_type, feature_type)
            for feature_seq in _extract_feature_seq(module_seq, start, end, strand):
                if min_length <= len(feature_seq) <= max_length:
                    candidates.append(
                        {
                            "sequence": feature_seq,
                            "canonical_id": canonical_id,
                            "length_bp": len(feature_seq),
                            "type_score": type_score,
                        }
                    )
        selected = _select_candidate(candidates)
        if selected:
            return selected

    return None


def main():
    log.info("Connecting to %s", DATABASE_URL.split("@")[-1])
    conn = connect_db()

    seed_records = []
    n_found = 0
    n_missing = 0

    for pipeline_cid, category, search_terms, min_len, display_name, db_hint in SEED_TARGETS:
        log.info("Querying for %s (%s)...", pipeline_cid, category)
        result = query_module_sequence(
            conn,
            pipeline_canonical_id=pipeline_cid,
            category=category,
            search_terms=search_terms,
            min_length=min_len,
        )
        if result:
            seq, db_cid = result
            n_found += 1
            log.info("  OK: %s len=%dbp (db: %s)", pipeline_cid, len(seq), db_cid)
            seed_records.append({
                "pipeline_canonical_id": pipeline_cid,
                "category": category,
                "display_name": display_name,
                "db_canonical_id": db_cid,
                "sequence": seq,
                "length_bp": len(seq),
            })
        else:
            n_missing += 1
            log.warning("  MISS: %s — no match for %s", pipeline_cid, search_terms)

    conn.close()

    # Write FASTA
    HERE.mkdir(parents=True, exist_ok=True)
    with open(SEED_FASTA, "w") as f:
        for rec in seed_records:
            seq_id = f"SPLICIFY_SEED_{rec['pipeline_canonical_id']}"
            desc = f"category={rec['category']} display_name={rec['display_name']}"
            f.write(f">{seq_id} {desc}\n")
            seq = rec["sequence"]
            for i in range(0, len(seq), 70):
                f.write(seq[i:i+70] + "\n")

    log.info("Wrote %d seed sequences to %s", len(seed_records), SEED_FASTA)

    # Write metadata
    meta = {
        "seed_targets_total": len(SEED_TARGETS),
        "found": n_found,
        "missing": n_missing,
        "records": seed_records,
    }
    with open(SEED_META, "w") as f:
        json.dump(meta, f, indent=2)
    log.info("Wrote seed metadata to %s", SEED_META)

    if n_missing > 0:
        log.warning("%d targets had no DB match — add manually to %s", n_missing, SEED_FASTA.name)

    return 0


if __name__ == "__main__":
    sys.exit(main())
