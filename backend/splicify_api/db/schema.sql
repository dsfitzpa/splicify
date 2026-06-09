-- Module Library Schema v2 (Rebuild)
-- Designed for granular biological module representation with rich intent metadata.
-- Compatible with module_library_db.py query layer.

CREATE TABLE IF NOT EXISTS plasmids (
    id           TEXT PRIMARY KEY,
    source_group TEXT,
    filename     TEXT,
    path         TEXT,
    topology     TEXT,        -- circular | linear | unknown
    length       INTEGER,
    seq_hash     TEXT,
    metadata     JSONB        -- {vector_type, has_guide_cassette, cloning, source_vector_context, ...}
);

-- Core modules table. module_type is now fine-grained:
--   Atomic:    promoter_pol2, promoter_pol3, polya_signal, cds_reporter, cds_nuclease,
--              cds_selection_marker, cds_2a_peptide, cds_nls, cds_tag, cds_generic,
--              guide_scaffold, replication_origin, enhancer_element, kozak_element,
--              intron_element, lentiviral_element, misc_element
--   Composite: pol2_expression_cassette, pol3_guide_cassette, nuclease_expression_cassette,
--              reporter_cassette, selection_cassette, bacterial_backbone,
--              lentiviral_backbone, lentiviral_expression_vector, expression_vector
CREATE TABLE IF NOT EXISTS modules (
    id            TEXT PRIMARY KEY,
    plasmid_id    TEXT REFERENCES plasmids(id) ON DELETE CASCADE,
    module_type   TEXT NOT NULL,
    payload_id    TEXT,        -- primary canonical ID (e.g. CDS_EGFP, PROMOTER_CMV)
    start         INTEGER,
    "end"         INTEGER,
    wraps         BOOLEAN      DEFAULT FALSE,
    length        INTEGER      NOT NULL,
    seq_hash      TEXT,
    sequence      TEXT         NOT NULL,
    end_inferred  BOOLEAN      DEFAULT FALSE,
    -- Rich metadata (all extra fields live here for backward compat):
    --   source_vector_context: standard|lentiviral|aav|retroviral
    --   is_atomic: true|false
    --   is_composite: true|false
    --   biological_role: promoter|reporter|nuclease|selection_marker|origin|backbone|linker|guide_expression|nls|2a_peptide
    --   vector_context: [mammalian, bacterial, any, ...]
    --   has_promoter, has_polya, has_nls: bool
    --   plannotate_canonical_ids: [...]
    --   complementary_module_types: [...]
    metadata      JSONB        DEFAULT '{}'::jsonb
);

CREATE TABLE IF NOT EXISTS module_features (
    module_id      TEXT    REFERENCES modules(id) ON DELETE CASCADE,
    canonical_id   TEXT,
    canonical_type TEXT,
    original_name  TEXT,
    feature_type   TEXT,
    start          INTEGER,
    "end"          INTEGER,
    strand         INTEGER,
    qualifiers     JSONB,
    -- Source of annotation: genbank | plannotate
    annotation_source TEXT DEFAULT 'genbank',
    PRIMARY KEY (module_id, canonical_id, original_name, start, "end")
);

CREATE TABLE IF NOT EXISTS module_text (
    module_id    TEXT PRIMARY KEY REFERENCES modules(id) ON DELETE CASCADE,
    display_name TEXT,
    description  TEXT,
    tokens       TEXT[],
    embedding    FLOAT8[],
    tsv          tsvector
);

CREATE TABLE IF NOT EXISTS module_kmers (
    module_id TEXT    REFERENCES modules(id) ON DELETE CASCADE,
    k         INTEGER,
    kmers     BIGINT[],
    PRIMARY KEY (module_id, k)
);

-- Indexes for fast retrieval
CREATE INDEX IF NOT EXISTS idx_modules_type        ON modules(module_type);
CREATE INDEX IF NOT EXISTS idx_modules_type_payload ON modules(module_type, payload_id);
CREATE INDEX IF NOT EXISTS idx_modules_plasmid      ON modules(plasmid_id);
CREATE INDEX IF NOT EXISTS idx_module_features_mod  ON module_features(module_id);
CREATE INDEX IF NOT EXISTS idx_module_features_cid  ON module_features(canonical_id);
CREATE INDEX IF NOT EXISTS idx_module_text_tsv      ON module_text USING gin(tsv);
CREATE INDEX IF NOT EXISTS idx_modules_meta_ctx     ON modules USING gin(metadata jsonb_path_ops);
