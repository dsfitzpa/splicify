# Hierarchical Feature Rules Sheet

This sheet is a derived view of pLannotate's `snapgene.csv`. It is intended to
replace the current hardcoded "most common cassette" heuristics with a
feature-driven grammar table.

## Output

Generated CSV:
- `backend/splicify_api/annotator/feature_db/schema/hierarchical_feature_rules.csv`

Generator:
- `backend/splicify_api/annotator/feature_db/pipeline/10_build_hierarchical_feature_sheet.py`

## Column meanings

- `host_specificity`: `exclusive`, `multi_host`, or `general`
- `taxonomic_scope`: the primary host context used for module grammar
- `taxonomy_basis`: why that host scope was chosen
- `functional_bucket`: normalized feature role used by the annotator
- `nested_module_role`: how the feature behaves in module construction
  - `hard_start`: definitive module start
  - `fallback_start`: module start only when a stronger upstream boundary is absent
  - `hard_end`: definitive module end
  - `internal`: internal feature expected inside a module
- `starts_module_types`: module families this feature can start
- `ends_module_types`: module families this feature can end
- `preferred_parent_modules`: module families this feature is expected to live inside
- `required_partner_features`: expected neighboring features for a valid module
- `dependency_notes`: plain-language summary of the rule
- `boundary_priority`: ordering hint when multiple boundaries overlap
- `validation_rules`: sequence or ordering checks that should be enforced
- `formatting_rules`: how to serialize the finished higher-order module

## Intended grammar

- Animal, plant, fungal, and insect Pol II cassettes should run from enhancer or
  promoter to polyA or the relevant transcript-terminating feature.
- Pol III guide/shRNA cassettes should run from the Pol III promoter to the
  terminator and include spacer/payload plus scaffold when applicable.
- Lentiviral payloads should run from 5' LTR or `psi` to 3' LTR.
- T-DNA payloads should run from right border to left border.
- Backbone modules should be host-specific. `pVS1` and T-DNA features belong to
  Agrobacterium/binary-vector context and should not be collapsed into a generic
  bacterial backbone.
- CDS-derived payloads should only be promoted to trusted modules when the full
  ORF is present and strand-aware start/stop checks pass.
