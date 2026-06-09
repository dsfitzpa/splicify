# Locked decisions

A running log of decisions that contributors should NOT re-litigate. If you disagree, open an issue — don't quietly change the implementation.

## v0.1.0 (pre-release, 2026-05-29 → ongoing)

### Project name + repo
- **Name:** Splicify (the agent + paper + repo). "AI Plasmid Design" remains the consumer-UI brand at aiplasmiddesign.com.
- **Repo:** monorepo at `github.com/dsfitzpa/splicify` (public, Apache-2.0).
- **Why:** single name reduces brand fragmentation; monorepo eliminates the `sys.path` shim between the backend and the agent.

### License
- **Code:** Apache License 2.0.
- **Paper + manifests:** CC-BY-4.0.
- **Reference data:** mixed per tier (MIT / CC0 / CC-BY-4.0) — see `data/MANIFEST.yml`.
- **Why:** Apache patent grant matters for a CRISPR-design tool. MIT loses the patent grant; AGPL would scare off pharma/biotech contributors.

### LLM backend
- **First-class:** Claude (Anthropic SDK).
- **Behind a protocol:** `LLMClient` (TBD) lets cheaper / local backends slot in later without rewriting the agent loop.
- **Why:** Claude Sonnet 4.6 is the model the system was tuned against; backend-agnostic now would over-generalise.

### Data hosting
- **Zenodo** for the three large data deposits, with DOIs minted before v0.1.0 launch. Manifests + SHAs ship in the repo.
- **Why:** cite-able artefacts; Zenodo handles long-term hosting; SHAs let users verify offline.

### Plasmid LM Haiku-generated descriptions
- **Decision:** ship the **regenerator script** (`scripts/regenerate_lm_descriptions.py`), not the original 31,662 pairs.
- **Why:** the annotation pipeline changed since the 2026-04 descriptions were generated; they would be stale anyway. Users regenerate locally against the current pipeline on their own Anthropic API key.

### Module_Library_gb redistribution
- **Decision (2026-06-05):** ship `scripts/rebuild_module_library.py` in v0.1.0 (default mode: local-strip; `--verify-source` opt-in re-fetches from snapgene.com/plasmids). SnapGene-curated annotations stripped; all annotations re-emitted via the MIT clean-room `feature_annotator.py`. Output is a clean derivative work and is redistributed on Zenodo.
- **Why:** SnapGene's plasmid library sequences are CC-licensed; their curated annotations are proprietary. Stripping + re-annotating produces a defensible clean release. Public LM corpus stays at ~7,338 plasmids (essentially original size).

### Cutover pattern (PR #1 + monorepo creation)
- **Bundled.** PR #1 (splicify-core extraction) lands as the initial commit in the new public `dsfitzpa/splicify`. Three private source repos (`splicify-agent`, `agent-v2-api`, `AI-Plasmid-Design`) get renamed to `*-legacy-*` (24h grace, not hard-deleted).
- **Why:** one cutover window, one "main now lives at github.com/dsfitzpa/splicify" moment, less coordination overhead.

## v0.2.0 (planned)

- Resume the halted NCBI engineered-plasmid fetch (4,357 remaining → public LM corpus to ~9,846).
- Migrate every cloning intent to `CanonicalCloningRequest` (currently only restriction is done; 7 to go).
- LLMClient protocol (cheap-LLM-for-reviewers, local-OSS-model evaluations).

## Open questions (NOT decided yet)

- **PyPI publish.** v0.1.0 is install-from-source only. Earliest opportunity is when the public API surface (`splicify_api/__init__.py`) is locked across at least two releases.
- **Wet-lab validation partner.** Paper credibility wants 5–10 agent-designed plasmids actually cloned. Lab partner TBD.
- **Whether to keep the v1 Classic backend long-term.** Currently the frontend toggles between Classic (v1) and Agent V2; if v2 stabilises, Classic could be deprecated.
