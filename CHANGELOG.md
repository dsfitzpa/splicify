# Changelog

All notable changes to Splicify will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] — 2026-06-XX (planned)

### Added
- Initial public release as a monorepo at `github.com/dsfitzpa/splicify` (Apache-2.0).
- `splicify_api` Python package with a 15-symbol public API surface (`__init__.py`).
- `splicify_api._data` resolver for reference data (`$SPLICIFY_DATA_DIR` → `~/.splicify/data` → bundled test fixtures fallback).
- `scripts/fetch_data.sh` — Zenodo downloader with SHA256 verification.
- `scripts/rebuild_module_library.py` — strip SnapGene-curated annotations from CC-licensed plasmid sequences and re-annotate via the MIT clean-room `feature_annotator.py`.
- `scripts/regenerate_lm_descriptions.py` — regenerate Plasmid LM (description, target_tokens) pairs against the current annotation pipeline.
- `data/MANIFEST.yml` — pinned reference data with SHA + Zenodo DOI per artefact.
- `docs/ARCHITECTURE.md` — contributor-facing system map (replaces the internal build journals).
- `docs/recipes/` — 8 worked examples (Gibson, Gateway, Golden Gate, restriction, SDM, sgRNA Golden Gate, Cas9 sgRNA on genomic locus, pegRNA).
- `LICENSE` (Apache-2.0), `DECISIONS.md`, `SECURITY.md`, `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`.

### Changed
- Backend imports `splicify_api` as an installable package (`pip install -e backend/`). The previous `sys.path` shim that injected `the source backend repo` into the agent's `sys.path` is removed.
- Both production systemd units drop `Environment=PYTHONPATH=...` and rely on the venv-installed package.
- Reference data (`feature_db_data/`, `Module_Library_gb/`) no longer ships inside the wheel; users run `scripts/fetch_data.sh` to populate `$SPLICIFY_DATA_DIR`.
- ~11 hardcoded `Path(__file__).parent / "feature_db_data"` (and `parent.parent`, and `Path("/root/...")`) literals across 12 files migrated to `_data.data_path(...)`.

### Removed
- Vestigial `Environment="DATABASE_URL=..."` from the production systemd unit (0 live consumers found across the repo; Postgres databases remain on disk for now).
- `psycopg` + `psycopg-binary` runtime dependency (followed `DATABASE_URL` out).
- 33 backup files + 1 `.broken` file + the dead parallel `~/python-libraries/splicify_api/` directory (PR #2; preserved in git history of the source private repos).
- The 31,662 stale Claude-Haiku-generated (description, target_tokens) pairs (regenerated against the current pipeline by `scripts/regenerate_lm_descriptions.py`).

### Security
- All in-line `Environment="ANTHROPIC_API_KEY=..."` / `OPENAI_API_KEY` / `NCBI_API_KEY` / `DATABASE_URL` values rotated as part of pre-release scrub.
- Both GitHub-webhook HMAC secrets rotated.
- Migration to `EnvironmentFile=/etc/splicify/api.env` in PR #1 so future rotations don't touch the unit file or leak into systemd backups.

### Source repos archived
- `dsfitzpa/splicify-agent` (backend) → `dsfitzpa/splicify-agent-legacy`.
- `dsfitzpa/agent-v2-api` → `dsfitzpa/agent-v2-api-legacy`.
- `dsfitzpa/AI-Plasmid-Design` (frontend) → `dsfitzpa/AI-Plasmid-Design-legacy`.
- `dsfitzpa/splicify` (the n8n-era Next.js frontend, ~5.9 MB) → `dsfitzpa/splicify-legacy-n8n`.

All legacy repos kept for 24h grace, then made read-only / hard-deleted.
