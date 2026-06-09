# Splicify

**Hierarchical Claude agent for plasmid cloning and CRISPR guide design.** Takes a natural-language prompt + optional GenBank uploads; emits an assembled-plasmid map, parts order, wet-lab protocol, workflow trace, and a plain-language reply (or guides.csv / guides.gb for CRISPR).

```
prompt + .gb upload  →  Triage  →  Explore × 3 (parallel)  →  Plan  →  Main ReAct  →  Summarizer
                                                                              │
                                                                              ▼
                       assembled.gb · parts_order.csv · protocol.csv · workflow_trace.txt · reply
```

87.9% on the LAB-Bench CloningScenarios public split (v1; v2 rerun in progress). KEAP1 R15C end-to-end case study at `docs/recipes/keap1-r15c.md`.

## Status

v0.1.0 — pre-release. Preprint pending. Public API surface (15 symbols) is stable; everything else is internal until v1.0.

## Install

```bash
git clone https://github.com/dsfitzpa/splicify.git
cd splicify

# Backend + agent (one venv)
python3 -m venv .venv
source .venv/bin/activate
pip install -e backend/ -e agent/

# Reference data (635 MB + 128 MB + 268 MB from Zenodo)
SPLICIFY_DATA_DIR=~/.splicify/data ./scripts/fetch_data.sh

# Frontend (separate)
cd frontend && npm install && npm run dev
```

Env vars (see `.env.example`):
- `ANTHROPIC_API_KEY` (required for the agent)
- `NCBI_API_KEY` + `NCBI_EMAIL` (recommended for genomic record fetch)
- `REDIS_URL` (defaults to `redis://localhost:6379/1` — agent session state)
- `SPLICIFY_DATA_DIR` (defaults to `~/.splicify/data`)

## Local-dev quickstart

```bash
# Start the backend
uvicorn splicify_api.main:app --reload --port 8000

# Start the agent
uvicorn main:app --reload --port 8002 --app-dir agent/

# Run tests
pytest backend/tests/ agent/tests/
```

## Architecture

Three layers:

1. **Annotation pipeline** (`backend/splicify_api/`) — 6 deterministic stages turning raw DNA into typed features + modules + SBO-typed interactions.
2. **Cloning dispatcher** (`backend/splicify_api/chat.py`) — 9 cloning intents (Gibson, Gateway, Golden Gate, restriction, SDM, sgRNA Golden Gate, ...) with auto-routing.
3. **Hierarchical agent** (`agent/splicify_agent/`) — Triage → 3 Explore subagents (parallel) → Plan → Main ReAct → Summarizer. 21-tool roster.

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the full map. The annotation + dispatcher layers also run as a plain REST API; the agent uses them as tools.

## Citation

If you use Splicify in your research, please cite the preprint (link added on launch).

## License

[Apache License 2.0](LICENSE). The reference data deposits on Zenodo are CC-BY-4.0 / CC0 / MIT depending on tier — see `data/MANIFEST.yml`.

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md). Issues + PRs welcome.

## Acknowledgements

The annotation pipeline's Step 1 (`feature_annotator.py`) is a clean-room reimplementation seeded by [GenoLIB](https://github.com/...). Plasmid LM training corpus draws sequences from [SnapGene's CC-licensed plasmid sets](https://www.snapgene.com/plasmids) (excluding coronavirus subset), [Addgene](https://www.addgene.org/), [FPbase](https://www.fpbase.org/), [UniProt/SwissProt](https://www.uniprot.org/), and [Rfam](https://rfam.org/).
