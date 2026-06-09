# Splicify — Architecture

**Audience:** new contributor reading this for ~20 minutes before opening a PR.
**Not in scope:** every refactor since 2026-04. For history, see `docs/internal/HISTORY.md` (the consolidated `PROJECT_SUMMARY.md` + `agent_v2_summary.md` build journals).
**Companion docs:** `docs/AGENT_DESIGN.md` for tool roster + system prompts; `docs/internal/LLM_ANNOTATION_WORKFLOW.md` for pipeline internals; `docs/recipes/` for 8 worked examples.

---

## 1. What Splicify is

Splicify takes a natural-language prompt + optional GenBank uploads, and returns:

- an assembled-plasmid map (`assembled.gb`), parts order (`parts_order.csv`), wet-lab protocol (`protocol.csv`), workflow trace (`workflow_trace.txt`), and a plain-language reply — for plasmid cloning
- a guide-RNA design table (`guides.csv`), annotated genomic locus (`guides.gb`), plus the same parts order + protocol + trace — for CRISPR (sgRNA, pegRNA, primers)

It is a hierarchical Claude agent on top of a deterministic 6-stage annotation pipeline and a 9-workflow cloning dispatcher. The annotation + dispatcher layers also run as a plain REST API; the agent uses them as tools.

---

## 2. The 30-second mental model

```
┌────────────────────────────────────────────────────────────┐
│  USER  (Next.js UI + chat composer + GenBank upload)       │
└────────────────────────────┬───────────────────────────────┘
                             │
        ┌────────────────────┴────────────────────┐
        │                                         │
        ▼                                         ▼
┌────────────────────┐                  ┌───────────────────────┐
│  POST /api/chat    │                  │  POST /agent_v2/chat   │
│  (Classic backend) │                  │  (Hierarchical agent)  │
└─────────┬──────────┘                  └──────────┬────────────┘
          │                                        │
          │   deterministic dispatch               │  triage → 3 Explore (parallel)
          │   (intent → 1 of 9 handlers)           │  → Plan → Main ReAct → Summarizer
          │                                        │
          ▼                                        ▼
┌──────────────────────────────────────────────────────────────┐
│  Splicify core (Python package `splicify_api`)               │
│  • 6-stage annotation pipeline                               │
│  • 9 cloning operators (Gibson / Gateway / Golden Gate / …)  │
│  • PartResolver / TargetPlasmidBuilder / CloningRouter       │
│  • Verifier + KB-part auto-orientation correction            │
│  • 21-tool agent roster (annotate, simulate, design, emit)   │
└──────────────────────────────────────────────────────────────┘
```

Three layers:

1. **Annotation pipeline (deterministic).** Raw DNA → typed features → modules → SBO-typed interactions. Steps 1–2.75; see §4.
2. **Workflow dispatcher (deterministic).** Intent classification + 9 cloning operators + auto-routing across methods. See §5.
3. **Agent (LLM-driven).** Claude tool-use loop orchestrating layers 1+2 via a 21-tool roster. See §6.

Layers 1 and 2 are usable standalone; the agent depends on both.

---

## 3. A single request, end-to-end

For an Agent-V2 request like _"Design a Cas9 sgRNA targeting residue 33 of KEAP1 and primer pairs for NGS + Sanger verification"_ with `KEAP1.gb` attached:

1. `POST /agent_v2/chat` (multipart). Router builds an `AttachmentRegistry`, classifies each upload as plasmid / genomic via `file_kind.py`, eagerly annotates plasmids on a background task (cache-warming), stashes `(kind, gb_text)` in `attachment_kinds`.
2. **Triage** (1 Sonnet call). Returns `intent=CRISPR_GUIDE`, shorthand summary, `is_new_topic`.
3. **Two Explore subagents run in parallel** (`asyncio.gather`, fresh contexts, narrow tool subsets):
   - `TargetLocator` (annotate_attachment + resolve_feature_position) → walks the KEAP1 spliced CDS to find residue 33's genomic coords.
   - `GuideStrategist` (analyze_design_intent + lookup_kb_part) → picks Cas9 NGG / Doench 2014 / NGS + Sanger readout.
4. **Plan** (1 Sonnet call, no tools, sees only the digested ExploreFindings) → writes `plan.md` checklist.
5. **Main ReAct loop** (up to 24 iterations, full 21-tool roster, plan.md crossing-off):
   - `design_guides` → 5 sgRNA candidates with Doench scores + off-target predictions
   - `design_primers(application=illumina)` → NGS amplicon primers
   - `design_primers(application=sanger)` → Sanger amplicon primers
   - `emit_guides_csv` → guides.csv (5 guides + 4 primer pairs + 2 cloning oligos)
   - `emit_guides_gb` → guides.gb (KEAP1 source preserved + sgRNA + primer_bind features appended)
   - `emit_parts_order`, `emit_protocol` → 2 more files (dispatcher backfills args from the cache)
   - `auto_emit_workflow_trace` runs server-side after Main returns (avoids Vercel's 5-min SSE timeout)
6. **Summarizer** (1 Sonnet call) → polishes the Main draft into a ≤350-word reply.
7. **Envelope** returned: `{ok, reply, files, viz, viz_list, agent_trace, session_id, intent, workflow}`.

Steady-state cost on this prompt: ~3 min wall clock, ~8 tool calls, ~188 KB envelope.

---

## 4. Annotation pipeline (6 stages, deterministic)

`POST /plannotate/annotate_sequence_llm` is the engine every other system reads from.

| Step | Module | Emits |
|---|---|---|
| 1 | `feature_annotator.py` (clean-room, MIT, GenoLIB-seeded) | Flat features across 6 reference tiers: feature_reference (1,062 GenoLIB + RefSeq), feature_motifs (224 short motifs), feature_protein (706 GenoLIB CDS), fpbase (721 FPs), swissprot (66,221 PE-1), Rfam_curated (1,737 families) |
| 2 / 2a / 2b | `orf_finder.py` + `module_extractor.resolve_cds_submodules` | `cds_orf` parents + protein / NLS / tag / linker / 2A / gap submodules. Gap-fill scans hardcoded constants AND ~80 motifs from `feature_motifs_kb.json` (P2A/T2A/E2A/F2A in 18-aa "to-G" form, SV40 NLS variants, 9 FLAG variants, etc.). Iterative re-scan splits residual gaps until no more motifs found |
| 2.5 | `rule_based_module_detector.py` (~2,884 LOC, 55+ rules) | Boundary-defined modules: lentiviral / AAV / T-DNA / Gateway / floxed / FRT / transposons / ori / selection / lacZα blue-white / Pol III guide cassettes / etc. |
| 2.6 | `mammalian_pol2_detector.py` | Pol II expression cassettes + lentiviral 3-module split (upstream / payload / downstream) with Gen 1/2/3 classification |
| 2.75 | `cloning_feature_annotator.py` (~917 LOC) | Type II / Type IIs restriction sites + Gateway att sites + PCR-feasibility warnings on `layer="cloning_feature"`. ~0.09 s on a 15 kb plasmid |
| 3 (disabled) | `llm_module_parser.py` | High-level LLM module ID — off since 2026-04-15 |

After Step 6 every annotation is enriched with `so_role` / `sbo_role` / `sbo_participation` URIs (`so_sbo_mapping.py`), and `interaction_builder.py` emits SBO-typed `Interaction` records consumed by `sbol_io.py` for SBOL3 round-trip.

**Caching contract.** `annotation_cache.py` keys by SHA1(sequence) + depth (`full` vs `modules_only`); `annotate_cached` hits the hierarchy endpoint, `annotate_llm_cached` hits the LLM endpoint (with rule-based modules + interactions). The same sequence is never annotated twice in a process.

---

## 5. Cloning dispatcher (`/api/chat`, deterministic)

`chat.py` classifies user intent then runs one of 9 handlers:

`annotate_gb`, `gateway_cloning`, `gibson_design`, `plasmid_design`, `sdm_design`, `sgrna_golden_gate`, `golden_gate_primer_design`, `restriction_cloning`, `unknown`.

Three layers of pre-dispatch processing, in order:

1. **Intent classifier** (`intent.py`) — deterministic regex + keyword matcher. Pre-resolves KB-known parts via `extract_part_candidates` + `identify_features_from_kb`; downstream handlers read pre-resolved hits on `intent_result["kb_resolved"]`.
2. **Canonical-request normaliser** (`request_normalizer.py`) — collapses prompt + uploads into a `CanonicalCloningRequest`. Currently only restriction-cloning is migrated; other intents fall back to the legacy shape.
3. **Unified predesign pipeline** (`_execute_unified_predesign`):
   - `PartResolver` → `annotate_cached(part, "full")` per resolved part
   - `TargetPlasmidBuilder.build_from_parts()` → `annotate_cached(target, "modules_only")`
   - Build `PlasmidSpec` + `diff_spec_against_target`
   - `CloningRouter.route()` (cost / time / risk)
   - `target_from_inventory_router.verify_target_design()` on the assembled product
   - `auto_correct_kb_part_orientation()` if the verifier flags KB-part orientation issues
   - Optional `llm_orchestrator.review_design` (no-op unless `PLASMID_LLM_ORCHESTRATOR=1`)

**Auto-routing.** When target + inventory are provided and the prompt is method-agnostic, `target_from_inventory_router.py` runs 9 `assess_*` methods, each returning a uniform `FeasibilityReport`, sorts by `(score, -work_estimate, success_estimate)`, dispatches the winner. Always prepends a per-workflow scorecard to the reply via `build_audit_markdown()`.

---

## 6. The agent (`/agent_v2/chat`, hierarchical)

See `docs/AGENT_DESIGN.md` for the full design. Headline:

- **Triage** (Sonnet, no tools) → classifies into `PLASMID_CLONING | CRISPR_GUIDE | REJECT`.
- **Three Explore subagents per plasmid pipeline** (`PartScout` / `TargetBuilder` / `MethodRouter`) or **two per CRISPR pipeline** (`TargetLocator` / `GuideStrategist`), all running in parallel via `asyncio.gather`, each on a fresh context with a narrow tool subset.
- **Plan agent** (Sonnet, no tools) — sees only the digested `ExploreFinding` summaries; emits `plan.md`.
- **Main agent** — Sonnet ReAct loop with the full 21-tool roster, crossing items off `plan.md` after each tool call.
- **Summarizer** — polishes the Main draft into a ≤350-word reply.

**Tool roster (21 tools):**

| Category | Tools | Source |
|---|---|---|
| Annotation + KB | `annotate_attachment`, `lookup_kb_part`, `analyze_design_intent` | wraps `splicify_api` |
| Simulation | `simulate_assembly`, `golden_gate_assemble`, `digest_plasmid`, `find_primer_binding_sites` | wraps `splicify_api` |
| Verification + routing | `verify_assembly`, `route_workflow`, `score_sanger_primer`, `compare_to_choice` | wraps `splicify_api` |
| Design | `design_guides`, `design_pegrnas`, `design_primers`, `design_primers_batch` | wraps `guide_designer` / `pegrna_designer` / `pcr` |
| External | `find_genomic_record` (NCBI fetch), `find_external_part` (Addgene, capped at 6/session) | `external_search.py` + rate gate |
| Coords | `resolve_feature_position` (plasmid + genomic, walks spliced CDS) | `feature_resolver.py` |
| Output emitters | `emit_assembled_gb`, `emit_parts_order`, `emit_protocol`, `emit_guides_csv`, `emit_guides_gb` | `splicify_agent/outputs/` |

**Server-side trace emission.** `auto_emit_workflow_trace()` runs AFTER `run_main_agent` returns — the LLM does NOT call `emit_workflow_trace` directly. Avoids Vercel's 5-min SSE timeout on heavy runs.

---

## 7. Key invariants (do not violate)

1. **The agent never sees raw DNA.** `AttachmentRegistry` registers sequences server-side; every tool takes `attachment_id`. `_strip_sequences` recursively redacts any field named `sequence` / `seq` / `dna` / `template` before returning a tool result to the LLM. If you add a tool that needs DNA: resolve the attachment server-side and return only a digest.
2. **Every annotation is SHA1-keyed.** Adding a code path that bypasses `annotate_cached` / `annotate_llm_cached` defeats the cache and triggers re-annotation storms.
3. **Output envelope shape is stable.** `{ok, reply, files, viz, viz_list, agent_trace, session_id, intent, workflow}`. Adding fields is fine; renaming or removing requires a major version bump.
4. **Triage is deterministic for cache stability.** Prompts live in `prompts/` as `.md` files; temperature stays at 0; the prefix is wrapped in `cache_control` so first-token latency stays sub-second.
5. **No raw-DNA logging.** `SPLICIFY_LOG_SEQUENCES=false` is the default. Don't `logger.info(plasmid_seq)`.
6. **Subagent boundaries are explicit.** Subagents pass only their `ExploreFinding(role, summary_md, key_facts, references, trace)` digest up — never the raw tool transcripts. Keeps Plan + Main context windows lean.

---

## 8. Where to add things

| You want to add… | Touch… |
|---|---|
| A new restriction enzyme | `feature_db_data/feature_motifs_kb.json` + tests under `tests/test_cloning_feature_annotator.py` |
| A new motif (NLS / tag / linker / 2A) | `feature_db_data/feature_motifs_kb.json` — `_load_motif_kb_protein_catalog()` picks it up automatically |
| A new module type (e.g. PiggyBac cassette) | `rule_based_module_detector.py` (new `_detect_*` + entry in `_DETECTORS`) + `interaction_builder.py` (new `_MODULE_BUILDERS` entry) + add to `MODULE_VALIDATION_RULES` in `target_from_inventory_router.py` if it has interaction-graph constraints |
| A new cloning intent (e.g. RecET recombineering) | `intent.py` (regex + extractor) + `chat.py` (dispatcher branch) + new `<intent>_designer.py` + adapter in `workflow_input_adapters.py` + `assess_<intent>_feasibility` in `target_from_inventory_router.py` |
| A new agent tool | Handler module + schema in `splicify_agent/tools.py` + register in `EMITTER_HANDLERS` + update `make_full_tool_roster()` + (if relevant) update system prompt in `tool_schemas.py` to teach the LLM when to call it |
| A new output file (e.g. `vector_map.svg`) | New emitter in `splicify_agent/outputs/` + schema in `tools.py` + register in `dispatch_with_emitters` + add to canonical emitter sequence in the relevant system prompt |
| A new viewer track | `frontend/app/components/{Circular,Linear}*Viewer.tsx` + viz shape in `plannotate_router.py` |
| A new reference DB tier | Build script in `splicify_api/feature_db/`, ingestion in `feature_annotator.py`, `data/MANIFEST.yml` entry, Zenodo upload |
| A new Explore subagent | New module in `splicify_agent/subagents/` (use `run_explore_subagent` shared loop) + register in the relevant pipeline (`_run_plasmid_pipeline` or `_run_crispr_pipeline`) + add system prompt to `prompts/` |

---

## 9. Deterministic vs LLM (at a glance)

| Stage | Deterministic | LLM |
|---|---|---|
| Annotation pipeline (Steps 1–2.75) | ✅ | — |
| Intent classification | ✅ (regex) | — |
| KB part resolution | ✅ | — |
| PartResolver / TargetPlasmidBuilder / CloningRouter | ✅ | — |
| Target verification + KB-part auto-correct | ✅ | — |
| 9 cloning operators (Gibson, Gateway, …) | ✅ | — |
| Triage / Explore / Plan / Main / Summarizer | — | ✅ Claude Sonnet 4.6 |
| Tool execution (the 21 tools) | ✅ | — |
| Output emitters | ✅ | — |
| Auto workflow trace | ✅ (server-side) | — |

**Implication:** ~85% of the system is deterministic Python. The LLM does planning + tool selection + summarization. Adding a new feature usually means a deterministic Python module + telling the LLM about it via a tool schema + system prompt nudge.

---

## 10. Repo layout (post Phase 1)

```
splicify/
├── backend/
│   ├── pyproject.toml
│   ├── splicify_api/            # the public Python package
│   │   ├── feature_annotator.py        # Step 1
│   │   ├── orf_finder.py               # Step 2a
│   │   ├── rule_based_module_detector.py  # Step 2.5 (55+ rules)
│   │   ├── mammalian_pol2_detector.py  # Step 2.6
│   │   ├── cloning_feature_annotator.py # Step 2.75
│   │   ├── interaction_builder.py      # SBO-typed interactions
│   │   ├── sbol_io.py                  # SBOL3 export/import
│   │   ├── intent.py                   # deterministic intent classifier
│   │   ├── chat.py                     # /api/chat dispatcher (9 intents)
│   │   ├── annotation_cache.py         # SHA1-keyed cache
│   │   ├── predesign/                  # PartResolver / TargetPlasmidBuilder / CloningRouter
│   │   ├── cloning/                    # operators + designers + canonical_request
│   │   ├── target_from_inventory_router.py  # auto-routing + verifier
│   │   ├── agent/                      # 11-tool v1 agent (kept for comparison)
│   │   └── plasmid_lm/                 # generative plasmid language model
│   └── tests/
├── agent/
│   ├── pyproject.toml
│   ├── splicify_agent/          # was agent_v2
│   │   ├── orchestrator.py             # triage → Explore → Plan → Main → Summarizer
│   │   ├── triage.py
│   │   ├── main_agent.py               # ReAct loop + plan.md crossing-off
│   │   ├── subagents/                  # PartScout / TargetBuilder / MethodRouter / TargetLocator / GuideStrategist / Plan / Summarizer
│   │   ├── outputs/                    # 6 emitters (assembled_gb, parts_order, protocol, workflow_trace, guides_csv, guides_gb)
│   │   ├── tools.py                    # 21-tool schemas + dispatch_with_emitters
│   │   ├── crispr_tools.py             # design_guides / design_pegrnas / design_primers wrappers
│   │   ├── crispr_pipeline.py          # _run_crispr_pipeline
│   │   ├── file_kind.py                # plasmid vs genomic classifier
│   │   ├── genomic_annotator.py        # Biopython-driven CDS upgrade + spliced translation
│   │   ├── feature_resolver.py         # plasmid + genomic codon/aa-position resolver
│   │   ├── memory.py                   # Redis SessionState (30d TTL)
│   │   └── prompts/                    # versioned system prompts (cache-stable)
│   └── tests/
├── frontend/
│   ├── package.json
│   └── app/
│       ├── components/
│       │   ├── InteractiveSequenceViewer.tsx   # unified panel parent
│       │   ├── CircularPlasmidViewer.tsx       # ~1,748 LOC
│       │   ├── LinearSequenceViewer.tsx
│       │   └── Chat.tsx                        # composer + AI Agent toggle
│       └── api/
│           ├── chat/route.ts                   # proxy → /api/chat
│           └── agent_v2/chat/route.ts          # proxy → /agent_v2/chat
├── data/
│   ├── MANIFEST.yml             # SHA + Zenodo DOI per DB
│   └── README.md
├── deploy/
│   ├── vps/                     # systemd + Traefik + GitHub webhook
│   ├── docker/                  # docker-compose for local dev
│   └── vercel/                  # frontend deploy notes
├── docs/
│   ├── ARCHITECTURE.md          # this file
│   ├── AGENT_DESIGN.md          # tool roster + system prompts deep dive
│   ├── recipes/                 # 8 worked examples
│   └── internal/
│       └── HISTORY.md           # consolidated PROJECT_SUMMARY + agent_v2_summary
├── paper/                       # preprint source
└── scripts/
    ├── fetch_data.sh
    └── verify_data.sh
```

---

## 11. Cost + latency expectations

Measured on the KEAP1 R15C single-target run (see `docs/recipes/keap1-r15c.md`):

| Metric | Value |
|---|---|
| Wall clock | ~3 min |
| Tool calls (Main loop) | 6–10 |
| Files emitted | 5–6 |
| Envelope size | ~188 KB |
| Anthropic cost (Sonnet 4.6, with prompt caching) | ~$X.YZ |
| LAB-Bench CloningScenarios (v1 baseline) | 29/33 = 87.9% |
| LAB-Bench CloningScenarios (v2 — pending Phase 5 rerun) | TBD |

Prompt caching cuts per-iteration TPM ~3–4× on a 6-iteration loop. Explore phase is hard-capped at 120 s (subagents return placeholder findings on timeout). `find_external_part` is capped at 6 calls per session.

---

## 12. Reading order for new contributors

1. **This file** — overall map.
2. `docs/AGENT_DESIGN.md` — the 21-tool roster + system prompts + dependency-injection contract.
3. `docs/recipes/keap1-r15c.md` — one worked example end-to-end.
4. `docs/internal/LLM_ANNOTATION_WORKFLOW.md` — pipeline internals (if you're touching annotation).
5. `docs/internal/MODULE_DETECTION_AND_SBOL3.md` — module hierarchy + SBOL3 (if you're adding a module type).
6. `docs/internal/CLONING_WORKFLOWS.md` — per-intent handlers (if you're adding a cloning workflow).
7. `docs/internal/HISTORY.md` — why things are the way they are.

