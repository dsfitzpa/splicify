#!/usr/bin/env python3
"""
LLM plasmid + part annotator — pilot driver.

Reads:
  - /var/data/plasmid_lm_corpus/token_corpus.jsonl  (rotation-0 examples)
  - /var/data/plasmid_lm_corpus/part_cooccurrence.jsonl (for part-level prompts)

Writes:
  - plasmid_descriptions.jsonl  — one row per plasmid with 4 description styles + tags
  - part_annotations.jsonl      — one row per canonical part with short/long + typical hosts + use cases

Pilot mode (default --sample 50 plasmids, --part-sample 30 parts) to validate
prompts before full batch run. Uses the Anthropic Messages API directly; batch
API wiring added in a later pass.

Hallucination-grounding filter (PLASMID_LLM_ANNOTATION_PLAN.md §6):
  - Every canonical feature name mentioned in a description must exist in the
    plasmid's token stream. Up to 2 retries on violation.
"""
from __future__ import annotations
import argparse
import json
import logging
import os
import random
import re
import sys
import time
from pathlib import Path

logger = logging.getLogger("llm_annotator")

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

DEFAULT_MODEL = "claude-haiku-4-5"
ANTHROPIC_KEY_ENV = "ANTHROPIC_API_KEY"
MAX_RETRIES = 2
REQUEST_TIMEOUT = 120


# ---------------------------------------------------------------------------
# Prompt templates
# ---------------------------------------------------------------------------

PLASMID_PROMPT_SYSTEM = """You annotate plasmids for a machine-learning training corpus. \
You must stay grounded in the provided TOKEN STREAM: do not mention any \
feature or payload that does not appear in the tokens, and do not invent \
names. Produce valid JSON only, no surrounding prose."""

PLASMID_PROMPT_TEMPLATE = """TOKEN STREAM (brackets indicate module boundaries; role:payload is the part):
{tokens}

INTERACTIONS (functional edges detected between the tokens above):
{interactions}

METADATA:
host: {host}
length: {length} bp
source: {source}
topology: {topology}

Produce JSON with EXACTLY these keys:

{{
  "user_request_short":      "string, <= 8 words, imperative or noun-phrase request",
  "user_request_functional": "string, <= 25 words, specific payload names",
  "lab_slack_question":      "string, <= 30 words, casual conversational ask",
  "methods_spec":            "string, <= 50 words, technical paper-Methods-section tone",
  "tags":                    ["3-6 short lowercase tags: host, delivery, category, key markers"]
}}

Rules:
- Every feature you mention must correspond to a token in the stream (e.g. Cas9, EGFP, EF-1α, bGH, pUC, etc.).
- Prefer concrete payload names over generic role words.
- "methods_spec" may reference architectural facts visible in the tokens (Gen 3 SIN lentivirus, 2A-skip, WPRE, etc.).
- Lead with intent: what the plasmid does, not what it IS.

Return only the JSON object."""

PART_PROMPT_SYSTEM = """You annotate a single biological part (one canonical plasmid \
feature) for a machine-learning training corpus. You must stay grounded in the \
corpus observations provided. Do not invent mechanisms you are not confident about. \
Produce valid JSON only."""

PART_PROMPT_TEMPLATE = """Canonical part: {canonical_id}
Role:          {role}
Observed in {occurrence_count} plasmids ({corpus_coverage_pct:.1f}% of corpus)

Host distribution (how it's used in the training corpus):
{host_dist}

Top co-occurring features (most-frequent neighbors across plasmids):
{neighbors}

Typical enclosing modules:
{modules}

Example plasmids containing it:
{example_plasmids}

Produce JSON with EXACTLY these keys:

{{
  "short":                   "string, <= 15 words, practical description",
  "long":                    "string, <= 60 words, mechanism + typical context",
  "typical_hosts":           ["list of lowercase hosts from: mammalian, bacterial, insect, plant, yeast"],
  "common_pairings_observed":["5-10 neighbor payload names taken from the list above"],
  "use_cases":               ["3-5 short use-case phrases"],
  "tags":                    ["3-5 lowercase tags"]
}}

Rules:
- Describe MECHANISM before opinion. Use active voice.
- Reuse neighbor names provided — do not pull in unrelated parts.
- If the part name is unfamiliar, describe the role category generically rather than inventing.

Return only the JSON object."""


# ---------------------------------------------------------------------------
# Grounding filter
# ---------------------------------------------------------------------------

# Match only strings that LOOK like biological payload names, not ordinary English:
#   - contain at least one uppercase letter AND one other non-alpha-lower char
#     (digit OR uppercase OR hyphen), OR
#   - be >=3 chars all-uppercase
# This deliberately excludes plain lowercase English words ("annotate", "ready",
# "screening") which are valid prose, not payload references.
PAYLOAD_NAME_RE = re.compile(
    r"\b("
    r"[A-Z][A-Za-z0-9]*[0-9][A-Za-z0-9\-]*"       # has digit somewhere, starts cap (Cas9, SV40, EF1a)
    r"|[A-Z]{2,}[A-Za-z]*(?:-[A-Za-z0-9]+)*"       # all-caps run (CMV, EGFP, AmpR-010, MCS-328)
    r"|[a-z]+[A-Z][A-Za-z0-9\-]*"                   # camelCase starting lowercase (lacZ, mCherry, dCas9)
    r"|[A-Z][a-z]+[A-Z][A-Za-z0-9\-]*"              # CamelCase with multiple caps (BamHI, HindIII, EcoRI)
    r")\b"
)


def _norm(s: str) -> str:
    """Lowercase + strip trailing numeric suffixes like -009, -001, _(2)."""
    s = re.sub(r"[-_]\(?\d{1,4}\)?$", "", s)
    return s.strip().lower()


def extract_payload_tokens(tokens: list[str]) -> set[str]:
    """
    Payload vocabulary for the grounding filter: every name emitted as a
    feature, submodule, module_open value, or CLN payload. Each is stored
    both normalized and as its sub-pieces (split on / _ -).
    """
    out: set[str] = set()
    for t in tokens:
        if ":" not in t:
            continue
        body = t.strip("<>")
        parts = body.split(":", 2)  # allow CLN:unique_cutter:EcoRI
        if len(parts) < 2:
            continue
        role = parts[0]
        # Skip pure-metadata header tokens
        if role in ("BOS", "EOS", "TOPOLOGY", "LEN_BIN", "HOST", "SOURCE",
                    "ROTATION_IDX", "ROTATION_IDX"):
            continue
        # Candidate payload strings: take every non-initial segment
        for payload in parts[1:]:
            payload = payload.strip()
            if not payload:
                continue
            n = _norm(payload)
            if not n:
                continue
            out.add(n)
            for piece in re.split(r"[-_/ ]+", n):
                if len(piece) >= 3:
                    out.add(piece)
    return out


def check_grounding(description: str, payload_set: set[str]) -> list[str]:
    """
    Only flag strings that look like biological payload names AND don't match
    anything in the payload vocabulary (even via substring). Plain English words
    aren't matched by the regex in the first place.
    """
    violations: list[str] = []
    for m in PAYLOAD_NAME_RE.finditer(description):
        raw = m.group(1)
        if len(raw) < 3:
            continue
        norm = _norm(raw)
        if not norm or norm.isdigit():
            continue
        if norm in payload_set:
            continue
        # Substring either direction
        hit = False
        for p in payload_set:
            if norm in p or p in norm:
                hit = True
                break
        if hit:
            continue
        # Also match against common biological "class" terms that are acceptable
        # without being in the token vocabulary (host types, delivery types).
        CLASS_TERMS = {
            "lentiviral", "lentivirus", "retroviral", "retrovirus",
            "adenovirus", "aav", "baculovirus", "crispr",
            "mammalian", "bacterial", "insect", "plant", "yeast",
            "transgene", "selectable", "inducible", "constitutive",
            "bluescribe", "bluescript",
        }
        if norm in CLASS_TERMS:
            continue
        violations.append(raw)
    return violations


# ---------------------------------------------------------------------------
# Anthropic client
# ---------------------------------------------------------------------------

def anthropic_call(system: str, user: str, model: str) -> str | None:
    try:
        import anthropic
    except ImportError:
        logger.error("anthropic package not installed; run: pip install anthropic")
        return None
    key = os.environ.get(ANTHROPIC_KEY_ENV)
    if not key:
        logger.error("%s env var not set", ANTHROPIC_KEY_ENV)
        return None
    client = anthropic.Anthropic(api_key=key, timeout=REQUEST_TIMEOUT)
    try:
        msg = client.messages.create(
            model=model,
            max_tokens=800,
            system=system,
            messages=[{"role": "user", "content": user}],
        )
    except Exception as exc:
        logger.warning("Anthropic API error: %s", exc)
        return None
    text = "".join(
        block.text for block in msg.content if getattr(block, "type", "") == "text"
    )
    return text.strip()


def parse_json_reply(text: str) -> dict | None:
    if not text:
        return None
    # Strip code fences if model wrapped output
    t = text.strip()
    if t.startswith("```"):
        t = re.sub(r"^```(?:json)?\s*", "", t)
        t = re.sub(r"\s*```$", "", t)
    # Grab the outermost { ... } if there's leading/trailing prose
    first = t.find("{")
    last = t.rfind("}")
    if first == -1 or last == -1:
        return None
    try:
        return json.loads(t[first:last + 1])
    except json.JSONDecodeError:
        return None


# ---------------------------------------------------------------------------
# Prompt builders
# ---------------------------------------------------------------------------

def _compact_interactions(interactions: list[dict]) -> str:
    if not interactions:
        return "(none detected)"
    lines = []
    for ix in interactions[:12]:  # cap
        rid = ix.get("rule_id", "?")
        src = ix.get("source_module", "")
        parts = ", ".join(ix.get("participants") or [])
        lines.append(f"  {rid}  ({src}):  {parts}")
    return "\n".join(lines) or "(none detected)"


def _compact_tokens(tokens: list[str]) -> str:
    # Skip rotation index (per-call noise), keep everything else
    kept = [t for t in tokens if not t.startswith("<ROTATION_IDX:")]
    return " ".join(kept)


def build_plasmid_prompt(ex: dict) -> str:
    return PLASMID_PROMPT_TEMPLATE.format(
        tokens=_compact_tokens(ex["tokens"]),
        interactions=_compact_interactions(ex.get("interaction_tokens") or []),
        host=ex.get("host", "unknown"),
        length=ex.get("length", "?"),
        source=ex.get("source", "unknown"),
        topology=ex.get("topology", "circular"),
    )


def build_part_prompt(row: dict) -> str:
    host_dist = row.get("host_distribution") or {}
    host_lines = "\n".join(f"  {h}: {n}" for h, n in host_dist.items()) or "  (none)"
    nb = row.get("top_neighbors") or []
    nb_lines = "\n".join(f"  {n['canonical_id']} (×{n['count']})" for n in nb[:12]) or "  (none)"
    mods = row.get("top_enclosing_modules") or []
    mod_lines = "\n".join(f"  {m['module_type']} (×{m['count']})" for m in mods[:5]) or "  (none)"
    examples = row.get("example_plasmids") or []
    ex_lines = ", ".join(examples) or "(none)"
    return PART_PROMPT_TEMPLATE.format(
        canonical_id=row["canonical_id"],
        role=row["role"],
        occurrence_count=row["occurrence_count"],
        corpus_coverage_pct=row["corpus_coverage"] * 100,
        host_dist=host_lines,
        neighbors=nb_lines,
        modules=mod_lines,
        example_plasmids=ex_lines,
    )


# ---------------------------------------------------------------------------
# Per-item annotate loop
# ---------------------------------------------------------------------------

def annotate_plasmid(ex: dict, model: str) -> dict | None:
    payload_set = extract_payload_tokens(ex["tokens"])
    prompt = build_plasmid_prompt(ex)
    for attempt in range(1 + MAX_RETRIES):
        reply = anthropic_call(PLASMID_PROMPT_SYSTEM, prompt, model=model)
        parsed = parse_json_reply(reply or "")
        if not parsed:
            logger.warning("attempt %d: JSON parse failed for %s", attempt, ex["plasmid_id"])
            continue
        # Grounding check on the four description fields
        all_text = " ".join(
            (parsed.get(k, "") or "") for k in
            ("user_request_short", "user_request_functional",
             "lab_slack_question", "methods_spec")
        )
        violations = check_grounding(all_text, payload_set)
        if not violations:
            return {
                "plasmid_id": ex["plasmid_id"],
                "source_corpus": ex.get("source", "unknown"),
                "host": ex.get("host", "unknown"),
                "descriptions": [
                    {"style": "user_request_short",      "text": parsed.get("user_request_short", "")},
                    {"style": "user_request_functional", "text": parsed.get("user_request_functional", "")},
                    {"style": "lab_slack_question",      "text": parsed.get("lab_slack_question", "")},
                    {"style": "methods_spec",            "text": parsed.get("methods_spec", "")},
                ],
                "tags": parsed.get("tags", []),
                "hallucination_check": {
                    "attempts": attempt + 1,
                    "unknown_tokens_mentioned": 0,
                },
            }
        logger.info("attempt %d: %s had %d suspicious names: %s",
                    attempt, ex["plasmid_id"], len(violations),
                    ", ".join(sorted(set(violations))[:5]))
        prompt += (
            "\n\nPrevious attempt mentioned names NOT present in the tokens: "
            + ", ".join(sorted(set(violations))[:10])
            + ". Regenerate without referring to these."
        )
    logger.warning("FAILED: %s — grounding violations survived %d retries",
                   ex["plasmid_id"], MAX_RETRIES)
    return {
        "plasmid_id": ex["plasmid_id"],
        "source_corpus": ex.get("source", "unknown"),
        "host": ex.get("host", "unknown"),
        "descriptions": [],
        "tags": [],
        "hallucination_check": {
            "attempts": MAX_RETRIES + 1,
            "failed": True,
            "last_violations": sorted(set(violations))[:10],
        },
    }


def annotate_part(row: dict, model: str) -> dict | None:
    prompt = build_part_prompt(row)
    reply = anthropic_call(PART_PROMPT_SYSTEM, prompt, model=model)
    parsed = parse_json_reply(reply or "")
    if not parsed:
        return {"canonical_id": row["canonical_id"], "failed": True}
    return {
        "canonical_id": row["canonical_id"],
        "role": row["role"],
        "occurrence_count": row["occurrence_count"],
        "short": parsed.get("short", ""),
        "long": parsed.get("long", ""),
        "typical_hosts": parsed.get("typical_hosts", []),
        "common_pairings_observed": parsed.get("common_pairings_observed", []),
        "use_cases": parsed.get("use_cases", []),
        "tags": parsed.get("tags", []),
    }


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def load_resume(path: Path, key: str) -> set[str]:
    if not path.exists():
        return set()
    done: set[str] = set()
    with open(path) as fh:
        for line in fh:
            try:
                done.add(json.loads(line).get(key, ""))
            except Exception:
                continue
    return done


def iter_rotation0_plasmids(corpus: Path, skip: set[str], sample: int | None, seed: int):
    # First pass: list eligible plasmid_ids at rotation 0
    eligible: list[tuple[int, dict]] = []
    with open(corpus) as fh:
        for offset, line in enumerate(fh):
            try:
                ex = json.loads(line)
            except Exception:
                continue
            if ex.get("rotation_idx", 0) != 0:
                continue
            if ex["plasmid_id"] in skip:
                continue
            eligible.append((offset, ex))
    if sample and len(eligible) > sample:
        rng = random.Random(seed)
        eligible = rng.sample(eligible, sample)
    for _off, ex in eligible:
        yield ex


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/token_corpus.jsonl"))
    ap.add_argument("--part-stats", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/part_cooccurrence.jsonl"))
    ap.add_argument("--output-plasmid", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/plasmid_descriptions.jsonl"))
    ap.add_argument("--output-part", type=Path,
                    default=Path("/var/data/plasmid_lm_corpus/part_annotations.jsonl"))
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--sample-plasmids", type=int, default=50)
    ap.add_argument("--sample-parts", type=int, default=30)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--skip-plasmids", action="store_true")
    ap.add_argument("--skip-parts", action="store_true")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s - %(levelname)s - %(message)s")

    if not args.skip_plasmids:
        skip = load_resume(args.output_plasmid, "plasmid_id") if args.resume else set()
        logger.info("plasmid annotation: skipping %d already-done", len(skip))
        n_ok = n_fail = 0
        args.output_plasmid.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output_plasmid, "a" if args.resume else "w") as fh:
            for i, ex in enumerate(iter_rotation0_plasmids(
                    args.corpus, skip, args.sample_plasmids, args.seed)):
                t0 = time.time()
                row = annotate_plasmid(ex, args.model)
                if row is None:
                    n_fail += 1
                    continue
                if row.get("hallucination_check", {}).get("failed"):
                    n_fail += 1
                else:
                    n_ok += 1
                fh.write(json.dumps(row, ensure_ascii=False) + "\n")
                fh.flush()
                logger.info("[%d] %s — %.1fs ok=%d fail=%d",
                            i + 1, ex["plasmid_id"], time.time() - t0, n_ok, n_fail)
        logger.info("PLASMID PILOT DONE: ok=%d fail=%d", n_ok, n_fail)

    if not args.skip_parts:
        # Pick top-N most-frequent parts not already annotated
        skip = load_resume(args.output_part, "canonical_id") if args.resume else set()
        rows: list[dict] = []
        with open(args.part_stats) as fh:
            for line in fh:
                row = json.loads(line)
                if row["canonical_id"] in skip:
                    continue
                rows.append(row)
        rows.sort(key=lambda r: -r["occurrence_count"])
        pilot_rows = rows[:args.sample_parts] if args.sample_parts else rows
        logger.info("part annotation: %d rows in pilot (out of %d stats)", len(pilot_rows), len(rows))
        n_ok = n_fail = 0
        args.output_part.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output_part, "a" if args.resume else "w") as fh:
            for i, row in enumerate(pilot_rows):
                t0 = time.time()
                out = annotate_part(row, args.model)
                if out is None or out.get("failed"):
                    n_fail += 1
                else:
                    n_ok += 1
                fh.write(json.dumps(out, ensure_ascii=False) + "\n")
                fh.flush()
                logger.info("[%d] %s — %.1fs ok=%d fail=%d",
                            i + 1, row["canonical_id"], time.time() - t0, n_ok, n_fail)
        logger.info("PART PILOT DONE: ok=%d fail=%d", n_ok, n_fail)

    return 0


if __name__ == "__main__":
    sys.exit(main())
