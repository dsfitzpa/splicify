#!/usr/bin/env python3
"""
Precompute per-canonical-part statistics from the tokenized corpus.

For every (role, payload) pair that appears as a feature token, compute:
  - occurrence_count
  - host_distribution (share by mammalian/bacterial/insect/plant/yeast/unknown)
  - top-K co-occurring feature tokens (same plasmid, not self)
  - top-K enclosing module_types (MOD_OPEN tokens whose span contains this part)
  - source_distribution (module_library / ncbi_engineered / etc.)
  - example plasmid IDs (first N)

Output: JSONL, one row per canonical part. Feeds the part-annotation prompt
builder (PLASMID_LLM_ANNOTATION_PLAN.md §4.2, §5.2).
"""
from __future__ import annotations
import argparse
import json
import logging
import re
from collections import Counter, defaultdict
from pathlib import Path

logger = logging.getLogger("cooccur")

# Role prefixes that count as "feature parts" (i.e. tokens we want to annotate)
FEATURE_ROLE_PREFIXES = {
    "CDS", "PROMOTER", "PROMOTER_POL2", "PROMOTER_POL3", "POLYA", "TERMINATOR",
    "ORI", "ORIT", "LTR", "ITR", "RBS", "KOZAK", "OPERATOR", "ENHANCER",
    "INSULATOR", "RECOMB", "PROTEIN_BIND", "PRIMER_BIND", "SIG_PEP", "MAT_PEP",
    "INTRON", "EXON", "NCRNA", "TRNA", "RRNA", "REG", "MISC", "MCS",
    "NLS", "TAG", "LINKER", "LINKER_2A", "GAP", "MARKER",
}

# Noise tokens to skip when building co-occurrence lists
SKIP_PREFIXES_FOR_NEIGHBORS = {
    "BOS", "EOS", "PAD", "UNK", "TOPOLOGY", "LEN_BIN", "HOST", "SOURCE",
    "ROTATION_IDX", "MOD_OPEN", "MOD_CLOSE", "UPSTREAM_REG", "/UPSTREAM_REG",
    "CDS_MOD", "/CDS_MOD", "DOWNSTREAM_REG", "/DOWNSTREAM_REG",
    "INT", "CLN",
    "VB_CAS", "/VB_CAS", "VB_FAMILY", "VB_SYS", "DRIVES", "BICISTRONIC",
}

TOKEN_RE = re.compile(r"<(?P<role>[^:>]+)(?::(?P<payload>[^>]*))?>")


def parse_token(t: str) -> tuple[str, str | None]:
    m = TOKEN_RE.match(t)
    if not m:
        return t, None
    role = m.group("role")
    payload = m.group("payload")
    return role, payload


def canonical_id(role: str, payload: str | None) -> str:
    """Canonical key for a part. Payload is normalized (strip trailing -NNN numeric suffixes)."""
    if payload is None:
        return f"{role}"
    # Normalize trailing suffixes like -009, _(2), -001
    p = re.sub(r"[-_]\(?\d{1,4}\)?$", "", payload)
    p = re.sub(r"-\d{3}$", "", p)
    p = p.strip().replace(" ", "_")
    return f"{role}:{p}"


def walk_tokens_with_module_context(tokens: list[str]):
    """Yield (role, payload, enclosing_module_type_or_None) for each feature token."""
    stack: list[str] = []
    for t in tokens:
        role, payload = parse_token(t)
        if role == "MOD_OPEN":
            stack.append(payload or "unknown_module")
            continue
        if role == "MOD_CLOSE":
            if stack:
                stack.pop()
            continue
        if role in SKIP_PREFIXES_FOR_NEIGHBORS:
            continue
        if role in FEATURE_ROLE_PREFIXES:
            yield role, payload, (stack[-1] if stack else None)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--min-occurrences", type=int, default=2,
                    help="Drop canonical_ids seen in fewer than N plasmids")
    ap.add_argument("--top-neighbors", type=int, default=20)
    ap.add_argument("--top-modules", type=int, default=5)
    ap.add_argument("--example-plasmids", type=int, default=5)
    ap.add_argument("--rotation", type=int, default=0)
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s - %(levelname)s - %(message)s")

    # Per-canonical_id aggregators
    counts: Counter = Counter()
    hosts: dict[str, Counter] = defaultdict(Counter)
    sources: dict[str, Counter] = defaultdict(Counter)
    neighbors: dict[str, Counter] = defaultdict(Counter)
    modules: dict[str, Counter] = defaultdict(Counter)
    examples: dict[str, list[str]] = defaultdict(list)
    roles_seen: dict[str, Counter] = defaultdict(Counter)

    n_plasmids = 0
    with open(args.corpus) as fh:
        for line in fh:
            ex = json.loads(line)
            if ex.get("rotation_idx", 0) != args.rotation:
                continue
            n_plasmids += 1
            pid = ex["plasmid_id"]
            host = ex.get("host", "unknown")
            src = ex.get("source", "unknown")
            # First pass: collect canonical_ids and enclosing module per feature token
            parts_here: list[tuple[str, str | None, str | None]] = list(
                walk_tokens_with_module_context(ex["tokens"])
            )
            if not parts_here:
                continue
            canon_here = [(canonical_id(r, p), r, p, m) for r, p, m in parts_here]
            # Count per-plasmid uniqueness (don't over-count a canonical_id that appears twice in one plasmid)
            seen_in_this_plasmid: set[str] = set()
            for cid, role, payload, enc in canon_here:
                if cid in seen_in_this_plasmid:
                    continue
                seen_in_this_plasmid.add(cid)
                counts[cid] += 1
                hosts[cid][host] += 1
                sources[cid][src] += 1
                roles_seen[cid][role] += 1
                if enc:
                    modules[cid][enc] += 1
                if len(examples[cid]) < args.example_plasmids:
                    examples[cid].append(pid)
            # Neighbor co-occurrence — per-plasmid, excluding self
            unique_cids = list(seen_in_this_plasmid)
            for cid in unique_cids:
                for other in unique_cids:
                    if other != cid:
                        neighbors[cid][other] += 1

    logger.info("scanned %d plasmids; %d unique canonical parts", n_plasmids, len(counts))

    # Write output
    kept = 0
    dropped = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        for cid, n in counts.most_common():
            if n < args.min_occurrences:
                dropped += 1
                continue
            # Dominant role for this canonical_id (ties broken by count)
            role = roles_seen[cid].most_common(1)[0][0]
            # Payload (extract from cid)
            payload = cid.split(":", 1)[1] if ":" in cid else None
            row = {
                "canonical_id": cid,
                "role": role,
                "payload": payload,
                "occurrence_count": n,
                "corpus_coverage": round(n / max(1, n_plasmids), 4),
                "host_distribution": dict(hosts[cid].most_common()),
                "source_distribution": dict(sources[cid].most_common()),
                "top_neighbors": [
                    {"canonical_id": k, "count": v}
                    for k, v in neighbors[cid].most_common(args.top_neighbors)
                ],
                "top_enclosing_modules": [
                    {"module_type": k, "count": v}
                    for k, v in modules[cid].most_common(args.top_modules)
                ],
                "example_plasmids": examples[cid],
            }
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
            kept += 1
    logger.info("wrote %d canonical parts → %s (dropped %d below --min-occurrences=%d)",
                kept, args.output, dropped, args.min_occurrences)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
