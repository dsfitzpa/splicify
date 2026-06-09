#!/usr/bin/env python3
"""
Respectful VectorBuilder GenBank fetcher.

Pulls the per-vector GenBank file from en.vectorbuilder.com/vector/{id}.gb for every
`representative_vectors[*].vector_id` in vector_systems.json. Cached to disk; skips
already-downloaded files.

ToS posture: this hits the same public per-vector page a user would browse to manually,
retrieving the GenBank format the UI exposes via its download menu. 2-second delay
between requests matches the existing scraper.py rate. Not a bulk-scrape of protected
endpoints. If VectorBuilder objects, delete the cache and stop.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import requests

BASE_URL = "https://en.vectorbuilder.com/vector/{vid}.gb"
USER_AGENT = "Mozilla/5.0 (compatible; plasmid-design-research/1.0)"
DELAY_SECONDS = 2.0
TIMEOUT = 30

logger = logging.getLogger("vectorbuilder_fetch")


def iter_vector_ids(systems_json: Path):
    with open(systems_json) as fh:
        data = json.load(fh)
    seen = set()
    for sys in data.get("vector_systems", []):
        for rv in sys.get("representative_vectors", []):
            vid = rv.get("vector_id") or rv.get("name")
            if not vid or vid in seen:
                continue
            seen.add(vid)
            yield vid, rv.get("description", ""), sys.get("id", "")


def fetch_one(vid: str, out_path: Path, session: requests.Session) -> bool:
    if out_path.exists() and out_path.stat().st_size > 100:
        logger.debug("cache hit: %s", vid)
        return True
    url = BASE_URL.format(vid=vid)
    try:
        resp = session.get(url, timeout=TIMEOUT, headers={"User-Agent": USER_AGENT})
    except requests.RequestException as exc:
        logger.warning("request failed for %s: %s", vid, exc)
        return False
    if resp.status_code != 200:
        logger.warning("HTTP %d for %s", resp.status_code, vid)
        return False
    body = resp.text
    if not body.lstrip().startswith("LOCUS"):
        logger.warning("non-GenBank body for %s (first 80 chars: %r)", vid, body[:80])
        return False
    out_path.write_text(body)
    logger.info("fetched %s → %s (%d bytes)", vid, out_path.name, len(body))
    return True


def fetch_all(systems_json: Path, out_dir: Path) -> dict[str, bool]:
    out_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, bool] = {}
    session = requests.Session()
    ids = list(iter_vector_ids(systems_json))
    logger.info("fetching %d VectorBuilder GenBanks into %s", len(ids), out_dir)
    first = True
    for vid, desc, sys_id in ids:
        if not first:
            time.sleep(DELAY_SECONDS)
        first = False
        out_path = out_dir / f"{vid}.gb"
        ok = fetch_one(vid, out_path, session)
        results[vid] = ok
    # Also emit a manifest so the tokenizer can pair each id back to its shorthand description
    manifest_path = out_dir / "manifest.json"
    manifest = {
        vid: {"description": desc, "system_id": sys_id, "fetched": ok}
        for (vid, desc, sys_id), ok in zip(ids, (results[v] for v, _, _ in ids))
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))
    logger.info("wrote manifest → %s", manifest_path)
    return results


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Fetch VectorBuilder representative GenBanks")
    ap.add_argument("--systems-json", type=Path,
                    default=Path("backend/splicify_api/vectorbuilder_db/vector_systems.json"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("backend/splicify_api/vectorbuilder_db/genbank_files"))
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args(argv)
    logging.basicConfig(level=args.log_level,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    results = fetch_all(args.systems_json, args.out_dir)
    n_ok = sum(1 for v in results.values() if v)
    logger.info("DONE: %d/%d succeeded", n_ok, len(results))
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
