#!/usr/bin/env python3
"""
Re-compute the host field + <HOST:...> header token on every example in an
existing token JSONL, using the priority-based rules in host_inference.py.

For module_library examples, derive the folder hint from the plasmid_id by
globbing Module_Library_gb/**/<stem>.gb and taking the top-level subdir.

Writes a new JSONL (leaving the original untouched) plus a host distribution
report.
"""
from __future__ import annotations
import argparse
import glob
import json
import logging
import os
import sys
from collections import Counter
from pathlib import Path
from .. import _data

sys.path.insert(0, str(Path(__file__).parent))
from host_inference import infer_host_priority, FOLDER_PRIOR

logger = logging.getLogger("rehost")


def build_folder_map(module_library_root: Path) -> dict[str, str]:
    """Map plasmid stem → top-level Module_Library subject folder."""
    out: dict[str, str] = {}
    for folder in FOLDER_PRIOR:
        folder_path = module_library_root / folder
        if not folder_path.exists():
            continue
        for gb in folder_path.rglob("*.gb"):
            out[gb.stem] = folder
    return out


def replace_host_token(tokens: list[str], new_host: str) -> list[str]:
    return [f"<HOST:{new_host}>" if t.startswith("<HOST:") else t for t in tokens]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--report", type=Path, default=None)
    ap.add_argument("--module-library-root", type=Path,
                    default=_data.data_path("Module_Library_gb"))
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s - %(levelname)s - %(message)s")

    logger.info("building folder map from %s", args.module_library_root)
    folder_map = build_folder_map(args.module_library_root)
    logger.info("folder map covers %d plasmid stems", len(folder_map))

    before = Counter()
    after = Counter()
    transitions = Counter()
    n = 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.input) as fh_in, open(args.output, "w") as fh_out:
        for line in fh_in:
            ex = json.loads(line)
            n += 1
            old_host = ex.get("host", "unknown")
            before[old_host] += 1
            folder = folder_map.get(ex["plasmid_id"]) if ex.get("source") == "module_library" else None
            new_host = infer_host_priority(ex["tokens"], folder_hint=folder, source=ex.get("source", ""))
            ex["host"] = new_host
            ex["tokens"] = replace_host_token(ex["tokens"], new_host)
            after[new_host] += 1
            transitions[(old_host, new_host)] += 1
            fh_out.write(json.dumps(ex, ensure_ascii=False) + "\n")

    logger.info("processed %d examples", n)

    report_lines = []
    def p(s: str):
        print(s); report_lines.append(s)

    p(f"# Rehost report\nTotal examples: {n}\n")
    p("## Before")
    for h, c in before.most_common(): p(f"  {h:20s} {c}  ({100*c/n:5.2f}%)")
    p("\n## After")
    for h, c in after.most_common(): p(f"  {h:20s} {c}  ({100*c/n:5.2f}%)")
    p("\n## Transitions (old → new) — top 20")
    for (a, b), c in transitions.most_common(20):
        p(f"  {a:15s} -> {b:15s} {c}")

    if args.report:
        args.report.write_text("\n".join(report_lines) + "\n")
        logger.info("wrote report → %s", args.report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
