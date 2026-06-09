#!/usr/bin/env bash
# Downloads Splicify reference data from Zenodo into $SPLICIFY_DATA_DIR.
#
# Usage:
#   SPLICIFY_DATA_DIR=~/.splicify/data ./scripts/fetch_data.sh [artefact_name]
#
# Without an artefact_name, downloads every artefact in data/MANIFEST.yml.
# With one, downloads only that artefact. Skips artefacts already present.
#
# Requires: python3, pyyaml (pip install pyyaml).
# Verifies SHA256 against the manifest after download.

set -euo pipefail

DEST="${SPLICIFY_DATA_DIR:-$HOME/.splicify/data}"
MANIFEST="$(cd "$(dirname "$0")/.." && pwd)/data/MANIFEST.yml"

mkdir -p "$DEST"

python3 - "$MANIFEST" "$DEST" "${1:-}" <<'PY'
import hashlib
import os
import pathlib
import sys
import tarfile
import urllib.request

import yaml

manifest_path, dest_str, want = sys.argv[1], sys.argv[2], sys.argv[3] or None
dest = pathlib.Path(dest_str)
dest.mkdir(parents=True, exist_ok=True)
manifest = yaml.safe_load(pathlib.Path(manifest_path).read_text())

for name, meta in manifest["artefacts"].items():
    if want and name != want:
        continue
    target = dest / meta["target_dir"]
    if target.exists() and any(target.iterdir()):
        print(f"[skip] {name} already present at {target}")
        continue

    doi_suffix = meta["zenodo_doi"].split(".")[-1]
    url = f"https://zenodo.org/record/{doi_suffix}/files/{name}.tar.gz"
    tmp = dest / f"{name}.tar.gz"
    print(f"[fetch] {name} <- {url}")
    urllib.request.urlretrieve(url, tmp)

    actual_sha = hashlib.sha256(tmp.read_bytes()).hexdigest()
    expected = meta["sha256"]
    if expected == "TBD" or expected == "TBD-on-zenodo-upload":
        print(f"[warn] {name}: SHA placeholder; skipping verification")
    elif actual_sha != expected:
        tmp.unlink()
        raise SystemExit(
            f"[error] SHA256 mismatch for {name}: got {actual_sha[:16]}..., "
            f"expected {expected[:16]}..."
        )
    else:
        print(f"[ok] {name}: SHA verified")

    print(f"[extract] {name} -> {target}")
    with tarfile.open(tmp, "r:gz") as tar:
        tar.extractall(dest)
    tmp.unlink()

print("[done]")
PY
