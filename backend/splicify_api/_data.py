"""Resolver for reference databases + curated libraries.

Replaces hardcoded ``Path(__file__).parent / "feature_db_data"`` and the
``.parent.parent / "Module_Library_gb"`` variants scattered across
``splicify_api``. After PR #1, those become::

    from splicify_api import _data
    ...
    feature_db = _data.data_path("feature_db_data")
    module_lib = _data.data_path("Module_Library_gb")

Resolution order:

1. ``$SPLICIFY_DATA_DIR`` if set and exists.
2. ``~/.splicify/data`` if exists (typical user install via fetch_data.sh).
3. ``<package>/feature_db_data/_test_fixtures`` (tiny smoke-test set; ships
   in the wheel so ``pytest`` works on a fresh checkout).

Test fixtures are intentionally last so production deployments without
env-var configuration don't silently fall back to a stub dataset.
"""
from __future__ import annotations

import os
import pathlib
from typing import Iterator

_PACKAGE_ROOT = pathlib.Path(__file__).resolve().parent
_TEST_FIXTURES = _PACKAGE_ROOT / "feature_db_data" / "_test_fixtures"


def _candidates() -> Iterator[pathlib.Path]:
    env = os.environ.get("SPLICIFY_DATA_DIR")
    if env:
        yield pathlib.Path(env).expanduser().resolve()
    yield pathlib.Path("~/.splicify/data").expanduser()
    yield _TEST_FIXTURES


def data_root() -> pathlib.Path:
    """Return the active data root, raising if no candidate exists."""
    for candidate in _candidates():
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        "No Splicify data directory found. Set $SPLICIFY_DATA_DIR or run "
        "scripts/fetch_data.sh to download the reference DBs from Zenodo."
    )


def data_path(*parts: str) -> pathlib.Path:
    """Resolve a path inside the active data root.

    Raises FileNotFoundError if the resolved path does not exist — fail
    loudly at lookup time rather than silently returning a non-existent
    Path that downstream callers will misuse.

    Example::

        data_path("feature_db_data", "feature_reference.fna")
        data_path("Module_Library_gb")
    """
    resolved = data_root().joinpath(*parts)
    if not resolved.exists():
        raise FileNotFoundError(
            f"{resolved} not found. Check $SPLICIFY_DATA_DIR or re-run "
            "scripts/fetch_data.sh"
        )
    return resolved
