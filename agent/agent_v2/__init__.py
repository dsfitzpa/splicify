"""splicify-agent — hierarchical Sonnet 4.6 plasmid + CRISPR cloning agent.

After PR #1 (splicify-core extraction), this package depends on the installed
``splicify`` distribution. Sys.path mutations have been removed; install the
core package into the same environment::

    pip install -e ../backend       # from a monorepo checkout
    # or, once PyPI publishing lands:
    pip install splicify

Failure to install ``splicify`` raises a clean ImportError at module load
time rather than silently importing a stale/missing module path.
"""

try:
    import splicify_api  # noqa: F401
except ImportError as exc:  # pragma: no cover — install-time guidance
    raise ImportError(
        "splicify-agent requires the 'splicify' core package. "
        "From a monorepo checkout: 'pip install -e ../backend'. "
        "From PyPI (once published): 'pip install splicify'."
    ) from exc
