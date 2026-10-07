"""QC and technical annotations from one pass over mass spectrometry data."""

try:
    from ._version import __version__
except ImportError:  # pragma: no cover - source tree before the Hatch VCS hook runs
    from importlib.metadata import PackageNotFoundError, version

    try:
        __version__ = version("prideQC")
    except PackageNotFoundError:
        __version__ = "0+unknown"

__all__ = ["__version__"]
