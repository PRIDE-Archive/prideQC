"""Resolve the prideQC implementation from this checkout for documentation illustrations.

The figures use selected private helpers to generate *synthetic examples*. They are
maintainer/developer source only, not the public prideQC API. No source is copied or
installed by Read the Docs when building the published SVG files.

PRIDEQC_SRC can override the checkout source directory if required.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

_DEFAULT_SRC = Path(__file__).resolve().parents[2] / "src"
SRC = Path(os.environ.get("PRIDEQC_SRC", str(_DEFAULT_SRC))).expanduser().resolve()
if not (SRC / "prideqc" / "__init__.py").is_file():
    raise RuntimeError(
        f"prideQC source not found at {SRC}. Run inside a prideQC checkout "
        "or set PRIDEQC_SRC=/path/to/prideQC/src."
    )
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from prideqc import models, metrics, mass_error, mass_shift  # noqa: E402,F401
from prideqc.models import Spectrum, Precursor  # noqa: E402,F401
