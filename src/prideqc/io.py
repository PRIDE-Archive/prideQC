"""Atomic, standards-compliant text output."""

from __future__ import annotations

import json
import math
import os
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any, TextIO


def json_safe(value: Any) -> Any:
    """Convert NumPy scalars/arrays and undefined numbers without emitting NaN."""
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


@contextmanager
def atomic_text(path: Path) -> Iterator[TextIO]:
    """Write beside the destination, then replace it only after successful close."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            yield handle
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_json(path: Path, value: Any) -> None:
    with atomic_text(path) as handle:
        json.dump(json_safe(value), handle, indent=2, allow_nan=False, ensure_ascii=False)
        handle.write("\n")
