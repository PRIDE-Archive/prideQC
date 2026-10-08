"""Optional developer-only SVG to PNG preview rendering."""

from __future__ import annotations

import re
from pathlib import Path


def render(svg: Path, png: Path, scale: float = 1.0) -> None:
    import cairosvg

    source = Path(svg).read_text(encoding="utf-8")
    # Use installed Liberation fonts to match the generator's layout metric.
    source = re.sub(r'font-family="Menlo[^"]*"', 'font-family="Liberation Mono"', source)
    source = re.sub(r'font-family="Arial[^"]*"', 'font-family="Liberation Sans"', source)
    png = Path(png)
    png.parent.mkdir(parents=True, exist_ok=True)
    cairosvg.svg2png(
        bytestring=source.encode("utf-8"), write_to=str(png), scale=scale, background_color="white"
    )
