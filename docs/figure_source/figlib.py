"""Minimal SVG figure toolkit used by the prideQC documentation figures.

Design goals: editable SVG (real <text>, named groups, no rasterised content),
a single restrained palette, and a built-in layout linter (text/text overlaps and
connector/text crossings) so that figures can be verified programmatically.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
from PIL import ImageFont

# --------------------------------------------------------------------------- palette
INK = "#17212B"        # primary text
INK2 = "#445064"       # secondary text (>= 7:1 on white)
MUTED = "#5E6B7D"      # tertiary text (>= 5:1 on white)
FAINT = "#9AA6B5"      # graphic-only (never text)
RULE = "#CBD4DF"
GRID = "#E6EBF1"
WASH = "#F4F7FA"
PAPER = "#FFFFFF"

BLUE = "#2A6496"       # raw / MS1 / measurement
BLUE_L = "#D6E4F1"
TEAL = "#157A7A"       # QC evidence
TEAL_L = "#D3EBEA"
AMBER = "#B9650B"      # thresholds, support, caution
AMBER_L = "#F7E4CC"
PLUM = "#7A4A93"       # mass shift / PTM
PLUM_L = "#E8DCF0"
RED = "#B03A2E"        # abstain / hold
RED_L = "#F5D9D5"
GREEN = "#2F7A4F"      # accepted
GREEN_L = "#D5ECDD"
SLATE = "#4B5B70"
SLATE_L = "#DDE3EA"

SANS = "Arial, 'Helvetica Neue', Helvetica, 'Liberation Sans', 'Nimbus Sans', 'DejaVu Sans', sans-serif"
MONO = "Menlo, Consolas, 'Liberation Mono', 'Courier New', 'DejaVu Sans Mono', monospace"

def _font_dir() -> Path:
    """Locate installed metric fonts for consistent, portable layout linting.

    Override with PRIDEQC_FIGURE_FONT_DIR if Liberation is installed elsewhere.
    The font files are NOT part of this source archive or the generated assets.
    """
    import os

    override = os.environ.get("PRIDEQC_FIGURE_FONT_DIR")
    candidates = ([Path(override).expanduser()] if override else []) + [
        Path("/usr/share/fonts/truetype/liberation"),
        Path("/usr/share/fonts/truetype/liberation2"),
        Path("/usr/local/share/fonts"),
        Path.home() / ".local/share/fonts",
    ]
    for directory in candidates:
        if (directory / "LiberationSans-Regular.ttf").is_file() and (directory / "LiberationMono-Regular.ttf").is_file():
            return directory
        if directory.is_dir():
            matches = list(directory.rglob("LiberationSans-Regular.ttf"))
            for match in matches:
                if (match.parent / "LiberationMono-Regular.ttf").is_file():
                    return match.parent
    raise RuntimeError(
        "Liberation Sans and Liberation Mono metric fonts are required to lint "
        "the figure layout. Install the Liberation font package (typically "
        "fonts-liberation), or set PRIDEQC_FIGURE_FONT_DIR to its directory."
    )


_FONT_DIR: Path | None = None
_font_cache: dict[tuple[str, float], ImageFont.FreeTypeFont] = {}


def text_width(s: str, size: float, weight: str = "400", mono: bool = False, italic: bool = False) -> float:
    """Measure with Liberation (Arial/Helvetica metrics). Inter is ~3% wider; we add a margin."""
    if mono:
        name = "LiberationMono-Regular.ttf"
    else:
        bold = str(weight) in ("600", "700", "bold")
        name = ("LiberationSans-" + ("BoldItalic" if bold and italic else "Bold" if bold else "Italic" if italic else "Regular")) + ".ttf"
    key = (name, size)
    if key not in _font_cache:
        global _FONT_DIR
        if _FONT_DIR is None:
            _FONT_DIR = _font_dir()
        _font_cache[key] = ImageFont.truetype(str(_FONT_DIR / name), 100)
    f = _font_cache[key]
    return f.getlength(s) * size / 100.0 * 1.04


@dataclass
class _T:
    x0: float; y0: float; x1: float; y1: float; s: str; tag: str


@dataclass
class _Seg:
    pts: list[tuple[float, float]]; tag: str


class Fig:
    def __init__(self, w: float, h: float, title: str, desc: str):
        self.w, self.h = w, h
        self.title, self.desc = title, desc
        self.out: list[str] = []
        self.texts: list[_T] = []
        self.conns: list[_Seg] = []
        self.defs: list[str] = []
        self._gid = 0
        self._marker_ids: dict[str, str] = {}

    # ------------------------------------------------------------ primitives
    def raw(self, s: str):
        self.out.append(s)

    def g(self, id_: str | None = None, **attrs):
        a = "".join(f' {k.replace("_", "-")}="{v}"' for k, v in attrs.items())
        self.out.append(f'<g{" id=%s" % chr(34) + id_ + chr(34) if id_ else ""}{a}>')

    def end(self):
        self.out.append("</g>")

    def rect(self, x, y, w, h, fill="none", stroke="none", sw=1, rx=0, opacity=None, extra=""):
        o = f' opacity="{opacity}"' if opacity is not None else ""
        r = f' rx="{rx}"' if rx else ""
        self.out.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{r}{o}{extra}/>')

    def line(self, x1, y1, x2, y2, stroke=RULE, sw=1, dash=None, cap="butt", opacity=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        o = f' opacity="{opacity}"' if opacity is not None else ""
        self.out.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" stroke-width="{sw}" stroke-linecap="{cap}"{d}{o}/>')

    def circle(self, x, y, r, fill=INK, stroke="none", sw=1, opacity=None):
        o = f' opacity="{opacity}"' if opacity is not None else ""
        self.out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{o}/>')

    def path(self, d, fill="none", stroke=INK, sw=1.5, dash=None, opacity=None, cap="round", join="round", extra=""):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        o = f' opacity="{opacity}"' if opacity is not None else ""
        self.out.append(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" stroke-linecap="{cap}" stroke-linejoin="{join}"{da}{o}{extra}/>')

    def poly(self, pts, fill="none", stroke=INK, sw=1.5, opacity=None, close=False, dash=None):
        d = "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts) + (" Z" if close else "")
        self.path(d, fill=fill, stroke=stroke, sw=sw, opacity=opacity, dash=dash)

    def text(self, x, y, s, size=15, fill=INK, weight="400", anchor="start", mono=False, italic=False,
             tag="", opacity=None, spacing=None, check=True, baseline="alphabetic", halo=False):
        fam = MONO if mono else SANS
        st = ' font-style="italic"' if italic else ""
        o = f' opacity="{opacity}"' if opacity is not None else ""
        ls = f' letter-spacing="{spacing}"' if spacing else ""
        db = f' dominant-baseline="{baseline}"' if baseline != "alphabetic" else ""
        if halo:   # underlay copy (renderer-independent halo)
            self.out.append(
                f'<text x="{x:.2f}" y="{y:.2f}" font-family="{fam}" font-size="{size}" font-weight="{weight}" '
                f'fill="{PAPER}" stroke="{PAPER}" stroke-width="4.5" stroke-linejoin="round" text-anchor="{anchor}"{st}{ls}{db} aria-hidden="true">{escape(s)}</text>')
        hl = ""
        self.out.append(
            f'<text x="{x:.2f}" y="{y:.2f}" font-family="{fam}" font-size="{size}" font-weight="{weight}" '
            f'fill="{fill}" text-anchor="{anchor}"{st}{o}{ls}{db}{hl}>{escape(s)}</text>')
        if check and s.strip():
            w = text_width(s, size, weight, mono, italic)
            if anchor == "start": x0 = x
            elif anchor == "middle": x0 = x - w / 2
            else: x0 = x - w
            top = y - size * 0.80 if baseline == "alphabetic" else y - size * 0.5
            bot = y + size * 0.22 if baseline == "alphabetic" else y + size * 0.5
            self.texts.append(_T(x0, top, x0 + w, bot, s, tag))

    def lines(self, x, y, rows, size=15, lh=1.35, **kw):
        for i, r in enumerate(rows):
            self.text(x, y + i * size * lh, r, size=size, **kw)

    # rich text: list of (string, dict) runs on one baseline
    def rich(self, x, y, runs, size=15, anchor="start", tag="", mono=False):
        total = sum(text_width(s, size, o.get("weight", "400"), o.get("mono", mono), o.get("italic", False)) for s, o in runs)
        cx = x if anchor == "start" else x - total / 2 if anchor == "middle" else x - total
        parts = []
        for s, o in runs:
            fam = MONO if o.get("mono", mono) else SANS
            st = ' font-style="italic"' if o.get("italic") else ""
            parts.append(f'<tspan font-family="{fam}" font-weight="{o.get("weight","400")}" fill="{o.get("fill",INK)}"{st}>{escape(s)}</tspan>')
        self.out.append(f'<text x="{cx:.2f}" y="{y:.2f}" font-size="{size}" text-anchor="start" xml:space="preserve">{"".join(parts)}</text>')
        self.texts.append(_T(cx, y - size * 0.8, cx + total, y + size * 0.22, "".join(s for s, _ in runs), tag))
        return total

    # ------------------------------------------------------------ connectors
    def marker(self, color: str, size=7) -> str:
        key = f"{color}-{size}"
        if key not in self._marker_ids:
            mid = f"ah{len(self._marker_ids)}"
            self._marker_ids[key] = mid
            self.defs.append(
                f'<marker id="{mid}" viewBox="0 0 10 10" refX="8.2" refY="5" markerWidth="{size}" markerHeight="{size}" orient="auto-start-reverse">'
                f'<path d="M0.8,1.2 L9,5 L0.8,8.8 Z" fill="{color}" stroke="none"/></marker>')
        return self._marker_ids[key]

    def arrow(self, pts, color=SLATE, sw=2, head=True, dash=None, tag="", curve=False, hs=7, opacity=None):
        """Polyline (or smooth Catmull-like cubic if curve) connector with arrow head; registered for crossing checks."""
        pts = [(float(a), float(b)) for a, b in pts]
        if curve and len(pts) == 4:
            (x0, y0), (x1, y1), (x2, y2), (x3, y3) = pts
            d = f"M{x0:.2f},{y0:.2f} C{x1:.2f},{y1:.2f} {x2:.2f},{y2:.2f} {x3:.2f},{y3:.2f}"
            t = np.linspace(0, 1, 40)
            samp = [((1-u)**3*x0+3*(1-u)**2*u*x1+3*(1-u)*u**2*x2+u**3*x3,
                     (1-u)**3*y0+3*(1-u)**2*u*y1+3*(1-u)*u**2*y2+u**3*y3) for u in t]
        else:
            d = "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts)
            samp = pts
        mk = f' marker-end="url(#{self.marker(color, hs)})"' if head else ""
        da = f' stroke-dasharray="{dash}"' if dash else ""
        o = f' opacity="{opacity}"' if opacity is not None else ""
        self.out.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{sw}" stroke-linecap="round" stroke-linejoin="round"{da}{mk}{o}/>')
        self.conns.append(_Seg(samp, tag))

    # ------------------------------------------------------------ lint
    def lint(self, name: str, allow_text_pairs: set[tuple[str, str]] = frozenset()) -> list[str]:
        problems: list[str] = []
        T = self.texts
        for i in range(len(T)):
            a = T[i]
            if a.x0 < -1 or a.x1 > self.w + 1 or a.y0 < -1 or a.y1 > self.h + 1:
                problems.append(f"[{name}] text out of canvas: {a.s!r} ({a.x0:.0f},{a.y0:.0f})-({a.x1:.0f},{a.y1:.0f})")
            for j in range(i + 1, len(T)):
                b = T[j]
                pad = 1.0
                if a.x0 < b.x1 - pad and b.x0 < a.x1 - pad and a.y0 < b.y1 - pad and b.y0 < a.y1 - pad:
                    if (a.tag, b.tag) in allow_text_pairs or (b.tag, a.tag) in allow_text_pairs:
                        continue
                    problems.append(f"[{name}] text overlap: {a.s!r} x {b.s!r}")
        for seg in self.conns:
            for t in T:
                if t.tag and t.tag == seg.tag:
                    continue
                if _poly_hits_rect(seg.pts, t.x0 - 2, t.y0 - 1, t.x1 + 2, t.y1 + 1):
                    problems.append(f"[{name}] connector {seg.tag or '?'} crosses text {t.s!r}")
        for i in range(len(self.conns)):
            for j in range(i + 1, len(self.conns)):
                a, b = self.conns[i], self.conns[j]
                if _polys_cross(a.pts, b.pts):
                    problems.append(f"[{name}] connectors cross: {a.tag or i} x {b.tag or j}")
        return problems

    # ------------------------------------------------------------ output
    def svg(self) -> str:
        head = (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w:g} {self.h:g}" width="{self.w:g}" height="{self.h:g}" '
            f'role="img" aria-labelledby="t d">\n<title id="t">{escape(self.title)}</title>\n<desc id="d">{escape(self.desc)}</desc>\n'
            f'<defs>{"".join(self.defs)}</defs>\n'
            f'<rect id="background" width="{self.w:g}" height="{self.h:g}" fill="{PAPER}"/>\n')
        body = "\n".join(o for o in self.out if not re.search(r"<text[^>]*></text>", o))
        return head + body + "\n</svg>\n"

    def save(self, path: Path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.svg(), encoding="utf-8")


def _seg_rect(p, q, x0, y0, x1, y1) -> bool:
    # Liang–Barsky segment/rect clip test
    (px, py), (qx, qy) = p, q
    dx, dy = qx - px, qy - py
    t0, t1 = 0.0, 1.0
    for pp, qq in ((-dx, px - x0), (dx, x1 - px), (-dy, py - y0), (dy, y1 - py)):
        if pp == 0:
            if qq < 0: return False
        else:
            r = qq / pp
            if pp < 0:
                if r > t1: return False
                t0 = max(t0, r)
            else:
                if r < t0: return False
                t1 = min(t1, r)
    return True


def _poly_hits_rect(pts, x0, y0, x1, y1) -> bool:
    return any(_seg_rect(pts[i], pts[i + 1], x0, y0, x1, y1) for i in range(len(pts) - 1))


def _ccw(a, b, c):
    return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])


def _seg_cross(a, b, c, d) -> bool:
    return _ccw(a, c, d) != _ccw(b, c, d) and _ccw(a, b, c) != _ccw(a, b, d)


def _polys_cross(A, B) -> bool:
    for i in range(len(A) - 1):
        for j in range(len(B) - 1):
            # touching endpoints are fine (shared junctions); only proper crossings count
            if _seg_cross(A[i], A[i + 1], B[j], B[j + 1]):
                if min(math.dist(A[i], B[j]), math.dist(A[i], B[j+1]), math.dist(A[i+1], B[j]), math.dist(A[i+1], B[j+1])) < 1.5:
                    continue
                return True
    return False


# --------------------------------------------------------------------------- charts
class Scale:
    def __init__(self, d0, d1, r0, r1, log=False):
        self.d0, self.d1, self.r0, self.r1, self.log = d0, d1, r0, r1, log

    def __call__(self, v):
        if self.log:
            v, a, b = math.log10(v), math.log10(self.d0), math.log10(self.d1)
        else:
            a, b = self.d0, self.d1
        return self.r0 + (v - a) / (b - a) * (self.r1 - self.r0)


def axes(f: Fig, x, y, w, h, xs: Scale, ys: Scale, xticks, yticks, xlabel="", ylabel="",
         xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:g}", grid=True, size=12.5, left_axis=True,
         xlabel_dy=38, ylabel_dx=-48, tag="ax", ylabel_rot=True, xticklabels=True, yticklabels=True):
    """Draw light grid + left/bottom axis in a plot box whose top-left is (x, y)."""
    base = y + h
    if grid:
        for v in yticks:
            f.line(x, ys(v), x + w, ys(v), GRID, 1)
    f.line(x, base, x + w, base, INK2, 1.2)
    if left_axis:
        f.line(x, y, x, base, INK2, 1.2)
    for v in xticks:
        f.line(xs(v), base, xs(v), base + 5, INK2, 1)
        if xticklabels:
            f.text(xs(v), base + 5 + size + 3, xfmt(v), size=size, fill=INK2, anchor="middle", tag=tag)
    for v in yticks:
        f.line(x - 5, ys(v), x, ys(v), INK2, 1)
        if yticklabels:
            f.text(x - 9, ys(v) + size * 0.35, yfmt(v), size=size, fill=INK2, anchor="end", tag=tag)
    if xlabel:
        f.text(x + w / 2, base + xlabel_dy, xlabel, size=13.5, fill=INK2, anchor="middle", tag=tag)
    if ylabel:
        cx, cy = x + ylabel_dx, y + h / 2
        f.out.append(f'<text x="{cx:.2f}" y="{cy:.2f}" transform="rotate(-90 {cx:.2f} {cy:.2f})" font-family="{SANS}" font-size="13.5" fill="{INK2}" text-anchor="middle">{escape(ylabel)}</text>')


def smooth_path(pts, closed_to=None):
    d = "M" + " L".join(f"{a:.2f},{b:.2f}" for a, b in pts)
    if closed_to is not None:
        d += f" L{pts[-1][0]:.2f},{closed_to:.2f} L{pts[0][0]:.2f},{closed_to:.2f} Z"
    return d


def kde(vals, grid, bw):
    vals = np.asarray(vals)
    z = (grid[:, None] - vals[None, :]) / bw
    return np.exp(-0.5 * z * z).sum(axis=1) / (len(vals) * bw * math.sqrt(2 * math.pi))


def panel_label(f: Fig, x, y, letter, title, sub=None, color=INK):
    """Editorial panel header: bold letter, title, thin rule underneath (no boxes)."""
    f.text(x, y, letter, size=21, weight="700", fill=color, tag="ph")
    f.text(x + 28, y, title, size=17, weight="600", fill=INK, tag="ph")
    if sub:
        f.text(x + 28, y + 21, sub, size=13.5, fill=MUTED, tag="ph")


def chip(f: Fig, x, y, s, fg=INK, bg=WASH, size=13.5, weight="600", pad=9, h=24, mono=False, stroke="none", anchor="start", tag="chip"):
    """Small label with a tinted background (square-cornered tag, not a rounded flowchart box)."""
    w = text_width(s, size, weight, mono) + 2 * pad
    x0 = x if anchor == "start" else x - w / 2 if anchor == "middle" else x - w
    f.rect(x0, y - h / 2, w, h, fill=bg, stroke=stroke, sw=1)
    f.text(x0 + w / 2, y + size * 0.35, s, size=size, weight=weight, fill=fg, anchor="middle", mono=mono, tag=tag)
    return x0, w
