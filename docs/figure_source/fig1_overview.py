"""Figure 1 - prideQC end-to-end overview."""
from __future__ import annotations
import math
import numpy as np
from figlib import *

rng = np.random.default_rng(3)


def doc(f, x, y, w=30, h=38, fill=PAPER, stroke=SLATE, sw=1.6):
    c = 9
    f.poly([(x, y), (x + w - c, y), (x + w, y + c), (x + w, y + h), (x, y + h)], fill=fill, stroke=stroke, sw=sw, close=True)
    f.poly([(x + w - c, y), (x + w - c, y + c), (x + w, y + c)], stroke=stroke, sw=sw * 0.8)


def stage_head(f, x, y, num, title, sub, color):
    f.rect(x, y - 30, 34, 4, fill=color)
    f.text(x, y, f"{num}", size=15, weight="700", fill=color, tag="sh")
    f.text(x + 20, y, title, size=17, weight="700", tag="sh")
    for i, s in enumerate(sub):
        f.text(x, y + 24 + i * 18, s, size=13.2, fill=MUTED, tag="sh")


# ------------------------------------------------------------------ engine glyphs (62 x 44 box at x,y top-left)
def g_tic(f, x, y):
    t = np.linspace(0, 1, 40)
    v = 0.15 + 0.8 * (1 / (1 + np.exp(-(t - 0.18) * 18))) * (1 / (1 + np.exp((t - 0.86) * 18))) * (1 + 0.18 * np.sin(t * 22))
    pts = [(x + 62 * a, y + 44 - 40 * b) for a, b in zip(t, v)]
    f.path(smooth_path(pts, closed_to=y + 44), fill=BLUE_L, stroke="none")
    f.poly(pts, stroke=BLUE, sw=2)
    f.line(x, y + 44, x + 62, y + 44, INK2, 1.2)


def g_raster(f, x, y):
    base = y + 44
    for c in range(6):
        cx = x + 3 + c * 10.5
        f.rect(cx, base - 40, 2.6, 40, fill=BLUE)
        for k in range(4):
            h = 10 + 12 * rng.random()
            f.rect(cx + 3 + k * 1.6, base - h, 1.3, h, fill=[TEAL, TEAL, AMBER, TEAL][k])
    f.line(x, base, x + 62, base, INK2, 1.2)


def g_precision(f, x, y):
    g = np.linspace(-1, 1, 40)
    pts = [(x + 31 + 29 * a, y + 42 - 36 * math.exp(-0.5 * (a / 0.22) ** 2)) for a in g]
    f.rect(x + 31 - 29 * 0.8, y + 4, 29 * 1.6, 38, fill=GREEN_L, opacity=0.8)
    f.path(smooth_path(pts, closed_to=y + 42), fill=BLUE_L, stroke="none")
    f.poly(pts, stroke=BLUE, sw=2)
    f.line(x, y + 42, x + 62, y + 42, INK2, 1.2)
    f.line(x + 31 - 29 * 0.8, y + 4, x + 31 - 29 * 0.8, y + 42, GREEN, 1.6); f.line(x + 31 + 29 * 0.8, y + 4, x + 31 + 29 * 0.8, y + 42, GREEN, 1.6)


def g_diag(f, x, y):
    base = y + 44
    xs = np.sort(rng.uniform(2, 60, 16))
    for i, a in enumerate(xs):
        f.line(x + a, base, x + a, base - 6 - 20 * rng.random(), FAINT, 1.6)
    for a, h in ((10, 38), (17, 30), (24, 34), (31, 28)):
        f.line(x + a, base, x + a, base - h, PLUM, 2.6)
    f.line(x, base, x + 62, base, INK2, 1.2)


def g_shift(f, x, y):
    base = y + 44
    for a in np.linspace(2, 60, 60):
        f.line(x + a, base, x + a, base - 2 - 3 * rng.random(), FAINT, 1.2)
    for a, h in ((13, 17), (28, 24), (46, 40)):
        f.line(x + a, base, x + a, base - h, PLUM, 3)
    f.line(x, base, x + 62, base, INK2, 1.2)


def build() -> Fig:
    W, H = 1400, 730
    f = Fig(W, H, "prideQC overview: raw spectra to QC evidence to refined SDRF",
            "Raw mzML or vendor RAW files are read once per file by a streaming QC engine that computes scan, acquisition, mass-error precision, diagnostic and "
            "recurrent mass-shift evidence. Structured outputs feed pmultiqc reporting and a confidence-gated cohort synthesis that refines and validates SDRF metadata.")
    f.text(60, 54, "prideQC: from raw spectra to QC evidence to refined metadata", size=28, weight="700", tag="title")
    f.text(60, 82, "Quality-control metrics and technical evidence are computed directly from the raw data; sufficiently supported evidence can refine SDRF metadata.",
           size=15, fill=MUTED, tag="title")

    # =============================================================== 1 INPUT
    stage_head(f, 60, 140, "1", "Input", ["raw or converted spectra", "from local files or PRIDE"], BLUE)
    doc(f, 62, 235, fill=BLUE_L, stroke=BLUE)
    f.text(104, 252, "mzML", size=15, weight="700", mono=True, tag="in")
    f.text(104, 270, "open standard", size=12.8, fill=MUTED, tag="in")
    doc(f, 62, 325, fill=BLUE_L, stroke=BLUE)
    f.text(104, 342, "vendor RAW", size=15, weight="700", tag="in")
    f.text(104, 360, "supported formats", size=12.8, fill=MUTED, tag="in")
    # spectrum-stream hint under inputs
    base = 470
    for i, a in enumerate(np.linspace(66, 190, 30)):
        h = 6 + 28 * abs(math.sin(i * 1.7)) * rng.uniform(0.4, 1)
        f.line(a, base, a, base - h, BLUE if i % 5 == 0 else FAINT, 2)
    f.line(62, base, 196, base, INK2, 1.2)
    f.text(62, base + 24, "scan after scan", size=13, fill=MUTED, tag="in")

    # =============================================================== 2 ENGINE
    stage_head(f, 262, 140, "2", "Streaming QC engine", ["each spectrum is read once and", "dispatched to evidence collectors"], TEAL)
    rows = [("Scan and chromatography", "TIC, base peak, RT, scan rates", g_tic),
            ("Acquisition", "cycles, charge, isolation, targets", g_raster),
            ("Mass-error precision", "repeat observations, no IDs", g_precision),
            ("Diagnostic evidence", "reporter and glycan oxonium ions", g_diag),
            ("Recurrent mass shifts", "related-spectrum mass differences", g_shift)]
    ry = [240, 312, 384, 456, 528]
    bus_l, bus_r = 262, 596
    f.line(bus_l, 212, bus_l, 556, TEAL, 3)
    f.line(bus_r, 212, bus_r, 556, TEAL, 3)
    for (t, s, g), y in zip(rows, ry):
        f.line(bus_l, y, 280, y, TEAL, 1.6)
        g(f, 286, y - 22)
        f.text(364, y - 2, t, size=15, weight="700", tag="eng")
        f.text(364, y + 17, s, size=12.8, fill=MUTED, tag="eng")
    f.text(262, 600, "Built on OpenMS readers and calculations.", size=13, fill=MUTED, tag="eng")
    # input -> engine
    f.arrow([(206, 384), (254, 384)], color=SLATE, sw=2.4, tag="a_in")

    # =============================================================== 3 OUTPUTS
    stage_head(f, 640, 140, "3", "Structured outputs", ["one set per input file"], SLATE)
    outs = [(270, "mzQC", "per-file QC document"), (395, "TSV tables", "metrics and evidence"), (520, "summary.json", "evidence and provenance")]
    for y, t, s in outs:
        f.arrow([(bus_r + 2, y), (638, y)], color=TEAL, sw=2, tag=f"ao{y}", hs=6)
        doc(f, 646, y - 19, 26, 34, fill=WASH)
        f.text(682, y - 2, t, size=14.5, weight="700", mono=(t == "summary.json"), tag="out")
        f.text(682, y + 17, s, size=12.8, fill=MUTED, tag="out")

    # =============================================================== 4 REPORTING (top track)
    stage_head(f, 856, 140, "4", "Reporting", ["optional: pmultiqc / MultiQC aggregate", "many runs into one interactive report"], BLUE)
    f.arrow([(838, 270), (870, 270)], color=BLUE, sw=2.4, tag="a_r1", hs=6) if False else None
    # collector bar
    f.arrow([(838 + 0, 270), (864, 270)], color=BLUE, sw=2.2, tag="a_r1", hs=6)
    f.rect(866, 238, 10, 64, fill=BLUE)
    f.text(871, 322, "pmultiqc /", size=13.2, weight="700", fill=BLUE, anchor="middle", tag="rep")
    f.text(871, 339, "MultiQC", size=13.2, weight="700", fill=BLUE, anchor="middle", tag="rep")
    f.arrow([(884, 270), (936, 270)], color=BLUE, sw=2.2, tag="a_r2", hs=6)
    # mini report window
    bx, by, bw, bh = 944, 196, 396, 150
    f.rect(bx, by, bw, bh, fill=PAPER, stroke=RULE, sw=1.6)
    f.rect(bx, by, bw, 20, fill=WASH, stroke=RULE, sw=1.6)
    for i in range(3):
        f.circle(bx + 12 + i * 12, by + 10, 3.2, fill=FAINT)
    yb = by + 132
    px = bx + 22
    for i in range(9):
        v = 58 if i != 5 else 20
        f.rect(px + i * 11, yb - v, 8, v, fill=AMBER if i == 5 else BLUE)
    f.line(px - 4, yb, px + 100, yb, INK2, 1)
    px2 = bx + 152
    f.line(px2, yb, px2 + 100, yb, INK2, 1)
    f.rect(px2 + 12, yb - 62, 28, 30, fill=BLUE_L, stroke=BLUE, sw=1.4)
    f.line(px2 + 12, yb - 47, px2 + 40, yb - 47, BLUE, 2)
    for k, (dx, dy) in enumerate([(62, 40), (70, 52), (78, 46), (66, 62), (74, 36), (82, 58), (60, 50), (86, 44)]):
        f.circle(px2 + dx, yb - dy, 3, fill=TEAL)
    f.circle(px2 + 72, yb - 94, 3.2, fill=AMBER)
    px3 = bx + 282
    f.line(px3, yb, px3 + 100, yb, INK2, 1); f.line(px3, yb, px3, yb - 90, INK2, 1)
    for dx, dy in [(22, 30), (30, 40), (26, 24), (36, 34), (20, 44), (32, 20), (40, 28), (28, 52)]:
        f.circle(px3 + dx, yb - dy, 3.2, fill=TEAL)
    for dx, dy in [(76, 66), (84, 74), (70, 78), (88, 62)]:
        f.circle(px3 + dx, yb - dy, 3.2, fill=AMBER)
    f.text(bx + bw / 2, by + bh + 26, "Interactive cross-run QC report", size=14.5, weight="700", fill=BLUE, anchor="middle", tag="rep")

    # =============================================================== 5 REFINEMENT (bottom track)
    f.line(840, 392, W - 60, 392, RULE, 1)
    stage_head(f, 856, 432, "5", "Cohort synthesis and SDRF refinement", [], PLUM)
    cells = [908, 1006, 1104, 1202, 1298]
    cy = 530
    f.arrow([(838, 520), (868, 520)], color=PLUM, sw=2.4, tag="a_f0", hs=6)
    # glyph 1: cohort synthesis (two clusters)
    x = cells[0]
    for k in range(9): f.circle(x - 16 + rng.normal(0, 7), cy + rng.normal(0, 7), 3.6, fill=TEAL)
    for k in range(4): f.circle(x + 20 + rng.normal(0, 4), cy - 8 + rng.normal(0, 5), 3.6, fill=AMBER)
    # glyph 2: gates (three posts + passing dots)
    x = cells[1]
    for dx in (-18, 0, 18):
        f.rect(x + dx - 2, cy - 26, 4, 52, fill=PLUM)
    f.circle(x - 30, cy, 4, fill=TEAL); f.circle(x + 30, cy, 4, fill=GREEN)
    # glyph 3: refined SDRF table
    x = cells[2]
    for i in range(4): f.line(x - 30, cy - 24 + i * 16, x + 30, cy - 24 + i * 16, INK2, 1.2)
    for i in range(4): f.line(x - 30 + i * 20, cy - 24, x - 30 + i * 20, cy + 24, INK2, 1.2)
    f.rect(x - 9, cy - 7, 18, 14, fill=GREEN_L); f.rect(x + 11, cy + 9, 18, 14, fill=PLUM_L)
    # glyph 4: validation
    x = cells[3]
    f.circle(x, cy, 25, fill=GREEN)
    f.path(f"M{x-11:.1f},{cy+1:.1f} l8,9 l15,-18", stroke=PAPER, sw=4.5)
    # glyph 5: change record
    x = cells[4]
    doc(f, x - 18, cy - 28, 36, 52, fill=WASH)
    for i in range(4): f.line(x - 10, cy - 10 + i * 9, x + 10, cy - 10 + i * 9, SLATE, 1.6)
    labels = [("Cohort", "synthesis"), ("Confidence and", "scope gates"), ("Refined", "SDRF"), ("Validation", ""), ("Change and", "audit log")]
    for x, (a, b) in zip(cells, labels):
        f.text(x, cy + 52, a, size=14, weight="700", anchor="middle", tag="ref")
        if b: f.text(x, cy + 70, b, size=14, weight="700", anchor="middle", tag="ref")
    for i in range(4):
        f.arrow([(cells[i] + 36, cy), (cells[i + 1] - 36, cy)], color=PLUM, sw=2.2, tag=f"a_f{i+1}", hs=6)
    # abstain exit above the gates glyph
    f.arrow([(cells[1], cy - 34), (cells[1], cy - 62)], color=RED, sw=2, dash="5 4", tag="a_abs", hs=6)
    f.text(cells[1] + 14, cy - 66, "abstain or hold:", size=13.2, weight="700", fill=RED, tag="abs")
    f.text(cells[1] + 14, cy - 50, "value left unchanged", size=13.2, fill=RED, tag="abs")
    f.text(856, cy + 108, "Evidence becomes metadata only after group-level confidence and scope gates;", size=13.5, fill=INK2, tag="cap")
    f.text(856, cy + 127, "every proposal and every skipped candidate is recorded.", size=13.5, fill=INK2, tag="cap")

    f.line(60, 690, W - 60, 690, RULE, 1)
    f.text(60, 716, "Reporting and refinement both start from the same per-file evidence. pmultiqc displays results; prideQC computes them.", size=14, weight="600", tag="foot")
    return f


if __name__ == "__main__":
    # Use the common build driver for portable output paths and layout linting.
    from build_all import main
    raise SystemExit(main(["--figures", "1"]))
