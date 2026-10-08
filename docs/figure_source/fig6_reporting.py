"""Figure 6 - pmultiqc / MultiQC reporting layer (synthetic example)."""
from __future__ import annotations
import numpy as np
from figlib import *

rng = np.random.default_rng(14)
CHC = {"1+": "#8A97A8", "2+": TEAL, "3+": BLUE, "4+": AMBER, "≥5": PLUM}


def doc(f, x, y, w=26, h=32, fill=PAPER, stroke=SLATE):
    c = 8
    f.poly([(x, y), (x + w - c, y), (x + w, y + c), (x + w, y + h), (x, y + h)], fill=fill, stroke=stroke, sw=1.4, close=True)
    f.poly([(x + w - c, y), (x + w - c, y + c), (x + w, y + c)], stroke=stroke, sw=1.2)


def build() -> Fig:
    W, H = 1400, 830
    f = Fig(W, H, "Cross-run reporting with pmultiqc and MultiQC",
            "Many per-run mzQC files from prideQC are aggregated by pmultiqc and MultiQC into one interactive HTML report with run-to-run distributions, "
            "outlier runs, acquisition characteristics and metric comparisons. The reporting layer displays prideQC values and does not recompute them.")
    f.text(60, 54, "Reporting across runs", size=28, weight="700", tag="title")
    f.text(60, 82, "pmultiqc and MultiQC aggregate and display per-run mzQC files. The scientific calculations were already done by prideQC.", size=15, fill=MUTED, tag="title")

    # ---------------------------------------------------------- left: files -> collector
    f.text(60, 150, "prideQC output", size=17, weight="700", tag="L")
    f.text(60, 171, "one mzQC per input file", size=13.5, fill=MUTED, tag="L")
    names = ["run01.mzQC", "run02.mzQC", "run03.mzQC", "run04.mzQC", "run05.mzQC", "run06.mzQC", "run24.mzQC"]
    y0 = 214
    bar_x = 330
    for i, nme in enumerate(names):
        y = y0 + i * 50
        if i == 6:
            [f.circle(73, y - 30 + k * 7, 1.8, fill=MUTED) for k in range(3)]
        doc(f, 62, y - 16, 22, 28, fill=WASH)
        f.text(96, y + 4, nme, size=14, mono=True, tag="L")
        f.arrow([(222, y - 2), (bar_x - 6, y - 2)], color=SLATE, sw=1.8, tag=f"fa{i}", hs=6)
    yb0, yb1 = y0 - 16, y0 + 6 * 50 + 12
    f.rect(bar_x, yb0, 14, yb1 - yb0, fill=TEAL)
    f.text(bar_x + 7, yb1 + 30, "pmultiqc / MultiQC", size=16, weight="700", fill=TEAL, anchor="middle", tag="C")
    f.text(bar_x + 7, yb1 + 52, "generic mzQC module", size=13.5, fill=MUTED, anchor="middle", tag="C")
    f.text(60, yb1 + 100, "multiqc --module mzqc results/", size=14, mono=True, tag="C")
    f.text(60, yb1 + 124, "writes multiqc_report.html and multiqc_data/", size=13.5, fill=MUTED, tag="C")
    f.arrow([(bar_x + 22, (yb0 + yb1) / 2), (520, (yb0 + yb1) / 2)], color=TEAL, sw=3, tag="out", hs=7)

    # ---------------------------------------------------------- right: report window
    bx, by, bw, bh = 540, 120, 800, 620
    f.rect(bx, by, bw, bh, fill=PAPER, stroke=RULE, sw=1.6)
    f.rect(bx, by, bw, 34, fill=WASH, stroke=RULE, sw=1.6)
    for i in range(3):
        f.circle(bx + 20 + i * 18, by + 17, 5, fill=FAINT)
    f.text(bx + bw / 2, by + 22, "multiqc_report.html", size=13.5, fill=MUTED, anchor="middle", mono=True, tag="bar")
    # sidebar
    sx = bx + 18
    f.text(sx, by + 66, "Sections", size=13.5, weight="700", fill=INK2, tag="nav")
    nav = ["Experiment Group PCA", "Run Overview", "Estimated Search Tolerances", "Mass Accuracy & Tolerance Details", "Putative Modification Mass Shifts"]
    for i, s_ in enumerate(nav):
        yy = by + 94 + i * 46
        if i == 1:
            f.rect(sx - 8, yy - 17, 3, 24, fill=TEAL)
        words = s_.split(" ")
        # wrap at ~20 chars
        lines, cur = [], ""
        for w_ in words:
            if len(cur) + len(w_) + 1 > 21 and cur:
                lines.append(cur); cur = w_
            else:
                cur = (cur + " " + w_).strip()
        lines.append(cur)
        for j, l in enumerate(lines):
            f.text(sx, yy + j * 16, l, size=12.8, fill=INK if i == 1 else INK2, weight="600" if i == 1 else "400", tag="nav")
    f.line(bx + 190, by + 50, bx + 190, by + bh - 16, GRID, 1.2)

    # main area: 2 x 2 plots
    mx0 = bx + 236
    cw, ch = 230, 130
    cxs = [mx0, mx0 + 290]
    cys = [by + 96, by + 96 + 270]

    def title(x, y, t, sub):
        f.text(x, y - 36, t, size=14.5, weight="700", tag="pt")
        f.text(x, y - 18, sub, size=12.8, fill=MUTED, tag="pt")

    # ---- 1 outlier run: MS2 spectra per run
    x, y = cxs[0], cys[0]
    title(x, y, "Run-to-run distribution", "MS2 spectra per run")
    n = 12
    vals = rng.normal(26000, 1500, n); vals[7] = 9800
    ys = Scale(0, 32000, y + ch, y); xs = Scale(0, n, x, x + cw)
    axes(f, x, y, cw, ch, xs, ys, [], [0, 15000, 30000], yfmt=lambda v: f"{v/1000:g}k", xlabel="", ylabel_dx=-40, tag="r1", xlabel_dy=24, size=12.5)
    med = float(np.median(vals))
    f.line(x, ys(med), x + cw, ys(med), SLATE, 1.4, dash="5 4")
    for i, v in enumerate(vals):
        out = i == 7
        f.rect(xs(i) + 3, ys(v), xs(1) - xs(0) - 6, ys(0) - ys(v), fill=AMBER if out else BLUE)
    f.text(xs(7.5), ys(0) + 22, "outlier run", size=13, weight="700", fill=AMBER, anchor="middle", tag="r1")

    # ---- 2 chromatography duration distribution (box + strip)
    x, y = cxs[1], cys[0]
    title(x, y, "Chromatography", "run duration (min)")
    d = np.concatenate([rng.normal(90, 4.0, 11), [61]])
    ys2 = Scale(50, 100, y + ch, y)
    axes(f, x, y, cw, ch, Scale(0, 1, x, x + cw), ys2, [], [60, 80, 100], xticklabels=False, ylabel_dx=-36, tag="r2", size=12.5)
    q1, q2, q3 = np.percentile(d, [25, 50, 75])
    cx_ = x + cw * 0.38
    f.rect(cx_ - 26, ys2(q3), 52, ys2(q1) - ys2(q3), fill=BLUE_L, stroke=BLUE, sw=1.6)
    f.line(cx_ - 26, ys2(q2), cx_ + 26, ys2(q2), BLUE, 2.4)
    for v in d:
        jit = rng.uniform(-14, 14)
        out = v < 70
        f.circle(cx_ + 80 + jit * 0.8, ys2(v), 4.5, fill=AMBER if out else TEAL, stroke=PAPER, sw=1)
    f.text(cx_ + 94, ys2(61) + 4, "short run", size=13, weight="700", fill=AMBER, tag="r2", halo=True)

    # ---- 3 acquisition: charge fractions per run
    x, y = cxs[0], cys[1]
    title(x, y, "Acquisition characteristics", "precursor charge fractions")
    nr = 10
    ys3 = Scale(0, 1, y + ch, y); xs3 = Scale(0, nr, x, x + cw)
    axes(f, x, y, cw, ch, xs3, ys3, [], [0, 0.5, 1], yfmt=lambda v: f"{v:g}", xlabel="run", ylabel_dx=-36, tag="r3", xlabel_dy=24, size=12.5)
    for i in range(nr):
        base = np.array([0.04, 0.52, 0.30, 0.09, 0.05]) + rng.normal(0, 0.012, 5)
        if i == 4:
            base = np.array([0.03, 0.30, 0.38, 0.19, 0.10])
        base = np.clip(base, 0.01, None); base /= base.sum()
        acc = 0
        for (k, col), v in zip(CHC.items(), base):
            f.rect(xs3(i) + 3, ys3(acc + v), xs3(1) - xs3(0) - 6, ys3(acc) - ys3(acc + v), fill=col)
            acc += v
    lx = x
    for k, col in CHC.items():
        f.rect(lx, y + ch + 36, 10, 10, fill=col); f.text(lx + 14, y + ch + 46, k, size=12.5, fill=INK2, tag="r3")
        lx += 14 + text_width(k, 12.5) + 16

    # ---- 4 experiment group PCA
    x, y = cxs[1], cys[1]
    title(x, y, "Experiment Group PCA", "runs in robustly scaled feature space")
    ax_ = Scale(-3, 4.5, x, x + cw); ay_ = Scale(-2.5, 2.5, y + ch, y)
    axes(f, x, y, cw, ch, ax_, ay_, [], [], grid=False, tag="r4")
    f.text(x + cw / 2, y + ch + 22, "PC1", size=12.5, fill=INK2, anchor="middle", tag="r4")
    g1 = rng.normal([-0.8, 0.0], [0.7, 0.7], (14, 2)); g2 = rng.normal([2.8, 0.6], [0.5, 0.6], (6, 2))
    for p in g1: f.circle(ax_(p[0]), ay_(p[1]), 5, fill=TEAL, stroke=PAPER, sw=1.2)
    for p in g2: f.circle(ax_(p[0]), ay_(p[1]), 5, fill=AMBER, stroke=PAPER, sw=1.2)
    f.text(ax_(-0.8), ay_(-1.9), "group 1", size=13, weight="700", fill=TEAL, anchor="middle", tag="r4")
    f.text(ax_(2.8), ay_(-1.5), "group 2", size=13, weight="700", fill=AMBER, anchor="middle", tag="r4")

    f.text(bx + 18, by + bh - 20, "Synthetic example. Section names are those of the mzQC report module; further sections cover mass-shift families.", size=12.8, fill=MUTED, tag="bn")

    # ---------------------------------------------------------- footer
    f.line(60, 772, W - 60, 772, RULE, 1)
    f.text(60, 802, "The report is for inspection and cross-run comparison. SDRF refinement reads each run's summary.json, not this HTML.", size=14.5, weight="600", tag="foot")
    return f


if __name__ == "__main__":
    # Use the common build driver for portable output paths and layout linting.
    from build_all import main
    raise SystemExit(main(["--figures", "6"]))
