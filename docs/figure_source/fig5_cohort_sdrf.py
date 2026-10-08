"""Figure 5 - cohort synthesis and SDRF refinement (synthetic cohort; real gate helpers and UniMod lookups)."""

from __future__ import annotations

from math import ceil

import numpy as np
from figlib import (
    AMBER,
    AMBER_L,
    FAINT,
    GREEN,
    GREEN_L,
    GRID,
    INK2,
    MUTED,
    PAPER,
    PLUM,
    PLUM_L,
    RED,
    RED_L,
    RULE,
    SANS,
    SLATE,
    TEAL,
    TEAL_L,
    WASH,
    Fig,
    Scale,
    axes,
    chip,
    panel_label,
    text_width,
)
from impl import mass_shift as ms

from prideqc import cohort

rng = np.random.default_rng(5)


def mark(f, x, y, kind, r=8):
    """kind: ok | no | na"""
    if kind == "ok":
        f.circle(x, y, r, fill=GREEN)
        f.path(
            f"M{x - r * 0.45:.1f},{y:.1f} l{r * 0.32:.1f},{r * 0.38:.1f} l{r * 0.62:.1f},{-r * 0.75:.1f}",
            stroke=PAPER,
            sw=2,
        )
    elif kind == "no":
        f.circle(x, y, r, fill=RED)
        d = r * 0.4
        f.path(
            f"M{x - d:.1f},{y - d:.1f} L{x + d:.1f},{y + d:.1f} M{x + d:.1f},{y - d:.1f} L{x - d:.1f},{y + d:.1f}",
            stroke=PAPER,
            sw=2,
        )
    else:
        f.line(x - 5, y, x + 5, y, FAINT, 2)


def doc(f, x, y, w=26, h=32, fill=PAPER, stroke=SLATE):
    c = 8
    f.poly(
        [(x, y), (x + w - c, y), (x + w, y + c), (x + w, y + h), (x, y + h)],
        fill=fill,
        stroke=stroke,
        sw=1.4,
        close=True,
    )
    f.poly([(x + w - c, y), (x + w - c, y + c), (x + w, y + c)], stroke=stroke, sw=1.2)


def build() -> Fig:
    W, H = 1400, 1760
    f = Fig(
        W,
        H,
        "Cohort evidence synthesis and conservative SDRF refinement",
        "Runs are grouped by acquisition features. Within a group, a shared tolerance is proposed from the maximum supported per-run estimate when "
        "coverage and unit gates pass, and recurrent mass-shift families must pass prevalence, identity, category and independent evidence gates. "
        "Accepted proposals are applied to the SDRF, validated and logged; uncertain evidence abstains or is held.",
    )
    f.text(
        60,
        54,
        "From per-run evidence to conservative SDRF proposals",
        size=28,
        weight="700",
        tag="title",
    )
    f.text(
        60,
        82,
        "Synthetic cohort. Per-run evidence is never study-wide truth: it must pass group-level gates, and every gate can abstain.",
        size=15,
        fill=MUTED,
        tag="title",
    )

    # ================================================================ A: cohort evidence layer
    panel_label(
        f,
        60,
        138,
        "A",
        "Cohort evidence layer: runs are grouped by acquisition features first",
        "gates below run inside one inferred experiment group (group 1 is shown)",
    )
    # stacked per-run summaries
    for i in range(5):
        doc(f, 70 + i * 7, 196 - i * 6 + 12, fill=PAPER if i else WASH)
    f.text(60, 268, "26 per-run", size=13.5, weight="600", tag="A")
    f.text(60, 286, "summaries", size=13.5, weight="600", tag="A")
    f.arrow([(160, 210), (232, 210)], color=SLATE, sw=2, tag="A_in")
    # scatter in scaled feature space (real two-means + silhouette)
    n1, n2 = 20, 6
    g1 = rng.normal([0.0, 0.0], [0.45, 0.50], (n1, 2))
    g2 = rng.normal([3.4, 1.9], [0.40, 0.45], (n2, 2))
    M = np.vstack([g1, g2])
    labels = cohort._two_means(M)
    sil = cohort._silhouette(M, labels)
    sx0, sw, sy0, sh = 270, 440, 176, 150
    xs = Scale(-1.6, 5.0, sx0, sx0 + sw)
    ys = Scale(-1.6, 3.0, sy0 + sh, sy0)
    axes(f, sx0, sy0, sw, sh, xs, ys, [], [], grid=False, tag="A", left_axis=True)
    f.text(
        sx0 + sw / 2,
        sy0 + sh + 24,
        "chromatography duration (robustly scaled)",
        size=13.5,
        fill=INK2,
        anchor="middle",
        tag="A",
    )
    cx = sx0 - 22
    cy = sy0 + sh / 2
    f.raw(
        f'<text x="{cx}" y="{cy}" transform="rotate(-90 {cx} {cy})" font-family="{SANS}" font-size="13.5" fill="{INK2}" text-anchor="middle">precursor tol. (scaled)</text>'
    )
    lab_main = labels[0]
    for (x_, y_), group_label in zip(M, labels, strict=False):
        f.circle(
            xs(x_), ys(y_), 6, fill=TEAL if group_label == lab_main else AMBER, stroke=PAPER, sw=1.5
        )
    for pts, col in ((M[labels == lab_main], TEAL), (M[labels != lab_main], AMBER)):
        mx, my = pts.mean(0)
        f.path(
            f"M{xs(mx) - 62:.1f},{ys(my):.1f} a62,46 0 1,0 124,0 a62,46 0 1,0 -124,0",
            stroke=col,
            sw=1.4,
            dash="4 4",
            fill="none",
            opacity=0.9,
        )
    f.text(
        xs(0.0) + 76, ys(-0.85), "group 1 (20 runs)", size=13.5, weight="700", fill=TEAL, tag="A"
    )
    f.text(
        xs(3.4),
        ys(1.9) + 66,
        "group 2 (6 runs)",
        size=13.5,
        weight="700",
        fill=AMBER,
        anchor="middle",
        tag="A",
    )
    # criteria
    kx = 780
    f.text(kx, 186, "A split is accepted only if", size=14.5, weight="700", tag="A")
    crit = [
        ("silhouette ≥ 0.45", f"here {sil:.2f}"),
        ("both groups ≥ minimum size", "20 and 6 runs"),
        ("≥ 1 feature clearly separated", "median ratio above its limit"),
    ]
    for i, (a, b) in enumerate(crit):
        yy = 214 + i * 28
        mark(f, kx + 8, yy - 5, "ok", 7)
        f.text(kx + 26, yy, a, size=14, tag="A")
        f.text(kx + 330, yy, b, size=13, fill=MUTED, tag="A")
    f.text(
        kx,
        308,
        "Strong chromatography-duration regimes use a dedicated split (silhouette ≥ 0.60).",
        size=13,
        fill=MUTED,
        tag="A",
    )
    f.text(
        kx,
        328,
        "Features: tolerance estimates, fragment unit and regime, duration, isolation width, MS1/MS2 counts,",
        size=13,
        fill=MUTED,
        tag="A",
    )
    f.text(
        kx,
        346,
        "instrument and acquisition metadata when complete. A numeric feature needs ≥ 80% of runs observed,",
        size=13,
        fill=MUTED,
        tag="A",
    )
    f.text(
        kx,
        364,
        "then log10, median-centring and IQR scaling. Groups reflect acquisition, not biology.",
        size=13,
        fill=MUTED,
        tag="A",
    )

    # ================================================================ B: tolerance branch
    yB = 440
    f.line(60, yB - 36, W - 60, yB - 36, RULE, 1)
    f.rect(60, yB - 22, 6, 26, fill=TEAL)
    f.text(78, yB, "B", size=21, weight="700", fill=TEAL, tag="B")
    f.text(106, yB, "Tolerance branch", size=19, weight="700", fill=TEAL, tag="B")
    f.text(
        296,
        yB,
        "one shared proposal per group, for precursor (ppm) and fragment (ppm or Da)",
        size=14,
        fill=MUTED,
        tag="B",
    )
    n = 20
    supported = 17
    need = max(1, ceil(0.8 * n))
    heads = [
        ("1", "Coverage", 60),
        ("2", "Compatible unit", 372),
        ("3", "Maximum, not mean", 630),
        ("4", "Round upward", 1000),
        ("5", "Shared proposal", 1186),
    ]
    hy = yB + 50
    for _k, (num, title, x) in enumerate(heads):
        f.circle(x + 11, hy - 5, 11, fill=TEAL)
        f.text(x + 11, hy, num, size=13.5, weight="700", fill=PAPER, anchor="middle", tag="Bh")
        f.text(x + 30, hy, title, size=15.5, weight="700", tag="Bh")
    # arrows between headers
    for k in range(4):
        xa = heads[k][2] + 30 + text_width(heads[k][1], 15.5, "700") + 10
        xb = heads[k + 1][2] - 6
        f.arrow([(xa, hy - 5), (xb, hy - 5)], color=TEAL, sw=2, tag=f"Br{k}", hs=6)
    # 1 coverage dots
    cx0, cy0 = 72, hy + 36
    for i in range(n):
        r_, c_ = divmod(i, 10)
        x_, y_ = cx0 + c_ * 27, cy0 + r_ * 28
        if i < supported:
            f.circle(x_ + 8, y_ + 8, 8, fill=TEAL)
        else:
            f.circle(x_ + 8, y_ + 8, 7, fill=PAPER, stroke=RED, sw=2)
    f.text(
        cx0,
        cy0 + 74,
        f"{supported} of {n} runs have a supported estimate",
        size=13.5,
        weight="600",
        tag="B1",
    )
    f.text(
        cx0,
        cy0 + 94,
        f"need ≥ 80%: max(1, ceil(0.8 × {n})) = {need} runs",
        size=13,
        fill=INK2,
        tag="B1",
    )
    mark(f, cx0 + 8, cy0 + 118, "ok", 7)
    f.text(
        cx0 + 24,
        cy0 + 123,
        f"{supported} ≥ {need}: pass",
        size=13.5,
        weight="700",
        fill=GREEN,
        tag="B1",
    )
    # 2 unit
    ux = 384
    f.text(ux, cy0 + 12, "all contributing runs", size=13.5, fill=INK2, tag="B2")
    f.text(ux, cy0 + 30, "report the same unit", size=13.5, fill=INK2, tag="B2")
    chip(f, ux, cy0 + 62, f"ppm × {supported}", fg=TEAL, bg=TEAL_L, tag="B2c")
    mark(f, ux + 86, cy0 + 62, "ok", 7)
    f.text(ux, cy0 + 100, "mixed ppm + Da → abstain", size=13.5, weight="600", fill=RED, tag="B2")
    f.text(ux, cy0 + 120, "(precursor tolerance is ppm only)", size=13, fill=MUTED, tag="B2")
    # 3 MAX vs MEAN strip
    tol = np.array(
        [
            6.1,
            6.8,
            7.4,
            8.0,
            8.3,
            8.9,
            9.4,
            9.9,
            10.5,
            11.0,
            11.8,
            12.4,
            13.0,
            13.9,
            15.2,
            16.8,
            18.4,
        ]
    )
    assert len(tol) == supported
    mean_t, max_t = tol.mean(), tol.max()
    chx0, chw = 640, 330
    cxs = Scale(4, 22, chx0, chx0 + chw)
    cby = cy0 + 86
    f.rect(
        cxs(mean_t), cy0 - 8, cxs(max_t) - cxs(mean_t), cby - cy0 + 8, fill=AMBER_L, opacity=0.55
    )
    f.line(chx0, cby, chx0 + chw, cby, INK2, 1.2)
    for v in (4, 8, 12, 16, 20):
        f.line(cxs(v), cby, cxs(v), cby + 5, INK2, 1)
        f.text(cxs(v), cby + 21, f"{v}", size=12.5, fill=INK2, anchor="middle", tag="B3")
    f.text(
        chx0 + chw,
        cby + 40,
        "per-run suggested tolerance (ppm)",
        size=13,
        fill=INK2,
        anchor="end",
        tag="B3",
    )
    placed = []
    for v in tol:
        lvl = 0
        while any(abs(cxs(v) - px) < 13 and level == lvl for px, level in placed):
            lvl += 1
        placed.append((cxs(v), lvl))
        f.circle(cxs(v), cby - 12 - lvl * 14, 5.5, fill=TEAL, stroke=PAPER, sw=1.2)
    f.line(cxs(mean_t), cy0 - 8, cxs(mean_t), cby, SLATE, 1.8, dash="5 4")
    f.line(cxs(max_t), cy0 - 8, cxs(max_t), cby, GREEN, 2.4)
    f.text(
        cxs(mean_t) - 7,
        cy0 + 6,
        f"mean {mean_t:.1f}",
        size=13.5,
        weight="600",
        fill=SLATE,
        anchor="end",
        tag="B3",
        halo=True,
    )
    f.text(
        cxs(max_t) + 7,
        cy0 + 6,
        f"max {max_t:.1f}",
        size=13.5,
        weight="700",
        fill=GREEN,
        tag="B3",
        halo=True,
    )
    cov_mean = int((tol <= mean_t).sum())
    f.text(
        chx0,
        cby + 64,
        f"a mean-based value leaves {supported - cov_mean} of {supported} supported runs above it",
        size=13.5,
        weight="600",
        fill=AMBER,
        tag="B3",
    )
    f.text(
        chx0,
        cby + 84,
        f"the maximum covers all {supported}: a shared value must cover the noisiest run",
        size=13.5,
        weight="600",
        fill=GREEN,
        tag="B3",
    )
    # 4 rounding
    rx = 1000
    f.text(rx, cy0 + 12, "ppm → whole ppm", size=13.5, fill=INK2, tag="B4")
    f.text(
        rx,
        cy0 + 40,
        f"{max_t:.1f}  →  {cohort._format_tolerance(max_t, 'ppm')}",
        size=17,
        mono=False,
        weight="700",
        fill=TEAL,
        tag="B4",
    )
    f.text(rx, cy0 + 74, "Da → next 0.001 Da", size=13.5, fill=INK2, tag="B4")
    f.text(
        rx,
        cy0 + 100,
        f"0.0137  →  {cohort._format_tolerance(0.0137, 'Da')}",
        size=17,
        weight="700",
        fill=TEAL,
        tag="B4",
    )
    f.text(rx, cy0 + 126, "always upward, never down", size=13, fill=MUTED, tag="B4")
    # 5 output
    ox = 1196
    f.text(ox, cy0 + 12, "comment[precursor", size=13, mono=True, tag="B5")
    f.text(ox, cy0 + 30, "mass tolerance]", size=13, mono=True, tag="B5")
    f.text(
        ox,
        cy0 + 66,
        cohort._format_tolerance(max_t, "ppm"),
        size=26,
        weight="700",
        fill=GREEN,
        tag="B5",
    )
    f.text(ox, cy0 + 92, "for group 1 runs", size=13, fill=MUTED, tag="B5")
    # abstain note
    f.text(
        60,
        cy0 + 206,
        "Abstain (no proposal): < 80% of runs supported, or mixed units. A proposed value is a reanalysis setting, not a record of the original search.",
        size=13.5,
        weight="600",
        fill=RED,
        tag="B",
    )

    # ================================================================ C: PTM branch matrix
    yC = cy0 + 270
    f.line(60, yC - 36, W - 60, yC - 36, RULE, 1)
    f.rect(60, yC - 22, 6, 26, fill=PLUM)
    f.text(78, yC, "C", size=21, weight="700", fill=PLUM, tag="C")
    f.text(106, yC, "Modification branch", size=19, weight="700", fill=PLUM, tag="C")
    f.text(
        298,
        yC,
        "recurrent mass-shift families, gated in order; a family that fails a gate stops there",
        size=14,
        fill=MUTED,
        tag="C",
    )
    colx = [380, 610, 810, 1010]
    heads = [
        (
            "1  Recurrence and prevalence",
            [
                "≥ 3 runs; prevalence ≥ 0.90",
                "≥ 0.80 high-support runs",
                "P(prevalence > 0.10) ≥ 0.99",
            ],
        ),
        (
            "2  Unique identity",
            [
                "exactly one candidate",
                "(decoys, substitutions excluded)",
                "seen in ≥ 0.80 of hit runs",
            ],
        ),
        ("3  Category", ["candidate is a", "biological PTM"]),
        ("4  Independent evidence", ["study evidence supports", "this exact UniMod entry"]),
    ]
    hy2 = yC + 40
    for (t, lines), x in zip(heads, colx, strict=False):
        f.text(x, hy2, t, size=14.5, weight="700", tag="Ch")
        for i, label in enumerate(lines):
            f.text(x, hy2 + 20 + i * 17, label, size=12.8, fill=MUTED, tag="Ch")
    f.text(1222, hy2, "Outcome", size=14.5, weight="700", tag="Ch")
    # rows with real lookups
    recs = ms.load_openms_modifications()

    def lookup(m):
        c = ms._matching_modifications(m, 0.02, recs, maximum_candidates=24)
        return cohort._sdrf_candidates({"unimod_candidates": c})

    def P(h, r):
        return cohort._beta_prevalence_probability(h, r)

    rows = []
    for mass, hits, ev in (
        (14.01565, 19, "supported"),
        (43.00581, 19, "not-evaluated"),
        (79.96633, 20, ""),
        (71.03711, 18, ""),
        (57.02146, 6, ""),
    ):
        cands = lookup(mass)
        rows.append((mass, hits, cands, ev))
    ry0 = hy2 + 100
    RH = 62
    for i, (mass, hits, cands, ev) in enumerate(rows):
        y = ry0 + i * RH
        f.line(60, y - 24, W - 60, y - 24, GRID, 1)
        prev = hits / 20
        p = P(hits, 20)
        ok1 = prev >= 0.90 and p >= 0.99
        uniq = len(cands) == 1
        acc, (nm, cat) = next(iter(cands.items())) if uniq else ("", ("", ""))
        ok2 = ok1 and uniq
        ok3 = ok2 and cat == "biological-ptm"
        ok4 = ok3 and ev == "supported"
        f.text(60, y + 2, f"+{mass:.4f} Da", size=15.5, weight="700", tag="Cm")
        f.text(
            60,
            y + 21,
            (f"{nm} ({acc})" if uniq else f"{len(cands)} catalogue candidates"),
            size=12.8,
            fill=MUTED,
            tag="Cm",
        )
        # gate 1
        mark(f, colx[0] + 8, y - 4, "ok" if ok1 else "no")
        f.text(
            colx[0] + 24, y, f"{hits}/20 runs, prev {prev:.2f}", size=13.5, weight="600", tag="Cg"
        )
        f.text(
            colx[0] + 24,
            y + 18,
            f"P = {p:.3f}" if p >= 0.0005 else "P < 0.001",
            size=12.8,
            fill=MUTED,
            tag="Cg",
        )
        # gate 2
        if ok1:
            mark(f, colx[1] + 8, y - 4, "ok" if uniq else "no")
            f.text(
                colx[1] + 24,
                y,
                (acc if uniq else f"{len(cands)} candidates"),
                size=13.5,
                weight="600",
                tag="Cg",
            )
            f.text(
                colx[1] + 24,
                y + 18,
                "unique" if uniq else "same mass window",
                size=12.8,
                fill=MUTED,
                tag="Cg",
            )
        else:
            mark(f, colx[1] + 8, y - 4, "na")
        # gate 3
        if ok2:
            mark(f, colx[2] + 8, y - 4, "ok" if ok3 else "no")
            f.text(
                colx[2] + 24,
                y,
                cat.replace("-or-artifact", "/artefact").replace(
                    "biological-ptm", "biological PTM"
                ),
                size=13.5,
                weight="600",
                tag="Cg",
            )
        else:
            mark(f, colx[2] + 8, y - 4, "na")
        # gate 4
        if ok3:
            mark(f, colx[3] + 8, y - 4, "ok" if ok4 else "no")
            f.text(colx[3] + 24, y, ev, size=13.5, weight="600", tag="Cg")
        else:
            mark(f, colx[3] + 8, y - 4, "na")
        # outcome
        if ok4:
            col, t1, t2 = GREEN, "ELIGIBLE", "standard column"
        elif ok3:
            col, t1, t2 = AMBER, "HOLD", "putative column only"
        else:
            col, t1, t2 = RED, "ABSTAIN", "no write"
        f.rect(1210, y - 20, 4, 44, fill=col)
        f.text(1224, y, t1, size=14.5, weight="700", fill=col, tag="Co")
        f.text(1224, y + 18, t2, size=12.8, fill=MUTED, tag="Co")
    yEnd = ry0 + 5 * RH - 24
    f.line(60, yEnd, W - 60, yEnd, GRID, 1)
    f.text(
        60,
        yEnd + 26,
        "Masses and candidate names are real catalogue look-ups; run counts are synthetic. Up to the 5 best-supported families per group are kept.",
        size=13,
        fill=MUTED,
        tag="C",
    )
    f.text(
        60,
        yEnd + 46,
        "Hold and eligible families are both written to the putative column with their evidence fields (runs, prevalence, support, status); only eligible ones reach the standard column.",
        size=13,
        fill=MUTED,
        tag="C",
    )

    # ================================================================ D: outcome
    yD = yEnd + 112
    f.line(60, yD - 36, W - 60, yD - 36, RULE, 1)
    panel_label(
        f,
        60,
        yD,
        "D",
        "Apply, validate and keep a record",
        "accepted proposals only; the original SDRF is the starting point",
    )
    ty = yD + 66
    cols = ["run", "precursor tol.", "modification", "putative"]
    cw = [46, 92, 92, 70]

    def table(x, y, rows_, kinds, title):
        f.text(x, y - 12, title, size=14.5, weight="700", tag="Dt")
        xx = x
        for h, w in zip(cols, cw, strict=False):
            f.text(xx + 4, y + 14, h, size=12.5, weight="600", fill=INK2, tag="Dt")
            xx += w
        f.line(x, y + 22, x + sum(cw), y + 22, INK2, 1.3)
        for r, (cells, ks) in enumerate(zip(rows_, kinds, strict=False)):
            yy = y + 24 + r * 34
            xx = x
            for _c, (txt, k, w) in enumerate(zip(cells, ks, cw, strict=False)):
                if k:
                    colr = {
                        "filled": GREEN_L,
                        "replaced": AMBER_L,
                        "appended": PLUM_L,
                        "conflict": RED_L,
                        "hold": PLUM_L,
                    }[k]
                    f.rect(xx + 1, yy + 2, w - 2, 30, fill=colr)
                    if k == "conflict":
                        f.rect(xx + 1, yy + 2, w - 2, 30, fill="none", stroke=RED, sw=1.4)
                f.text(xx + 4, yy + 22, txt, size=12.8, tag="Dt")
                xx += w
            f.line(x, yy + 34, x + sum(cw), yy + 34, GRID, 1)
        return y + 24 + len(rows_) * 34

    orig = [("run1", "—", "—", "—"), ("run2", "10 ppm", "—", "—"), ("run3", "19 ppm", "—", "—")]
    new = [
        ("run1", "19 ppm", "Methylation", "putative"),
        ("run2", "19 ppm", "Methylation", "putative"),
        ("run3", "19 ppm", "Methylation", "putative"),
    ]
    kinds = [
        (None, "filled", "appended", "appended"),
        (None, "replaced", "appended", "appended"),
        (None, None, "appended", "appended"),
    ]
    x_o = 60
    table(x_o, ty, orig, [(None,) * 4] * 3, "original SDRF")
    # "+" then proposals
    f.text(
        x_o + sum(cw) + 22,
        ty + 70,
        "+",
        size=26,
        weight="700",
        fill=SLATE,
        anchor="middle",
        tag="Dp",
    )
    px = x_o + sum(cw) + 34
    f.text(px, ty - 12, "accepted proposals", size=14.5, weight="700", tag="Dp")
    f.circle(px + 6, ty + 20, 6, fill=TEAL)
    f.text(px + 20, ty + 25, "19 ppm (group 1)", size=13, tag="Dp")
    f.circle(px + 6, ty + 48, 6, fill=PLUM)
    f.text(px + 20, ty + 53, "Methylation: eligible", size=13, tag="Dp")
    f.circle(px + 6, ty + 76, 6, fill=PAPER, stroke=PLUM, sw=2)
    f.text(px + 20, ty + 81, "putative: kept as evidence", size=13, tag="Dp")
    f.text(px + 20, ty + 99, "abstained: not written", size=13, fill=MUTED, tag="Dp")
    f.circle(px + 6, ty + 94, 0.1, fill=PAPER)
    ax_ = px + 176
    f.arrow([(ax_, ty + 66), (ax_ + 36, ty + 66)], color=SLATE, sw=2.2, tag="Dar1")
    xr = ax_ + 52
    yb2 = table(xr, ty, new, kinds, "refined SDRF")
    # change-kind legend
    lg = [
        ("no colour", PAPER, "already equal (match)"),
        ("filled", GREEN_L, "was empty"),
        ("replaced", AMBER_L, "tolerance updated"),
        ("appended", PLUM_L, "new value"),
        ("conflict", RED_L, "differs: left unchanged, flagged"),
    ]
    lx, lyy = x_o, yb2 + 30
    for name, colr, lab in lg:
        f.rect(
            lx,
            lyy - 11,
            14,
            14,
            fill=colr,
            stroke=RED if name == "conflict" else (FAINT if name == "no colour" else "none"),
            sw=1.4,
        )
        f.text(
            lx + 20,
            lyy,
            f"{name}: {lab}" if name != "no colour" else lab,
            size=13,
            fill=INK2,
            tag="Dl",
        )
        lx += 20 + text_width(f"{name}: {lab}" if name != "no colour" else lab, 13) + 26
    # validation + audit log
    vx = xr + sum(cw) + 40
    f.arrow([(xr + sum(cw) + 8, ty + 66), (vx - 6, ty + 66)], color=SLATE, sw=2.2, tag="Dar2")
    mark(f, vx + 14, ty + 56, "ok", 14)
    f.text(vx + 40, ty + 50, "validation", size=15, weight="700", tag="Dv")
    f.text(vx + 40, ty + 69, "SDRF-pipelines", size=13, fill=MUTED, tag="Dv")
    f.text(vx + 40, ty + 87, "on the refined file", size=13, fill=MUTED, tag="Dv")
    lx2 = vx + 196
    f.arrow([(vx + 164, ty + 66), (lx2 - 10, ty + 66)], color=SLATE, sw=2.2, tag="Dar3")
    f.text(lx2, ty - 12, "change and audit record", size=14.5, weight="700", tag="Dg")
    logs = [
        "sdrf-changes.tsv",
        "sdrf-refinement.log.txt",
        "cohort-refinement.json",
        "manifest.json",
    ]
    for i, s_ in enumerate(logs):
        doc(f, lx2, ty + 8 + i * 30, 18, 22, fill=WASH)
        f.text(lx2 + 28, ty + 25 + i * 30, s_, size=13, mono=True, tag="Dg")
    f.text(lx2, ty + 140, "row, column, old, new, status,", size=13, fill=MUTED, tag="Dg")
    f.text(lx2, ty + 158, "evidence and support", size=13, fill=MUTED, tag="Dg")
    # optional adjudication path (below, dashed, no crossing)
    ya = lyy + 62
    f.line(60, ya - 26, W - 60, ya - 26, GRID, 1)
    f.text(
        60, ya, "Uncertain evidence has a visible exit.", size=15, weight="700", fill=RED, tag="Dx"
    )
    f.text(
        60,
        ya + 24,
        "Abstain and hold outcomes leave the SDRF value unchanged (or as a labelled putative entry). An optional bounded packet of the same evidence can be",
        size=13.5,
        fill=INK2,
        tag="Dx",
    )
    f.text(
        60,
        ya + 44,
        "reviewed for accept / reject / abstain decisions; it does not change any QC computation (llm-refinement-packet.json).",
        size=13.5,
        fill=INK2,
        tag="Dx",
    )
    f.h = ya + 80
    return f


if __name__ == "__main__":
    # Use the common build driver for portable output paths and layout linting.
    from build_all import main

    raise SystemExit(main(["--figures", "5"]))
