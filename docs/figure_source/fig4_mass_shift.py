"""Figure 4 - recurrent mass-shift scouting, computed with prideQC's mass_shift.py on synthetic spectra."""
from __future__ import annotations
import math
import numpy as np
from figlib import *
import synth_shift as S
from impl import mass_shift as ms

X0, X1 = 330, 1340


def step_head(f, y, n, title, sub):
    f.text(60, y, str(n), size=30, weight="700", fill=PLUM, tag="sh")
    f.text(92, y - 6, title, size=16.5, weight="700", tag="sh")
    for i, l in enumerate(sub):
        f.text(92, y + 16 + i * 18, l, size=13.2, fill=MUTED, tag="sh")


def stems(f, xs, mz, inten, base, up, hmax, color, sw=2.0, opacity=None):
    for m, i in zip(mz, inten):
        y2 = base - hmax * i if up else base + hmax * i
        f.line(xs(m), base, xs(m), y2, color, sw, opacity=opacity)


def build() -> Fig:
    d = S.build()
    obs, clusters, info = S.family_data()
    W, H = 1400, 2170
    f = Fig(W, H, "Recurrent mass-shift scouting",
            "Two related synthetic MS2 spectra: strongest peaks, shared coarse bins, unchanged and shifted fragment matches, and the accepted neutral mass difference; "
            "then a recurrent mass-shift family, isotope and adduct check, and tight UniMod mass-compatible lookup. Mass compatibility is not identity.")
    f.text(60, 54, "Scouting recurrent mass shifts without identifications", size=28, weight="700", tag="title")
    f.text(60, 82, "Synthetic spectra scored with prideQC's own mass-shift code. The scout looks for related spectra whose precursor masses differ by a recurring amount.",
           size=15, fill=MUTED, tag="title")
    f.text(60, 128, "I   From two related spectra to one accepted pair", size=19, weight="700", fill=INK, tag="part")
    f.line(60, 140, W - 60, 140, RULE, 1)

    A, B = d["A"], d["B"]
    xs = Scale(100, 1500, X0, X1)

    # --------------------------------------------------------- step 1: strongest peaks
    step_head(f, 190, 1, "Strongest peaks", ["Top 60 peaks per spectrum;", "±2 Da of precursor removed;", "L2-normalised intensities."])
    base = 300
    hmax = 92
    rawA, rawB = d["rawA"], d["rawB"]
    nA = rawA[1] / rawA[1].max(); nB = rawB[1] / rawB[1].max()
    stems(f, xs, rawA[0], nA * 0.5, base - 3, True, hmax, FAINT, 1.2, 0.55)      # all raw peaks, faint
    stems(f, xs, rawB[0], nB * 0.5, base + 3, False, hmax, FAINT, 1.2, 0.55)
    stems(f, xs, A.mz, A.intensity / A.intensity.max(), base - 3, True, hmax, BLUE, 2.0)
    stems(f, xs, B.mz, B.intensity / B.intensity.max(), base + 3, False, hmax, TEAL, 2.0)
    f.line(X0, base, X1, base, INK2, 1.2)
    for v in (200, 400, 600, 800, 1000, 1200, 1400):
        f.line(xs(v), base + hmax + 8, xs(v), base + hmax + 13, INK2, 1)
        f.text(xs(v), base + hmax + 30, f"{v}", size=12.5, fill=INK2, anchor="middle", tag="s1")
    f.line(X0, base + hmax + 8, X1, base + hmax + 8, INK2, 1.2)
    f.text(X1, base + hmax + 50, "fragment m/z", size=13.5, fill=INK2, anchor="end", tag="s1")
    f.text(X1, base - hmax + 2, "spectrum A: no modification", size=14, weight="700", fill=BLUE, anchor="end", tag="s1")
    f.text(X1, base + hmax - 2, "spectrum B: precursor mass shifted", size=14, weight="700", fill=TEAL, anchor="end", tag="s1")
    for pm, up, col, lab in ((d["pA"], True, BLUE, f"precursor {d['pA']:.2f} (2+)"), (d["pB"], False, TEAL, f"precursor {d['pB']:.2f} (2+)")):
        y = base - 18 if up else base + 18
        f.path(f"M{xs(pm)-5:.1f},{y - (6 if up else -6):.1f} l5,{12 if up else -12} l5,{-12 if up else 12} Z", fill=col, stroke="none")
        f.text(xs(pm) - 10, y + (-2 if up else 6), lab, size=12.5, fill=col, anchor="end", tag="s1", halo=True)
    f.text(X0, base + hmax + 50, "faint stems: weaker peaks outside the top 60", size=12.5, fill=MUTED, tag="s1")

    # --------------------------------------------------------- step 2: shared coarse bins
    y2 = 500
    step_head(f, y2 + 20, 2, "Shared coarse bins", ["1 Da fragment bins index", "spectra; ≥ 5 shared bins", "(within 900 s) get scored."])
    for row, (bins, col) in enumerate(((A.bins, BLUE), (B.bins, TEAL))):
        yy = y2 + 4 + row * 26
        for b_ in bins:
            shared = b_ in d["shared"]
            f.rect(xs(b_), yy, 2.2, 18, fill=AMBER if shared else col, opacity=1.0 if shared else 0.45)
    f.text(X0 - 12, y2 + 17, "A bins", size=13, fill=BLUE, anchor="end", tag="s2")
    f.text(X0 - 12, y2 + 43, "B bins", size=13, fill=TEAL, anchor="end", tag="s2")
    f.rect(X1 - 360, y2 + 62, 10, 10, fill=AMBER)
    f.text(X1 - 342, y2 + 72, f"{len(d['shared'])} shared bins ≥ 5: B retrieves A as a candidate", size=13.5, weight="600", fill=AMBER, tag="s2")

    # --------------------------------------------------------- step 3: unchanged / full / half
    y3 = 620
    step_head(f, y3 + 20, 3, "Unchanged and shifted matches",
              ["One-to-one matches (0.05 Da):", "unchanged, shifted by Δ,", "or shifted by Δ/2 (2+ ions)."])
    rows = [("unchanged", "offset 0", TEAL), ("full", f"A + Δ  (Δ = {d['signed']:.3f})", PLUM), ("half", f"A + Δ/2  ({d['signed']/2:.3f})", AMBER)]
    r_base = [y3 + 118, y3 + 118 + 150, y3 + 118 + 300]
    for (key, lab, col), yb in zip(rows, r_base):
        r = d["res"][key]
        off = r["offset"]
        hm = 46
        matchedA = {i for i, j in r["pairs"]}; matchedB = {j for i, j in r["pairs"]}
        for idx, (m, i) in enumerate(zip(A.mz, A.intensity)):
            if xs(m + off) > X1 + 1:
                continue
            hh = hm * i / A.intensity.max()
            f.line(xs(m + off), yb - 2, xs(m + off), yb - 2 - hh, col if idx in matchedA else FAINT, 2.0 if idx in matchedA else 1.2, opacity=1 if idx in matchedA else 0.5)
        for idx, (m, i) in enumerate(zip(B.mz, B.intensity)):
            hh = hm * i / B.intensity.max()
            f.line(xs(m), yb + 2, xs(m), yb + 2 + hh, col if idx in matchedB else FAINT, 2.0 if idx in matchedB else 1.2, opacity=1 if idx in matchedB else 0.5)
        f.line(X0, yb, X1, yb, INK2, 1)
        f.text(X0 - 12, yb - 22, lab.split("  ")[0], size=14, weight="700", fill=col, anchor="end", tag="s3")
        if "(" in lab:
            f.text(X0 - 12, yb - 5, "(" + lab.split("(")[1], size=12.5, fill=MUTED, anchor="end", tag="s3")
        else:
            f.text(X0 - 12, yb - 5, "no shift", size=12.5, fill=MUTED, anchor="end", tag="s3")
        f.text(X0 - 12, yb + 36, f"{r['count']} matches", size=13.5, weight="600", fill=col, anchor="end", tag="s3")
        f.text(X0 - 12, yb + 54, f"dot = {r['dot']:.3f}", size=12.5, fill=MUTED, anchor="end", tag="s3")
        tag = {"unchanged": "always scored", "full": "retained: higher dot product", "half": "evaluated, lower dot: not retained"}[key]
        f.text(X1, yb - 56, tag, size=13, weight="600" if key != "half" else "400", fill=col if key != "half" else MUTED, anchor="end", tag="s3")
    f.text(X0, r_base[2] + 70, "A is drawn translated by the tested offset, so matched peaks sit directly above their partners in B. Only the better of the two shifted tests is added to the score.",
           size=13, fill=MUTED, tag="s3")

    # --------------------------------------------------------- step 4: accepted pair
    y4 = r_base[2] + 112
    step_head(f, y4 + 20, 4, "Accepted mass difference", ["Charge-aware neutral masses;", "absolute difference kept."])
    MA = A.neutral_mass; MB = B.neutral_mass
    f.text(X0, y4 + 44, f"{abs(MB - MA):.3f} Da", size=44, weight="700", fill=PLUM, tag="s4")
    f.text(X0, y4 + 70, "neutral precursor mass difference", size=13.5, fill=MUTED, tag="s4")
    f.text(X0, y4 + 92, "absolute value; direction unresolved", size=13, fill=MUTED, tag="s4")
    gx = 700
    gates = [
        (f"M = m/z × z − z × 1.00728      A {MA:.4f}   B {MB:.4f}", ""),
        (f"|M_B − M_A| within 0.5 – 500 Da", f"{abs(MB - MA):.3f}"),
        ("unchanged matches ≥ 4", f"{d['unchanged']}"),
        ("unchanged + shifted matches ≥ 10", f"{d['total']}"),
        ("similarity = min(1, dot₀ + dot_shifted) ≥ 0.25", f"{d['sim']:.2f}"),
    ]
    for i, (lab, val) in enumerate(gates):
        yy = y4 + 20 + i * 26
        f.text(gx, yy, lab, size=13.5, mono=True, tag="g4")
        if val:
            f.text(X1 - 36, yy, val, size=13.5, mono=True, weight="700", anchor="end", tag="g4")
            f.circle(X1 - 12, yy - 5, 7, fill=GREEN)
            f.path(f"M{X1-15.5:.1f},{yy-5:.1f} l2.6,3 l4.8,-5.8", stroke=PAPER, sw=1.9)
    f.text(gx, y4 + 20 + 5 * 26 + 6, "accepted pair: at most the 2 best-scoring pairs per spectrum are kept", size=13.5, weight="600", fill=PLUM, tag="g4")

    # =========================================================== Part II
    yP = y4 + 214
    f.text(60, yP, "II   From many accepted pairs to a recurrent mass-shift family", size=19, weight="700", tag="part")
    f.line(60, yP + 12, W - 60, yP + 12, RULE, 1)

    # --------------------------------------------------------- step 5: delta-mass distribution
    y5 = yP + 60
    step_head(f, y5 + 14, 5, "Recurrent shifts cluster", ["Accepted |Δ| are grouped;", "a family needs ≥ 8 pairs", "from ≥ 6 distinct spectra."])
    ow = 640
    oxs = Scale(0, 90, X0, X0 + ow)
    oh = 170
    oy = y5 + 96
    bins = np.arange(0, 90.0001, 0.25)
    vals = np.array([o.delta_mass_da for o in obs])
    cnt, _ = np.histogram(vals, bins=bins)
    ys5 = Scale(0, 100, oy + oh, oy)
    axes(f, X0, oy, ow, oh, oxs, ys5, [0, 20, 40, 60, 80], [0, 50, 100], xlabel="accepted neutral mass difference |Δ| (Da)", ylabel="pairs per 0.25 Da", ylabel_dx=-50, tag="s5")
    centers = {round(c.center_da, 3): c for c in clusters}
    f.line(X0, ys5(8), X0 + ow, ys5(8), AMBER, 1.3, dash="5 4")
    f.text(X0 + ow, ys5(8) - 6, "8-pair minimum", size=12.5, weight="600", fill=AMBER, anchor="end", tag="s5", halo=True)
    iso_adduct = {}
    for c in clusters:
        art = ms._artifact_classification(c.center_da, 0.02, 1.00335483507)
        iso_adduct[c.center_da] = art
    for c_, e0 in zip(cnt, bins[:-1]):
        if c_ > 0:
            f.rect(oxs(e0), ys5(c_), max(1.6, oxs(e0 + 0.25) - oxs(e0)), ys5(0) - ys5(c_), fill=FAINT if c_ < 6 else (SLATE if False else PLUM), opacity=0.9)
    # colour clusters by artefact check
    for c in clusters:
        art = iso_adduct[c.center_da]
        col = SLATE if art else PLUM
        n_here = int(((vals > c.minimum_da - 1e-9) & (vals < c.maximum_da + 1e-9)).sum())
        f.rect(oxs(c.center_da) - 1.2, ys5(c.pair_support), 2.4, ys5(0) - ys5(c.pair_support), fill=col)
    srt = sorted(clusters, key=lambda c: c.center_da)
    nlab = len(srt)
    for k, c in enumerate(srt):
        art = iso_adduct[c.center_da]
        xl = X0 + 36 + k * (ow - 72) / (nlab - 1)
        f.line(xl, oy - 30, oy and oxs(c.center_da), ys5(c.pair_support) - 3, FAINT, 1.1) if False else \
            f.line(xl, oy - 24, oxs(c.center_da), ys5(c.pair_support) - 3, FAINT, 1.1)
        f.text(xl, oy - 32, f"{c.center_da:.3f}", size=13, weight="700", fill=SLATE if art else PLUM, anchor="middle", tag="s5l")
    lgx = X0
    for lab, col in (("family, no isotope/adduct match", PLUM), ("isotope-like or adduct-like", SLATE), ("background pairs (no recurrence)", FAINT)):
        f.rect(lgx, oy - 70, 10, 10, fill=col); f.text(lgx + 16, oy - 61, lab, size=12.5, fill=INK2, tag="s5lg")
        lgx += 16 + text_width(lab, 12.5) + 26

    # zoom on the strongest family
    big = max(clusters, key=lambda c: c.pair_support)
    zx0, zw = X0 + ow + 70, X1 - (X0 + ow + 70)
    zxs = Scale(big.center_da - 0.07, big.center_da + 0.07, zx0, zx0 + zw)
    vz = vals[(vals > big.center_da - 0.07) & (vals < big.center_da + 0.07)]
    zb = np.linspace(big.center_da - 0.07, big.center_da + 0.07, 29)
    zc, _ = np.histogram(vz, bins=zb)
    zys = Scale(0, zc.max() * 1.15, oy + oh, oy)
    tolc = info["cluster_tol"]; tolu = 0.02
    axes(f, zx0, oy, zw, oh, zxs, zys, [big.center_da - 0.05, big.center_da, big.center_da + 0.05], [], grid=False,
         xlabel="zoom: |Δ| (Da)", xfmt=lambda v: f"{v:.2f}", tag="s5z")
    f.rect(zxs(big.center_da - tolc), oy - 4, zxs(big.center_da + tolc) - zxs(big.center_da - tolc), oh + 4, fill=AMBER_L, opacity=0.6)
    f.rect(zxs(big.center_da - tolu), oy + 20, zxs(big.center_da + tolu) - zxs(big.center_da - tolu), oh - 16, fill=PLUM_L, opacity=0.8)
    for c_, e0, e1 in zip(zc, zb[:-1], zb[1:]):
        f.rect(zxs(e0) + 0.8, zys(c_), zxs(e1) - zxs(e0) - 1.6, zys(0) - zys(c_), fill=PLUM)
    f.line(zxs(big.center_da), oy - 4, zxs(big.center_da), oy + oh, INK, 1.4, dash="4 3")
    f.text(zx0, oy - 12, "recurrence window", size=12.5, weight="700", fill=AMBER, tag="s5z")
    f.text(zx0, oy + 14, f"±{tolc:.3f} Da", size=12.5, fill=AMBER, tag="s5z")
    f.text(zx0 + zw, oy + 12 + 22, "annotation window", size=12.5, weight="700", fill=PLUM, anchor="end", tag="s5z")
    f.text(zx0 + zw, oy + 12 + 40, "±0.020 Da, fixed", size=12.5, fill=PLUM, anchor="end", tag="s5z")
    # stats under overview/zoom
    sy = oy + oh + 62
    f.text(X0, sy, f"family {big.center_da:.4f} Da   median of {big.pair_support} pairs · {big.unique_spectrum_support} spectra · σ {big.sigma_da:.4f} Da (1.4826 × MAD) · similarity {big.median_similarity:.2f}",
           size=13.5, mono=False, weight="600", tag="s5s")
    f.text(X0, sy + 22, f"recurrence window = max(0.01, 6 × √2 × σ_ppm × M / 10⁶) = {tolc:.3f} Da   (σ_ppm = {info['sigma_ppm']:.1f} from mass-error precision; representative M = {info['M_rep']:.0f} Da).",
           size=13, fill=INK2, tag="s5s")
    f.text(X0, sy + 42, "A wider recurrence window never loosens the chemical annotation window.", size=13, fill=MUTED, tag="s5s")

    # --------------------------------------------------------- step 6: annotation rail
    y6 = sy + 94
    step_head(f, y6 + 14, 6, "Candidate chemistry", ["Isotope or adduct first,", "then a tight UniMod lookup."])
    nodes = [("Recurrent family", f"{big.center_da:.3f} Da"), ("Isotope / adduct check", "±0.020 Da"), ("UniMod lookup", "OpenMS catalogue, ±0.020 Da"), ("Candidate category", "")]
    nx = [X0 + 20, X0 + 215, X0 + 450, X0 + 700]
    ny = y6 + 26
    for k, (t, s_) in enumerate(nodes):
        f.circle(nx[k], ny, 11, fill=PLUM)
        f.text(nx[k], ny + 5, str(k + 1) if False else "abcd"[k], size=13, weight="700", fill=PAPER, anchor="middle", tag="nn")
        f.text(nx[k] + 22, ny + 5, t, size=15, weight="700", tag="nt")
    for k in range(3):
        xa = nx[k] + 22 + text_width(nodes[k][0], 15, "700") + 12
        xb = nx[k + 1] - 20
        f.arrow([(xa, ny), (xb, ny)], color=PLUM, sw=2, tag=f"rail{k}", hs=6)
    # details under nodes
    f.text(nx[1] - 11, ny + 38, "C13 isotope ×1–4 (1.00335)", size=13, fill=INK2, tag="nd")
    f.text(nx[1] - 11, ny + 56, "Na / K-for-H adduct shifts", size=13, fill=INK2, tag="nd")
    f.text(nx[1] - 11, ny + 78, "no match for this family", size=13, weight="700", fill=GREEN, tag="nd")
    f.text(nx[2] - 11, ny + 38, "mass-compatible entries only;", size=13, fill=INK2, tag="nd")
    f.text(nx[2] - 11, ny + 56, "decoys and substitutions kept", size=13, fill=INK2, tag="nd")
    f.text(nx[2] - 11, ny + 74, "out of the default list", size=13, fill=INK2, tag="nd")
    recs = ms.load_openms_modifications()
    cand = ms._matching_modifications(big.center_da, 0.02, recs, maximum_candidates=24)
    vis = ms._display_modification_candidates(cand, maximum_candidates=8)
    cats = sorted({c["candidate_category"] for c in vis})
    cl = "putative-ptm" if "biological-ptm" in cats else "sample-prep-modification"
    f.text(nx[3] - 11, ny + 38, "biological PTM · sample-prep / artefact ·", size=13, fill=INK2, tag="nd")
    f.text(nx[3] - 11, ny + 56, "other · (substitution, decoy: suppressed)", size=13, fill=INK2, tag="nd")
    f.text(nx[3] - 11, ny + 78, f"here: {len(vis)} biological-PTM candidates → {cl}", size=13, weight="700", fill=PLUM, tag="nd")

    # number line of real catalogue entries near the family
    ly = ny + 190
    nl_x0, nl_x1 = X0 + 20, X1 - 20
    rxs = Scale(-0.12, 0.12, nl_x0, nl_x1)
    f.text(X0 - 12, ly - 20, "catalogue entries", size=13.5, weight="700", fill=INK2, anchor="end", tag="nl")
    f.text(X0 - 12, ly - 2, "near this mass", size=13, fill=MUTED, anchor="end", tag="nl")
    f.rect(rxs(-0.02), ly - 78, rxs(0.02) - rxs(-0.02), 100, fill=PLUM_L, opacity=0.8)
    f.line(nl_x0, ly + 22, nl_x1, ly + 22, INK2, 1.2)
    for v in (-0.10, -0.05, 0, 0.05, 0.10):
        f.line(rxs(v), ly + 22, rxs(v), ly + 27, INK2, 1)
        f.text(rxs(v), ly + 44, f"{v:+.2f}" if v else "0", size=12.5, fill=INK2, anchor="middle", tag="nl")
    f.text(nl_x1, ly + 64, "catalogue mass − observed family mass (Da)", size=13, fill=INK2, anchor="end", tag="nl")
    f.text(rxs(0), ly - 84, "±0.020 Da annotation window", size=13, weight="700", fill=PLUM, anchor="middle", tag="nl")
    near = []
    for r in recs:
        res = abs(r.delta_mass_da) - big.center_da
        if abs(res) <= 0.12:
            near.append((res, r))
    near.sort(key=lambda t: t[0])
    CAT = {"biological-ptm": PLUM, "sample-prep-or-artifact": AMBER, "other-modification": SLATE, "decoy": FAINT, "amino-acid-substitution": FAINT}
    LAY = {"UniMod:423": (24, "left"), "UniMod:40": (46, "left"), "UniMod:21": (46, "right"), "UniMod:99913": (70, "right"),
           "UniMod:1927": (24, "right"), "UniMod:1117": (46, "right"), "UniMod:1104": (70, "right")}
    for res, r in near:
        cat = ms._candidate_category(r)
        col = CAT.get(cat, SLATE)
        suppressed = cat in ("decoy", "amino-acid-substitution")
        h, side = LAY.get(r.accession, (24, "right"))
        x = rxs(res)
        f.line(x, ly + 22, x, ly + 22 - h, col, 1.6)
        if suppressed:
            f.circle(x, ly + 22 - h, 5, fill=PAPER, stroke=FAINT, sw=1.6)
        else:
            f.circle(x, ly + 22 - h, 5.5, fill=col)
        name = r.name if len(r.name) <= 26 else r.name[:24] + "…"
        inside = abs(res) <= 0.02
        txt = f"{name}  {res:+.4f}" if inside else name
        f.text(x + (10 if side == "right" else -10), ly + 22 - h + 5, txt, size=13 if inside else 12.5,
               weight="700" if (inside and not suppressed) else "400", fill=(INK if not suppressed else MUTED) if inside else MUTED,
               anchor="start" if side == "right" else "end", tag="nlab", halo=True)
    # category legend
    lx = X0 + 20
    ly2 = ly + 98
    for lab, col, hollow in (("biological PTM", PLUM, False), ("sample-prep or artefact", AMBER, False), ("other modification", SLATE, False), ("decoy / substitution (suppressed by default)", FAINT, True)):
        f.circle(lx + 5, ly2 - 4, 5.5, fill=PAPER if hollow else col, stroke=col if hollow else "none", sw=1.6)
        f.text(lx + 16, ly2, lab, size=13, fill=INK2, tag="leg")
        lx += 16 + text_width(lab, 13) + 28

    # --------------------------------------------------------- warning band
    wy = ly2 + 52
    f.rect(60, wy, 6, 92, fill=RED)
    f.path(f"M{92:.1f},{wy + 30:.1f} l17,-30 l17,30 Z", fill="none", stroke=RED, sw=2.4)
    f.text(109, wy + 27, "!", size=19, weight="700", fill=RED, anchor="middle", tag="warn")
    f.text(150, wy + 28, "Mass compatibility does not establish modification identity or site localization.", size=21, weight="700", fill=INK, tag="warn")
    f.text(150, wy + 56, "Several distinct entries can fit one mass window (here two biological PTMs). No peptide sequence is identified; no site is localized.", size=14, fill=INK2, tag="warn")
    f.text(150, wy + 78, "Fixed modifications with no unmodified partner spectrum can be invisible to a pair-based scout.", size=14, fill=INK2, tag="warn")
    f.h = wy + 120
    return f


if __name__ == "__main__":
    # Use the common build driver for portable output paths and layout linting.
    from build_all import main
    raise SystemExit(main(["--figures", "4"]))
