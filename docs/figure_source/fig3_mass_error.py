"""Figure 3 - identification-free mass-error precision (synthetic data; numbers from prideQC's _robust_error)."""
from __future__ import annotations
import math
import numpy as np
from figlib import *
from impl import mass_error

SIGMA_TRUE = 1.6   # single-measurement sigma (ppm) used to generate the synthetic data
rng = np.random.default_rng(21)


def build() -> Fig:
    W, H = 1400, 1368
    f = Fig(W, H, "Identification-free mass-error precision from repeated observations",
            "Synthetic repeated precursor and fragment observations give pairwise ppm errors; a robust median/MAD estimate, "
            "divided by the square root of two, gives single-measurement sigma; six sigma is the supported tolerance suggestion.")
    f.text(60, 54, "Mass-error precision without identifications", size=28, weight="700", tag="title")
    f.text(60, 82, "Synthetic data. The same ion is observed repeatedly; spread between repeats measures precision, not the search settings once used.",
           size=15, fill=MUTED, tag="title")

    # ================================================= synthetic data (real estimator for numbers)
    n_pairs = 1900
    d = rng.normal(0, math.sqrt(2) * SIGMA_TRUE, n_pairs)
    bad = rng.random(n_pairs) < 0.05
    d[bad] = rng.uniform(-20, 20, bad.sum())             # broad background matches inside the 20 ppm candidate window
    R = mass_error._robust_error(list(d))
    n_clusters = 262
    sd_naive = float(np.std(d, ddof=1))
    tol = 6 * R["single_measurement_sigma"]
    conf = "high" if n_pairs >= 1000 and n_clusters >= 250 else "moderate"

    # ================================================= A: repeated observations
    panel_label(f, 60, 138, "A", "Repeated observations of the same ion", "precursors (left) and fragment peak centres of likely repeated spectra (right)")
    # ---- A1 precursor cluster
    ax, ay, aw, ah = 118, 214, 530, 190
    xs = Scale(0, 120, ax, ax + aw)
    ys = Scale(-28, 28, ay + ah, ay)
    axes(f, ax, ay, aw, ah, xs, ys, [0, 30, 60, 90, 120], [-20, -10, 0, 10, 20], xlabel="retention time relative to first observation (s)",
         ylabel="m/z offset from cluster mean (ppm)", ylabel_dx=-52, tag="A1")
    f.rect(ax, ys(20), aw, ys(-20) - ys(20), fill=TEAL_L, opacity=0.45)
    f.line(ax, ys(20), ax + aw, ys(20), TEAL, 1.2, dash="5 4"); f.line(ax, ys(-20), ax + aw, ys(-20), TEAL, 1.2, dash="5 4")
    f.text(ax + aw - 6, ys(20) + 17, "±20 ppm candidate window", size=13, fill=TEAL, weight="600", anchor="end", tag="A1", halo=True)
    t_obs = np.array([4, 19, 33, 52, 67, 84, 99, 114.0])
    e_obs = np.round(rng.normal(0, SIGMA_TRUE, len(t_obs)), 2)
    e_obs -= e_obs.mean()
    for i in range(len(t_obs) - 1):
        f.line(xs(t_obs[i]), ys(e_obs[i]), xs(t_obs[i + 1]), ys(e_obs[i + 1]), TEAL, 1.6, opacity=0.8)
    for t, e in zip(t_obs, e_obs):
        f.circle(xs(t), ys(e), 5, fill=TEAL, stroke=PAPER, sw=1.5)
    # unrelated precursors (different charge / m/z) are not repeats
    for t, e in ((18, 25), (58, -25), (101, 24)):
        f.path(f"M{xs(t)-5:.1f},{ys(e):.1f} l5,-5 l5,5 l-5,5 Z", fill=PAPER, stroke=FAINT, sw=1.6)
    f.text(xs(18) + 12, ys(25) + 5, "different charge or m/z: not a repeat", size=13, fill=MUTED, tag="A1", halo=True)
    # annotate consecutive delta
    k = 3
    f.line(xs(t_obs[k]) + 13, ys(e_obs[k]), xs(t_obs[k]) + 13, ys(e_obs[k + 1]), INK, 1.4)
    f.text(xs(t_obs[k]) + 22, ys(10), "pairwise error (consecutive repeats)", size=13, weight="600", tag="A1", halo=True)
    f.line(xs(t_obs[k]) + 13, ys(e_obs[k]) - 6, xs(t_obs[k]) + 13, ys(10) + 4, INK, 1, dash="2 3")
    # constraint chips
    cx0, cw = chip(f, ax, ay - 14, "same charge", fg=TEAL, bg=TEAL_L, tag="A1c")
    cx1, cw1 = chip(f, ax + cw + 8, ay - 14, "≤ 120 s apart", fg=TEAL, bg=TEAL_L, tag="A1c")
    chip(f, ax + cw + 8 + cw1 + 8, ay - 14, "≤ 20 ppm of cluster mean", fg=TEAL, bg=TEAL_L, tag="A1c")

    # ---- A2 fragment stems + zoom inset
    fx0, fw = 760, 330
    fy_mid = 308
    mz = np.sort(rng.uniform(300, 1100, 34))
    inten = rng.lognormal(0, 0.9, 34); inten /= inten.max()
    ex = Scale(300, 1100, fx0, fx0 + fw)
    f.line(fx0, fy_mid, fx0 + fw, fy_mid, INK2, 1.2)
    for m, i in zip(mz, inten):
        f.line(ex(m), fy_mid - 2, ex(m), fy_mid - 2 - i * 78, BLUE, 1.8)
        f.line(ex(m), fy_mid + 2, ex(m), fy_mid + 2 + i * 78 * rng.uniform(0.7, 1.2), TEAL, 1.8)
    for m in (320, 560, 800, 1040):
        f.line(ex(m), fy_mid, ex(m), fy_mid + 5, INK2, 1)
        f.text(ex(m), fy_mid + 20 + 0, "", size=1, check=False)
    f.text(fx0 + fw / 2, ay + ah + 38 - 4, "m/z", size=13.5, fill=INK2, anchor="middle", tag="A2")
    for m in (400, 600, 800, 1000):
        f.line(ex(m), ay + ah, ex(m), ay + ah + 5, INK2, 1)
        f.text(ex(m), ay + ah + 20, f"{m}", size=12.5, fill=INK2, anchor="middle", tag="A2")
    f.line(fx0, ay + ah, fx0 + fw, ay + ah, INK2, 1.2)
    f.text(fx0, ay - 6, "spectrum 1", size=13, weight="600", fill=BLUE, tag="A2")
    f.text(fx0 + fw, ay - 6, "spectrum 2 (repeat)", size=13, weight="600", fill=TEAL, anchor="end", tag="A2")
    f.text(fx0 + fw / 2, ay + ah - 4 + 0, "", size=1, check=False)
    # top-50, one-to-one match
    f.text(fx0, ay + ah + 58, "top-50 peak centres, matched one-to-one (≤ 0.2 Da)", size=13, fill=MUTED, tag="A2")
    # zoom inset on one matched pair
    zx, zw = 1130, 210
    zxs = Scale(-0.0075, 0.0075, zx, zx + zw)
    zy = Scale(0, 1.1, ay + ah - 10, ay + 10)
    z_c = 645.3612
    g = np.linspace(-0.0075, 0.0075, 120)
    s_ppm = 2.6
    off = z_c * s_ppm * 1e-6 / 1.0
    sig_pk = 0.0021
    f.path(smooth_path([(zxs(a), zy(math.exp(-0.5 * (a / sig_pk) ** 2))) for a in g], closed_to=zy(0)), fill=BLUE_L, stroke="none", opacity=0.9)
    f.poly([(zxs(a), zy(math.exp(-0.5 * (a / sig_pk) ** 2))) for a in g], stroke=BLUE, sw=1.8)
    f.poly([(zxs(a), zy(0.8 * math.exp(-0.5 * ((a - off) / sig_pk) ** 2))) for a in g], stroke=TEAL, sw=1.8)
    f.line(zxs(0), zy(0), zxs(0), zy(1.02), BLUE, 1, dash="3 3"); f.line(zxs(off), zy(0), zxs(off), zy(0.84), TEAL, 1, dash="3 3")
    f.line(zx, zy(0), zx + zw, zy(0), INK2, 1.2)
    f.line(zxs(0), zy(1.07), zxs(off), zy(1.07), INK, 1.5)
    f.text(zx + zw / 2, ay - 6, f"zoom: one pair at m/z {z_c:.2f}", size=13, weight="600", fill=INK2, anchor="middle", tag="A2")
    f.text(zxs(off) + 8, zy(1.07) + 4, f"+{s_ppm:.1f} ppm", size=13.5, weight="600", tag="A2")
    for v in (-0.005, 0, 0.005):
        f.line(zxs(v), zy(0), zxs(v), zy(0) + 5, INK2, 1); f.text(zxs(v), zy(0) + 20, f"{v:+.3f}" if v else "0", size=12.5, fill=INK2, anchor="middle", tag="A2")
    f.text(zx + zw / 2, ay + ah + 40, "Δm/z (Da)", size=13, fill=INK2, anchor="middle", tag="A2")
    # formula
    f.rich(ax, ay + ah + 82, [("error (ppm) = (m₂ − m₁) / ((m₁ + m₂)/2) × 10⁶", {"mono": True})], size=14, tag="A")
    f.text(fx0 + 0, ay + ah + 78, "Profile spectra: ephemeral peak centres. Input spectra are never modified.", size=13, fill=MUTED, tag="A")

    # ================================================= B, C, D
    top = 562
    f.line(60, top - 52, W - 60, top - 52, RULE, 1)
    # ---- B histogram of pairwise errors
    bx0, bw, bh = 118, 330, 200
    panel_label(f, 60, top - 14, "B", "Pairwise error distribution", f"{n_pairs:,} precursor repeat pairs from {n_clusters} clusters")
    lim = 14
    xs = Scale(-lim, lim, bx0, bx0 + bw)
    edges = np.linspace(-lim, lim, 57)
    inside = d[(d > -lim) & (d < lim)]
    cnt, _ = np.histogram(inside, bins=edges)
    ys = Scale(0, cnt.max() * 1.12, top + 50 + bh, top + 50)
    axes(f, bx0, top + 50, bw, bh, xs, ys, [-12, -6, 0, 6, 12], [0, 100, 200], xlabel="pairwise error (ppm)", ylabel="pairs per bin", ylabel_dx=-46, tag="B")
    c0 = R["pairwise_median"]; th = R["robust_inlier_threshold_3sigma"]
    f.rect(xs(c0 - th), top + 50, xs(c0 + th) - xs(c0 - th), bh, fill=AMBER_L, opacity=0.55)
    for c, e0, e1 in zip(cnt, edges[:-1], edges[1:]):
        f.rect(xs(e0) + 0.5, ys(c), xs(e1) - xs(e0) - 1, ys(0) - ys(c), fill=BLUE)
    f.line(xs(c0), top + 44, xs(c0), top + 50 + bh, INK, 1.5, dash="5 4")
    n_out_far = int(np.sum(np.abs(d) >= lim))
    f.text(bx0 + bw, top + 38, "", size=1, check=False)
    f.text(xs(c0 - th) - 6, top + 66, "inlier band", size=13, weight="600", fill=AMBER, anchor="end", tag="B")
    f.text(xs(c0 - th) - 6, top + 83, "±3 pairwise σ", size=13, fill=AMBER, anchor="end", tag="B")
    f.text(bx0, top + 50 + bh + 56, f"{R['robust_inlier_fraction']*100:.1f}% inliers; {R['robust_outlier_count']} broad background pairs kept in the distribution",
           size=13, fill=INK2, tag="B")
    # ---- C robust estimate (CDF of |e - center|)
    cx0_, cw_ = 520, 330
    panel_label(f, 462, top - 14, "C", "Robust scale: median and MAD", "robust to background matches")
    a = np.sort(np.abs(d - c0))
    xmax = 8
    xs2 = Scale(0, xmax, cx0_, cx0_ + cw_)
    ys2 = Scale(0, 1, top + 50 + bh, top + 50)
    axes(f, cx0_, top + 50, cw_, bh, xs2, ys2, [0, 2, 4, 6, 8], [0, 0.5, 1], xlabel="|error − center| (ppm)", ylabel="cumulative fraction", ylabel_dx=-46, tag="C")
    cdf_pts = [(xs2(min(v, xmax)), ys2((i + 1) / len(a))) for i, v in enumerate(a) if v <= xmax][::6]
    f.poly(cdf_pts, stroke=BLUE, sw=2.4)
    mad = float(np.median(a))
    f.line(cx0_, ys2(0.5), xs2(mad), ys2(0.5), AMBER, 1.6, dash="5 4")
    f.line(xs2(mad), ys2(0.5), xs2(mad), ys2(0), AMBER, 1.6, dash="5 4")
    f.circle(xs2(mad), ys2(0.5), 5, fill=AMBER)
    f.text(xs2(mad) + 12, ys2(0.5) + 20, f"MAD = {mad:.2f} ppm", size=14, weight="600", fill=AMBER, tag="C", halo=True)
    # sample SD marker for contrast
    f.line(xs2(min(sd_naive, xmax)), ys2(0), xs2(min(sd_naive, xmax)), ys2(0.12), SLATE, 1.6)
    f.text(xs2(min(sd_naive, xmax)) + 7, ys2(0.12) - 4, f"plain SD = {sd_naive:.2f}", size=13, fill=SLATE, tag="C", halo=True)
    # chain of numbers (not boxes): rule-separated rows
    ry = top + 50 + bh + 62
    # ---- D tolerance
    dx0, dw = 920, 330
    panel_label(f, 862, top - 14, "D", "Supported tolerance suggestion", "single-measurement precision × 6")
    sg = R["single_measurement_sigma"]
    xs3 = Scale(-1.1 * tol, 1.1 * tol, dx0, dx0 + dw)
    gx = np.linspace(-1.1 * tol, 1.1 * tol, 200)
    dens = np.exp(-0.5 * (gx / sg) ** 2)
    ys3 = Scale(0, 1.38, top + 50 + bh, top + 50)
    axes(f, dx0, top + 50, dw, bh, xs3, ys3, [-tol, 0, tol], [], xlabel="error of one measurement (ppm)", grid=False, tag="D", xfmt=lambda v: f"{v:+.1f}" if abs(v) > 1e-9 else "0")
    f.rect(xs3(-tol), top + 50, xs3(tol) - xs3(-tol), bh, fill=GREEN_L, opacity=0.55)
    f.path(smooth_path([(xs3(a_), ys3(b_)) for a_, b_ in zip(gx, dens)], closed_to=ys3(0)), fill=BLUE_L, stroke="none")
    f.poly([(xs3(a_), ys3(b_)) for a_, b_ in zip(gx, dens)], stroke=BLUE, sw=2.2)
    f.line(xs3(-tol), top + 50, xs3(-tol), top + 50 + bh, GREEN, 2); f.line(xs3(tol), top + 50, xs3(tol), top + 50 + bh, GREEN, 2)
    f.line(xs3(sg), ys3(math.exp(-0.5)), xs3(sg), ys3(0), BLUE, 1.5, dash="3 3")
    f.text(xs3(sg) + 7, ys3(0.30) , f"σ = {sg:.2f}", size=13.5, weight="600", fill=BLUE, tag="D", halo=True)
    f.line(xs3(-tol), top + 50 + 14, xs3(tol), top + 50 + 14, GREEN, 2)
    f.text(dx0 + dw / 2, top + 50 + 38, f"suggested ±{tol:.1f} ppm = 6σ", size=14.5, weight="700", fill=GREEN, anchor="middle", tag="D")
    f.text(dx0 + dw / 2, ys3(0) - 8, "", size=1, check=False)
    f.text(dx0, ry, "Emitted only when support gates pass:", size=13.5, weight="600", tag="D")
    # gates under D
    g1 = [("≥ 200 pairs", n_pairs >= 200, f"{n_pairs:,}"), ("≥ 100 clusters", n_clusters >= 100, f"{n_clusters}")]

    # numeric chain across C bottom (rule-separated lines)
    y0 = ry + 24
    rows = [("center = median(e)", f"{R['pairwise_median']:+.2f} ppm"),
            ("MAD = median(|e − center|)", f"{mad:.2f} ppm"),
            ("pairwise σ = 1.4826 × MAD", f"{R['pairwise_sigma']:.2f} ppm"),
            ("single σ = pairwise σ / √2", f"{sg:.2f} ppm")]
    for i, (lab, val) in enumerate(rows):
        yy = y0 + i * 24 - 24
        if i == 0:
            continue
    # place numeric chain in column C under the rich line (replace the first row with a clean 4-row table)
    # (first row already printed above; print the remaining ones)
    for i, (lab, val) in enumerate(rows):
        f.text(cx0_, ry + i * 22, lab, size=13.5, mono=True, tag="Cr")
        f.text(cx0_ + cw_, ry + i * 22, val, size=13.5, mono=True, weight="700", anchor="end", tag="Cr", fill=AMBER if i == 3 else INK)
    # gate rows under D
    for i, (lab, ok, val) in enumerate(g1):
        yy = ry + 24 + i * 22
        f.circle(dx0 + 8, yy - 5, 6, fill=GREEN)
        f.path(f"M{dx0+4.5:.1f},{yy-5:.1f} l2.5,2.8 l4.5,-5.4", stroke=PAPER, sw=1.8)
        f.text(dx0 + 24, yy, lab, size=13.5, mono=True, tag="Dg")
        f.text(dx0 + dw, yy, val, size=13.5, mono=True, weight="700", anchor="end", tag="Dg")
    yy = ry + 24 + 2 * 22
    f.text(dx0, yy, f"6 × {sg:.2f} = {tol:.1f} ppm   (confidence: {conf}: ≥ 1000 pairs, ≥ 250 clusters)", size=13, mono=False, fill=INK2, tag="Dg")

    # ================================================= E: regimes, abstention, output units
    ey0 = 1062
    f.line(60, ey0 - 62, W - 60, ey0 - 62, RULE, 1)
    panel_label(f, 60, ey0 - 26, "E", "Fragment precision: resolution regime decides the unit, ambiguity abstains",
                "decision on the narrow 0.2 Da evidence; precursor tolerance is reported in ppm only")
    ax0, aw0 = 118, 640
    xs = Scale(1, 100, ax0, ax0 + aw0, log=True)
    base = ey0 + 100
    f.rect(xs(1), base - 40, xs(10) - xs(1), 40, fill=TEAL_L)
    f.rect(xs(20), base - 40, xs(100) - xs(20), 40, fill=PLUM_L)
    # abstain gap hatched
    f.defs.append('<pattern id="hatch" width="7" height="7" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="7" stroke="#B03A2E" stroke-width="2" opacity="0.55"/></pattern>')
    f.rect(xs(10), base - 40, xs(20) - xs(10), 40, fill="url(#hatch)")
    f.line(ax0, base, ax0 + aw0, base, INK2, 1.2)
    for v in (1, 2, 5, 10, 20, 50, 100):
        f.line(xs(v), base, xs(v), base + 5, INK2, 1); f.text(xs(v), base + 22, f"{v}", size=12.5, fill=INK2, anchor="middle", tag="E")
    f.text(ax0 + aw0 / 2, base + 44, "single-measurement fragment σ from 0.2 Da pairs (ppm)", size=13.5, fill=INK2, anchor="middle", tag="E")
    f.text(xs(3.1), base - 54, "high resolution", size=14, weight="700", fill=TEAL, anchor="middle", tag="E")
    f.text(xs(3.1), base - 72, "σ ≤ 10 ppm and ≤ 0.01 Da", size=13, fill=INK2, anchor="middle", tag="E")
    f.text(xs(14.1), base - 54, "abstain", size=14, weight="700", fill=RED, anchor="middle", tag="E")
    f.text(xs(14.1), base - 72, "intermediate", size=13, fill=INK2, anchor="middle", tag="E")
    f.text(xs(45), base - 54, "low resolution", size=14, weight="700", fill=PLUM, anchor="middle", tag="E")
    f.text(xs(45), base - 72, "σ ≥ 20 ppm and ≥ 0.01 Da", size=13, fill=INK2, anchor="middle", tag="E")
    f.text(xs(3.1), base - 17, "tolerance in ppm", size=13.5, weight="600", fill=TEAL, anchor="middle", tag="E")
    f.text(xs(45), base - 17, "tolerance in Da", size=13.5, weight="600", fill=PLUM, anchor="middle", tag="E")
    # (ppm,Da) discordant note
    f.text(ax0, base + 78, "Discordant ppm/Da evidence also abstains. Both conditions must hold for either regime.", size=13, fill=MUTED, tag="E")

    # censoring check
    cx2 = 800
    f.text(cx2, ey0 + 18, "Low-resolution censoring check", size=14.5, weight="600", tag="E2")
    wxs = Scale(0, 1.0, cx2, cx2 + 200)
    rows_w = [(0.5, 0.47, "0.5 Da window", "3σ core 0.47 ≥ 0.45 (90%): censored, widen", RED),
              (1.0, 0.62, "1.0 Da window", "3σ core 0.62 < 0.90: uncensored, use it", GREEN)]
    for i, (win, core, lab, note, col) in enumerate(rows_w):
        y_ = ey0 + 50 + i * 66
        f.rect(wxs(0), y_, wxs(win) - wxs(0), 16, fill=WASH, stroke=FAINT, sw=1.2)
        f.rect(wxs(0), y_ + 3, wxs(core) - wxs(0), 10, fill=col)
        f.line(wxs(0.9 * win), y_ - 4, wxs(0.9 * win), y_ + 20, INK, 1.5)
        f.text(wxs(win) + 8, y_ + 13, lab if i == 1 else "", size=13.5, weight="600", tag="E2")
        f.text(cx2, y_ + 38, note, size=13, fill=col, weight="600", tag="E2")
    f.text(cx2 + 108, ey0 + 50 + 13, "0.5 Da", size=13.5, weight="600", tag="E2")
    f.text(cx2, ey0 + 50 + 2 * 66 - 6, "bar: robust 3σ core · tick: 90% of window", size=13, fill=MUTED, tag="E2")
    f.text(cx2, ey0 + 50 + 2 * 66 + 14, "All windows censored: the Da suggestion abstains.", size=13, fill=MUTED, tag="E2")

    # native high-res mixture illustration
    mx0, mw = 1110, 230
    f.text(mx0, ey0 + 18, "High-resolution fragment fit", size=14.5, weight="600", tag="E3")
    gx = np.linspace(-20, 20, 160)
    sig_f = 2.0
    bgd = 0.10
    dens_f = (1 - bgd) * np.exp(-0.5 * (gx / sig_f) ** 2) / (sig_f * math.sqrt(2 * math.pi)) + bgd / 40
    mxs = Scale(-20, 20, mx0, mx0 + mw)
    mys = Scale(0, dens_f.max() * 1.1, ey0 + 140, ey0 + 52)
    f.path(smooth_path([(mxs(a_), mys(b_)) for a_, b_ in zip(gx, dens_f)], closed_to=mys(0)), fill=BLUE_L, stroke="none")
    f.poly([(mxs(a_), mys(b_)) for a_, b_ in zip(gx, dens_f)], stroke=BLUE, sw=2)
    f.line(mx0, mys(bgd / 40), mx0 + mw, mys(bgd / 40), RED, 1.6, dash="5 4")
    f.line(mx0, mys(0), mx0 + mw, mys(0), INK2, 1.2)
    for v in (-20, 0, 20):
        f.line(mxs(v), mys(0), mxs(v), mys(0) + 5, INK2, 1); f.text(mxs(v), mys(0) + 20, f"{v:+d}" if v else "0", size=12.5, fill=INK2, anchor="middle", tag="E3")
    f.text(mx0 + mw / 2, mys(0) + 40, "coarse-bin pair error (ppm)", size=13, fill=INK2, anchor="middle", tag="E3")
    f.text(mx0 + mw, mys(bgd / 40) - 22, "uniform background", size=13, fill=RED, anchor="end", tag="E3", halo=True)
    f.text(mxs(3.2), mys(dens_f.max()) + 4, "Gaussian signal", size=13, fill=BLUE, weight="600", tag="E3")
    f.text(mx0, ey0 + 212, "Mixture fit gives pairwise σ,", size=13, fill=MUTED, tag="E3")
    f.text(mx0, ey0 + 230, "then ÷ √2 and × 6 as in C and D.", size=13, fill=MUTED, tag="E3")

    # final strip: output line
    f.line(60, 1308, W - 60, 1308, RULE, 1)
    f.text(60, 1338, "This is a precision estimate: it cannot see a fixed calibration offset or an isotope-error policy, and it is not the historical search tolerance.",
           size=14.5, fill=INK, weight="600", tag="foot")
    return f


if __name__ == "__main__":
    # Use the common build driver for portable output paths and layout linting.
    from build_all import main
    raise SystemExit(main(["--figures", "3"]))
