"""Figure 2 - how prideQC reduces a synthetic LC-MS run to compact QC metrics.
Every number printed in the figure is produced by prideQC's own metrics.py on synthetic scans."""
from __future__ import annotations
import math
import numpy as np
from figlib import *
import synth_run

SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def sci(v: float, d=1) -> str:
    e = int(math.floor(math.log10(abs(v))))
    return f"{v / 10**e:.{d}f}×10{str(e).translate(SUP)}"


def build() -> Fig:
    D = synth_run.build()
    M = D["metrics"]
    fixed_M, fixed_grid = synth_run.fixed_target_metrics()
    W, H = 1400, 1486
    f = Fig(W, H, "How prideQC computes QC metrics from a spectrum stream",
            "A synthetic 90 minute DDA LC-MS/MS run is reduced to TIC, base-peak, scan-rate, charge, isolation-window and "
            "repeated-target metrics. Formulas are those used in prideQC's metric code.")

    # ---------------------------------------------------------------- header
    f.text(60, 54, "From a spectrum stream to compact QC metrics", size=28, weight="700", tag="title")
    f.text(60, 82, "Synthetic 90-minute DDA LC-MS/MS run. Every value shown is computed by prideQC's own metric code on the synthetic scans.",
           size=15, fill=MUTED, tag="title")

    # ---------------------------------------------------------------- A: stream
    panel_label(f, 60, 138, "A", "Raw stream, reduced scan by scan",
                f"{M['NumberOfSpectra_MS1']:,} MS1 and {M['NumberOfSpectra_MS2']:,} MS2 scans. For these metrics each scan is reduced to a few scalars; peak arrays are not retained.")
    px, pw = 100, 1240
    top, h = 192, 108
    base = top + h
    rt_max = 5400.0
    xs = Scale(0, rt_max, px, px + pw)
    ys = Scale(0, 2.7e9, base, top)
    # TIC curve, binned to 20 s medians for a clean trace
    bins = np.arange(0, rt_max + 1, 20)
    idx = np.digitize(D["ms1_rt"], bins) - 1
    tic_b = np.array([np.median(D["ms1_tic"][idx == i]) if np.any(idx == i) else np.nan for i in range(len(bins) - 1)])
    xc = (bins[:-1] + bins[1:]) / 2
    ok = ~np.isnan(tic_b)
    pts = [(xs(a), ys(min(b, 2.7e9))) for a, b in zip(xc[ok], tic_b[ok])]
    axes(f, px, top, pw, h, xs, ys, [0, 900, 1800, 2700, 3600, 4500, 5400], [0, 1e9, 2e9],
         xlabel="", ylabel="MS1 TIC (a.u.)", xfmt=lambda v: f"{v/60:g}", yfmt=lambda v: "0" if v == 0 else sci(v, 0), ylabel_dx=-66, tag="A")
    f.text(px + pw, base + 40, "retention time (min)", size=13.5, fill=INK2, anchor="end", tag="A")
    f.path(smooth_path(pts, closed_to=base), fill=BLUE_L, stroke="none", opacity=0.9)
    f.poly(pts, stroke=BLUE, sw=1.6)
    # zoom window
    z0, z1 = 2460.0, 2480.0
    zc = (xs(z0) + xs(z1)) / 2
    f.rect(zc - 3.5, top, 7, h, fill=AMBER, opacity=0.55)
    # wedge to the raster
    by = base + 52
    ry_top = by + 34
    rb = ry_top + 142
    f.line(zc - 3.5, base + 6, px, ry_top, AMBER, 1.3, dash="5 4"); f.line(zc + 3.5, base + 6, px + pw, ry_top, AMBER, 1.3, dash="5 4")
    f.rect(px, ry_top, pw, rb - ry_top, fill="none", stroke=AMBER, sw=1.3)
    f.text(zc + 10, top - 6, "20 s zoom (below)", size=13.5, fill=AMBER, weight="600", tag="A")

    # raster of scans inside the window
    rxs = Scale(z0, z1, px, px + pw)
    f.line(px, rb, px + pw, rb, INK2, 1.2)
    for v in range(int(z0), int(z1) + 1, 5):
        f.line(rxs(v), rb, rxs(v), rb + 5, INK2, 1)
        f.text(rxs(v), rb + 22, f"{v/60:.0f}:{v%60:02.0f}", size=12.5, fill=INK2, anchor="middle", tag="A")
    f.text(px + pw / 2, rb + 44, "retention time (min:s)", size=13.5, fill=INK2, anchor="middle", tag="A")
    CHC = {1: "#8A97A8", 2: TEAL, 3: BLUE, 4: AMBER, 5: PLUM, 6: RED, 0: "#BCC6D2"}
    sel1 = (D["ms1_rt"] >= z0 - 2.5) & (D["ms1_rt"] <= z1 + 0.01)
    ms1_rts = D["ms1_rt"][sel1]
    ms1_in = ms1_rts[(ms1_rts >= z0) & (ms1_rts <= z1)]
    for t, tic in zip(D["ms1_rt"][sel1], D["ms1_tic"][sel1]):
        if z0 <= t <= z1:
            f.rect(rxs(t) - 2.2, rb - 74, 4.4, 74, fill=BLUE)
    sel2 = (D["ms2_rt"] >= z0) & (D["ms2_rt"] <= z1)
    rng = np.random.default_rng(3)
    for t, z in zip(D["ms2_rt"][sel2], D["ms2_charge"][sel2]):
        hh = 10 + 26 * rng.random() ** 1.7
        f.rect(rxs(t) - 1.7, rb - hh, 3.4, hh, fill=CHC[min(z, 6)])
    # cycle-time dimension between two consecutive MS1 scans
    a, b = ms1_in[3], ms1_in[4]
    dy = rb - 74 - 14
    f.line(rxs(a), dy, rxs(b), dy, INK, 1.5)
    f.line(rxs(a), dy - 5, rxs(a), dy + 5, INK, 1.5); f.line(rxs(b), dy - 5, rxs(b), dy + 5, INK, 1.5)
    f.text((rxs(a) + rxs(b)) / 2, dy - 10, f"cycle time ≈ {b - a:.1f} s  (ΔRT between consecutive MS1 scans)", size=13.5, fill=INK, anchor="middle", weight="600", tag="A")
    # legend (left, under the axis label)
    ly = rb + 76
    f.rect(px, ly - 9, 5, 18, fill=BLUE); f.text(px + 14, ly + 5, "MS1 scan", size=13.5, fill=INK2, tag="A")
    lx = px + 120
    f.text(lx, ly + 5, "MS2 scan, coloured by precursor charge:", size=13.5, fill=INK2, tag="A")
    lx += 292
    for lab, z in (("1+", 1), ("2+", 2), ("3+", 3), ("4+", 4), ("5+", 5), ("≥6", 6), ("unknown", 0)):
        f.rect(lx, ly - 9, 12, 18, fill=CHC[z]); f.text(lx + 18, ly + 5, lab, size=13.5, fill=INK2, tag="A")
        lx += 18 + text_width(lab, 13.5) + 18
    f.text(px + pw, ly + 5, "MS2 bar height: relative TIC (illustrative)", size=13.5, fill=MUTED, anchor="end", tag="A")

    # ---------------------------------------------------------------- metric families
    f.line(60, 672, W - 60, 672, RULE, 1)
    f.text(60, 712, "Metric families derived from the stream", size=20, weight="700", tag="mf")
    cols = [60, 500, 940]
    pw_ = 330
    row_y = [770, 1110]
    ph = 130

    def callout(x, y, rows, color):
        f.rect(x, y - 15, 3.5, 17 * len(rows) + 8, fill=color)
        for i, r in enumerate(rows):
            f.text(x + 14, y + i * 20, r.replace("  ", "\u00a0\u00a0"), size=13.2, mono=True, fill=INK, tag="co")

    # ---- B: TIC area & CV (trapezoid geometry on 8 real synthetic scans)
    cx, cy = cols[0], row_y[0]
    panel_label(f, cx, cy, "B", "TIC, area and variability", "MS1 total ion current per scan")
    px2, top2 = cx + 70, cy + 44
    scans = [(t, v) for t, v in zip(D["ms1_rt"], D["ms1_tic"]) if 2461 <= t <= 2480][:8]
    tt = np.array([s[0] for s in scans]); vv = np.array([s[1] for s in scans])
    sx = Scale(tt[0], tt[-1], px2 + 10, px2 + pw_ - 20)
    sy = Scale(0, vv.max() * 1.25, top2 + ph, top2)
    axes(f, px2, top2, pw_ - 14, ph, sx, sy, [], [0, round(vv.max() / 1e9, 1) * 1e9], ylabel="TIC (a.u.)", yfmt=lambda v: "0" if v == 0 else sci(v, 1), ylabel_dx=-60, tag="B")
    for i in range(len(tt) - 1):
        poly = [(sx(tt[i]), sy(0)), (sx(tt[i]), sy(vv[i])), (sx(tt[i + 1]), sy(vv[i + 1])), (sx(tt[i + 1]), sy(0))]
        f.poly(poly, fill="#B7CFE6" if i == 3 else (BLUE_L if i % 2 == 0 else "#E8F0F8"), stroke="none", close=True)
    f.poly([(sx(a), sy(b)) for a, b in zip(tt, vv)], stroke=BLUE, sw=2)
    for a, b in zip(tt, vv):
        f.circle(sx(a), sy(b), 3.2, fill=BLUE)
    # annotate trapezoid k
    k = 3
    f.line(sx(tt[k]), sy(0) + 12, sx(tt[k + 1]), sy(0) + 12, INK, 1.4)
    f.line(sx(tt[k]), sy(0) + 7, sx(tt[k]), sy(0) + 17, INK, 1.4); f.line(sx(tt[k + 1]), sy(0) + 7, sx(tt[k + 1]), sy(0) + 17, INK, 1.4)
    f.text(sx(tt[k + 1]) + 10, sy(0) + 17, "ΔRT", size=13.5, anchor="start", weight="600", tag="B")
    f.text(sx(tt[k]) - 4, sy(vv[k]) - 12, "TIC_i", size=13.5, anchor="end", weight="600", fill=BLUE, tag="B", halo=True)
    f.text(sx(tt[k + 1]) + 4, sy(vv[k + 1]) - 12, "TIC_i+1", size=13.5, anchor="start", weight="600", fill=BLUE, tag="B", halo=True)
    f.text(px2 + (pw_ - 14) / 2, sy(0) + 40, "retention time (8 consecutive MS1 scans)", size=13.5, fill=INK2, anchor="middle", tag="B")
    callout(cx, cy + ph + 112, [
        "TIC area = Σ ΔRT × (TIC_i + TIC_i+1) / 2",
        f"  run total {sci(M['TIC_MS1_Area'])} a.u.·s",
        "CV = sample SD (n−1) / mean",
        f"  TIC_MS1_CV = {M['TIC_MS1_CV']:.2f}"], BLUE)

    # ---- C: base peak
    cx, cy = cols[1], row_y[0]
    panel_label(f, cx, cy, "C", "Base-peak intensity", "most intense peak of each MS1 scan")
    px2, top2 = cx + 70, cy + 44
    bx = Scale(0, rt_max, px2, px2 + pw_ - 14)
    bp_top = math.ceil(D["ms1_bp"].max() / 5e7) * 5e7
    by_ = Scale(0, bp_top, top2 + ph, top2)
    axes(f, px2, top2, pw_ - 14, ph, bx, by_, [0, 1800, 3600, 5400], [0, bp_top / 2, bp_top], xlabel="retention time (min)", ylabel="base peak (a.u.)",
         xfmt=lambda v: f"{v/60:g}", yfmt=lambda v: "0" if v == 0 else sci(v, 0), ylabel_dx=-60, tag="C")
    for t, v in list(zip(D["ms1_rt"], D["ms1_bp"]))[::2]:
        f.circle(bx(t), by_(v), 1.35, fill=BLUE, opacity=0.35)
    mean_bp = M["BasePeak_MS1_Mean"]
    f.line(px2, by_(mean_bp), px2 + pw_ - 14, by_(mean_bp), AMBER, 2, dash="6 4")
    f.text(px2 + pw_ - 14, by_(mean_bp) - 8, f"mean {sci(mean_bp)}", size=13.5, weight="600", fill=AMBER, anchor="end", tag="C", halo=True)
    mx_i = int(np.argmax(D["ms1_bp"]))
    f.circle(bx(D["ms1_rt"][mx_i]), by_(D["ms1_bp"][mx_i]), 4.5, fill=PAPER, stroke=INK, sw=1.6)
    f.text(bx(D["ms1_rt"][mx_i]) - 9, by_(D["ms1_bp"][mx_i]) + 5, f"max {sci(M['BasePeak_All_Max'])}", size=13, fill=INK, anchor="end", tag="C", halo=True)
    callout(cx, cy + ph + 112, [
        "base peak = max(intensity) per scan",
        f"  BasePeak_MS1_Mean = {sci(mean_bp)}",
        "BasePeak_All_Max = max over all levels",
        f"  = {sci(M['BasePeak_All_Max'])} a.u."], BLUE)

    # ---- D: scan frequency & cycle time
    cx, cy = cols[2], row_y[0]
    panel_label(f, cx, cy, "D", "Scan frequency and cycle time", "how fast the instrument sampled the run")
    px2, top2 = cx + 70, cy + 40
    rt2 = D["ms2_rt"]
    grid_t = np.arange(0, rt_max - 60, 10.0)
    freq = np.array([(np.searchsorted(rt2, t + 60, side="right") - np.searchsorted(rt2, t, side="left")) / 60 for t in grid_t])
    fx = Scale(0, rt_max, px2, px2 + pw_ - 14)
    fy = Scale(0, 6, top2 + 58, top2)
    axes(f, px2, top2, pw_ - 14, 58, fx, fy, [0, 1800, 3600, 5400], [0, 3, 6], ylabel="MS2 / s", xfmt=lambda v: f"{v/60:g}", xticklabels=False, ylabel_dx=-36, tag="D1")
    f.poly([(fx(a), fy(b)) for a, b in zip(grid_t, freq)], stroke=TEAL, sw=1.6)
    imax = int(np.argmax(freq))
    f.circle(fx(grid_t[imax]), fy(freq[imax]), 4.5, fill=AMBER)
    f.text(fx(grid_t[imax]) + (9 if fx(grid_t[imax]) < px2 + 150 else -9), fy(freq[imax]) + (18 if fy(freq[imax]) < top2 + 24 else -9), f"fastest 60 s window: {M['FastestFrequency_MS2']:.1f} Hz", size=13, weight="600", fill=AMBER, tag="D1", anchor="start" if fx(grid_t[imax]) < px2 + 150 else "end", halo=True)
    # cycle-time histogram
    cyc = np.diff(D["ms1_rt"])
    h0 = top2 + 58 + 36
    hx = Scale(1.2, 3.0, px2, px2 + pw_ - 14)
    counts, edges = np.histogram(cyc[(cyc > 1.2) & (cyc < 3.0)], bins=np.linspace(1.2, 3.0, 31))
    hy = Scale(0, counts.max() * 1.1, h0 + 58, h0)
    axes(f, px2, h0, pw_ - 14, 58, hx, hy, [1.2, 1.8, 2.4, 3.0], [], xlabel="MS1-to-MS1 cycle time (s)", xfmt=lambda v: f"{v:.1f}", grid=False, xlabel_dy=36, tag="D2")
    for c, e0, e1 in zip(counts, edges[:-1], edges[1:]):
        f.rect(hx(e0) + 0.8, hy(c), hx(e1) - hx(e0) - 1.6, hy(0) - hy(c), fill=TEAL_L, stroke=TEAL, sw=1)
    f.line(hx(M["AvgCycleTime_MS1"]), h0 - 4, hx(M["AvgCycleTime_MS1"]), h0 + 58, AMBER, 2, dash="6 4")
    f.text(hx(M["AvgCycleTime_MS1"]) + 7, h0 + 8, f"mean {M['AvgCycleTime_MS1']:.2f} s", size=13, weight="600", fill=AMBER, tag="D2", halo=True)
    callout(cx, cy + ph + 112, [
        "ScanRate = n scans / RT span × 60",
        f"  MS1 {M['ScanRate_MS1']:.1f} scans/min",
        "FastestFrequency = max scans in 60 s / 60",
        f"  MS2 {M['FastestFrequency_MS2']:.1f} Hz"], TEAL)

    # ---- E: charge
    cx, cy = cols[0], row_y[1]
    panel_label(f, cx, cy, "E", "Precursor charge distribution", "first precursor of each MS2 scan")
    px2, top2 = cx + 70, cy + 44
    fr = M["MS2_PrecursorCharge_Fractions"]
    labs = ["1+", "2+", "3+", "4+", "5+", "≥6", "n.a."]
    ex = Scale(0, 7, px2, px2 + pw_ - 14)
    ey = Scale(0, 0.6, top2 + ph, top2)
    axes(f, px2, top2, pw_ - 14, ph, ex, ey, [], [0, 0.2, 0.4, 0.6], ylabel="fraction of MS2", yfmt=lambda v: f"{v:g}", ylabel_dx=-46, tag="E")
    bw = (ex(1) - ex(0)) * 0.62
    for i, (lab, v) in enumerate(zip(labs, fr["fraction"])):
        z = [1, 2, 3, 4, 5, 6, 0][i]
        xc_ = ex(i + 0.5)
        f.rect(xc_ - bw / 2, ey(v), bw, ey(0) - ey(v), fill=CHC[z])
        f.text(xc_, ey(v) - 7, f"{v*100:.0f}%" if v >= 0.015 else f"{v*100:.1f}%", size=13, weight="600", anchor="middle", fill=INK, tag="E")
        f.text(xc_, ey(0) + 20, lab, size=13.5, anchor="middle", fill=INK2, tag="E")
    f.text(px2 + (pw_ - 14) / 2, ey(0) + 40, "charge state (n.a. = not assigned)", size=13.5, anchor="middle", fill=INK2, tag="E")
    callout(cx, cy + ph + 112, [
        "ChargeRatio_3over2 = n(3+) / n(2+)",
        f"  = {M['ChargeRatio_3over2']:.2f}",
        "unknown = charge ≤ 0 in the file"], TEAL)

    # ---- F: isolation window
    cx, cy = cols[1], row_y[1]
    panel_label(f, cx, cy, "F", "Isolation-window distribution", "width of the first precursor window, in Th")
    px2, top2 = cx + 70, cy + 52
    ixs = Scale(0.5, 50, px2, px2 + pw_ - 14, log=True)
    iys = Scale(0, 0.9, top2 + ph, top2)
    axes(f, px2, top2, pw_ - 14, ph, ixs, iys, [1, 3, 10, 30], [0, 0.4, 0.8], xlabel="isolation width (Th, log scale)", ylabel="fraction of MS2", ylabel_dx=-46, tag="F")
    f.rect(ixs(15), top2, ixs(50) - ixs(15), ph, fill=AMBER_L, opacity=0.8)
    f.rect(ixs(0.5), top2, ixs(3) - ixs(0.5), ph, fill=BLUE_L, opacity=0.45)
    for w_, p in synth_run.WIDTHS:
        frac = float(np.mean(D["ms2_width"] == w_))
        f.rect(ixs(w_) - 5, iys(frac), 10, iys(0) - iys(frac), fill=BLUE)
    f.text(ixs(50), top2 - 8, f"≥15 Th wide: {(1 - M['IsolationWidth_MS2_FractionLe15'])*100:.0f}% of scans", size=13, weight="600", fill=AMBER, anchor="end", tag="F")
    f.text(ixs(0.5), top2 - 8, f"≤3 Th narrow: {np.mean(D['ms2_width'] <= 3)*100:.0f}%", size=13, weight="600", fill=BLUE, tag="F2")
    f.text(ixs(1.6) + 12, iys(0.8) + 4, "1.6 Th: 80%", size=13, weight="600", fill=BLUE, tag="F")
    f.line(ixs(15), top2, ixs(15), top2 + ph, AMBER, 1.6, dash="5 4")
    callout(cx, cy + ph + 112, [
        "FractionLe15 = n(width ≤ 15 Th) / n(width)",
        f"  = {M['IsolationWidth_MS2_FractionLe15']:.2f}   median = {M['IsolationWidth_MS2_Median']:.1f} Th",
        "narrow ≤3 Th and wide ≥15 Th are the",
        "  limits the acquisition heuristic uses"], TEAL)

    # ---- G: repeated precursor targets
    cx, cy = cols[2], row_y[1]
    panel_label(f, cx, cy, "G", "Repeated precursor-target behaviour", "target m/z of every MS2, by MS1-delimited cycle")
    px2, top2 = cx + 70, cy + 44
    half = (pw_ - 14 - 26) / 2
    nC = 24
    # DDA
    first_cycles = np.unique(D["cycle_of_ms2"])[400:400 + nC]
    gx = Scale(0, nC, px2, px2 + half)
    gy = Scale(350, 1000, top2 + ph, top2)
    axes(f, px2, top2, half, ph, gx, gy, [], [400, 700, 1000], ylabel="precursor m/z", ylabel_dx=-46, tag="G1")
    for ci, c in enumerate(first_cycles):
        for mz in D["ms2_mz"][D["cycle_of_ms2"] == c]:
            if 350 <= mz <= 1000:
                f.circle(gx(ci + 0.5), gy(mz), 2.3, fill=TEAL, opacity=0.9)
    f.text(px2 + half / 2, top2 + ph + 22, "this run (DDA)", size=13.5, anchor="middle", weight="600", fill=TEAL, tag="G1")
    f.text(px2 + half / 2, top2 + ph + 40, "cycle →", size=13, anchor="middle", fill=MUTED, tag="G1")
    # fixed targets
    px3 = px2 + half + 26
    gx2 = Scale(0, nC, px3, px3 + half)
    axes(f, px3, top2, half, ph, gx2, gy, [], [400, 700, 1000], grid=True, yticklabels=False, tag="G2")
    for ci in range(nC):
        for mz in fixed_grid:
            f.circle(gx2(ci + 0.5), gy(mz), 2.3, fill=PLUM, opacity=0.9)
    f.text(px3 + half / 2, top2 + ph + 22, "contrast: fixed target list", size=13.5, anchor="middle", weight="600", fill=PLUM, tag="G2")
    f.text(px3 + half / 2, top2 + ph + 40, "cycle →", size=13, anchor="middle", fill=MUTED, tag="G2")
    callout(cx, cy + ph + 112, [
        "targets rounded to 0.1 Th, per cycle",
        "TargetSetModalFraction = most common",
        "  target set / eligible cycles",
        f"  this run {M['AcquisitionCycle_TargetSetModalFraction']:.4f}  vs  fixed list {fixed_M['AcquisitionCycle_TargetSetModalFraction']:.2f}"], PLUM)

    # ---------------------------------------------------------------- footer
    f.line(60, 1438, W - 60, 1438, RULE, 1)
    f.text(60, 1466, "Retention times are handled in seconds and plotted in minutes. Panel G contrasts the run with a separate synthetic fixed-target scheme to show what the metric detects.",
           size=13.5, fill=MUTED, tag="foot")
    return f


if __name__ == "__main__":
    # Use the common build driver for portable output paths and layout linting.
    from build_all import main
    raise SystemExit(main(["--figures", "2"]))
