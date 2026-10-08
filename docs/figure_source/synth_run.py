"""Synthetic DDA LC-MS/MS run -> real prideQC metrics (metrics.py)."""

from __future__ import annotations

import numpy as np
from impl import Precursor, Spectrum, metrics

SEED = 11
CHARGE_P = {1: 0.03, 2: 0.52, 3: 0.30, 4: 0.08, 5: 0.02, 6: 0.01, 0: 0.04}
WIDTHS = [(1.6, 0.80), (1.2, 0.08), (2.0, 0.08), (0.7, 0.04)]


def tic_envelope(rt_s: np.ndarray) -> np.ndarray:
    m = rt_s / 60.0
    ramp = 1 / (1 + np.exp(-(m - 7.5) / 1.6))
    fall = 1 / (1 + np.exp((m - 83.0) / 1.2))
    humps = 1 + 0.30 * np.sin(m / 9.0) + 0.18 * np.sin(m / 3.1 + 1.0) + 0.10 * np.sin(m / 1.3)
    return 1.6e9 * ramp * fall * humps + 3e6


def build(duration_s=5400.0, n_ms2=(8, 12)):
    rng = np.random.default_rng(SEED)
    run = metrics.RunSummary()
    ms1_rt, ms1_tic, ms1_bp = [], [], []
    ms2_rt, ms2_charge, ms2_width, ms2_mz = [], [], [], []
    cycle_of_ms2 = []
    t, cyc = 2.0, 0
    ch_vals = list(CHARGE_P)
    ch_p = np.array(list(CHARGE_P.values()))
    ch_p /= ch_p.sum()
    w_vals = [w for w, _ in WIDTHS]
    w_p = np.array([p for _, p in WIDTHS])
    while t < duration_s:
        env = float(tic_envelope(np.array([t]))[0])
        tic = env * float(np.exp(rng.normal(0, 0.14)))
        if rng.random() < 0.0015:  # rare transient (spray instability)
            tic *= rng.choice([0.05, 0.08, 11.0])
        bp = max(1e4, tic * float(rng.uniform(0.035, 0.075)))
        k = 30  # a spectrum of 30 peaks: one base peak plus 29 smaller peaks that complete the TIC
        inten = np.full(k, (tic - bp) / (k - 1))
        inten[0] = bp
        run.consume_spectrum(Spectrum(1, t, np.linspace(400.0, 1600.0, k), inten))
        ms1_rt.append(t)
        ms1_tic.append(tic)
        ms1_bp.append(bp)
        n = int(np.clip(round(3 + 9.5 * (env / 1.6e9) ** 0.9 * rng.uniform(0.85, 1.15)), 2, 12))
        tt = t + 0.28
        for _ in range(n):
            z = int(rng.choice(ch_vals, p=ch_p))
            mz = float(np.clip(rng.lognormal(np.log(640), 0.28), 360, 1500))
            width = float(rng.choice(w_vals, p=w_p))
            prec = Precursor(
                mz=mz, charge=z, intensity=float(bp * rng.uniform(0.01, 0.6)), isolation_width=width
            )
            ms2_tic = tic * 0.012 * float(rng.uniform(0.3, 1.7))
            run.consume_spectrum(
                Spectrum(
                    2,
                    tt,
                    np.array([300.0, 900.0]),
                    np.array([ms2_tic * 0.2, ms2_tic * 0.8]),
                    precursors=(prec,),
                )
            )
            ms2_rt.append(tt)
            ms2_charge.append(z)
            ms2_width.append(width)
            ms2_mz.append(mz)
            cycle_of_ms2.append(cyc)
            tt += 0.17
        t = tt + float(rng.uniform(0.02, 0.10))
        cyc += 1
    vals = {m.key: m.value for m in metrics.QCMetricCalculator().calculate(run)}
    return dict(
        metrics=vals,
        ms1_rt=np.array(ms1_rt),
        ms1_tic=np.array(ms1_tic),
        ms1_bp=np.array(ms1_bp),
        ms2_rt=np.array(ms2_rt),
        ms2_charge=np.array(ms2_charge),
        ms2_width=np.array(ms2_width),
        ms2_mz=np.array(ms2_mz),
        cycle_of_ms2=np.array(cycle_of_ms2),
    )


def fixed_target_metrics(n_cycles=40, n_targets=12, seed=5):
    """Contrast case: the same targets in the same order in every MS1-delimited cycle (e.g. a scheduled/fixed-window scheme)."""
    rng = np.random.default_rng(seed)
    grid = np.round(np.linspace(420, 960, n_targets) + rng.normal(0, 3, n_targets), 1)
    run = metrics.RunSummary()
    t = 0.0
    for _c in range(n_cycles):
        run.consume_spectrum(Spectrum(1, t, np.array([400.0]), np.array([1e8])))
        t += 0.3
        for g in grid:
            run.consume_spectrum(
                Spectrum(
                    2,
                    t,
                    np.array([300.0]),
                    np.array([1e6]),
                    precursors=(
                        Precursor(mz=float(g), charge=2, intensity=1e6, isolation_width=25.0),
                    ),
                )
            )
            t += 0.12
    vals = {m.key: m.value for m in metrics.QCMetricCalculator().calculate(run)}
    return vals, grid


if __name__ == "__main__":
    d = build()
    m = d["metrics"]
    for k in (
        "NumberOfSpectra_MS1",
        "NumberOfSpectra_MS2",
        "ChromatographyDuration",
        "TIC_MS1_Area",
        "TIC_MS1_CV",
        "TIC_MS1_Median",
        "BasePeak_MS1_Mean",
        "BasePeak_All_Max",
        "ScanRate_MS1",
        "ScanRate_MS2",
        "FastestFrequency_MS1",
        "FastestFrequency_MS2",
        "AvgCycleTime_MS1",
        "ChargeRatio_3over2",
        "MS2_PrecursorCharge_Fractions",
        "IsolationWidth_MS2_FractionLe15",
        "IsolationWidth_MS2_Median",
        "IsolationPrecursorMz_MS2_UniqueFraction",
        "IsolationPrecursorMz_MS2_MaxRepeatFraction",
        "IsolationPrecursorMz_MS2_RepeatedFraction",
        "AcquisitionCycle_Count",
        "AcquisitionCycle_MS2Count_Mode",
        "AcquisitionCycle_TargetSetModalFraction",
        "AcquisitionCycle_TargetSetDistinctCount",
        "MS1_to_MS2_Ratio",
    ):
        print(k, m[k])
    fv, grid = fixed_target_metrics()
    print("--- fixed")
    for k in (
        "IsolationPrecursorMz_MS2_UniqueFraction",
        "IsolationPrecursorMz_MS2_MaxRepeatFraction",
        "IsolationPrecursorMz_MS2_RepeatedFraction",
        "AcquisitionCycle_TargetSetModalFraction",
        "AcquisitionCycle_TargetSetDistinctCount",
    ):
        print(k, fv[k])
