"""Synthetic related MS2 pair + recurrent-family data, scored with prideQC's real mass_shift.py functions."""

from __future__ import annotations

import math

import numpy as np
from impl import mass_shift as ms

PROTON = 1.007276466621
RES = {
    "G": 57.02146,
    "A": 71.03711,
    "S": 87.03203,
    "P": 97.05276,
    "V": 99.06841,
    "T": 101.04768,
    "C": 103.00919,
    "L": 113.08406,
    "I": 113.08406,
    "N": 114.04293,
    "D": 115.02694,
    "Q": 128.05858,
    "K": 128.09496,
    "E": 129.04259,
    "M": 131.04049,
    "H": 137.05891,
    "F": 147.06841,
    "R": 156.10111,
    "Y": 163.06333,
    "W": 186.07931,
}
H2O = 18.010565
PHOSPHO = 79.966331


def fragments(seq: str, mod_pos: int | None, delta: float):
    """b/y 1+ and 2+ ions; residues >= mod_pos carry the modification."""
    masses = [
        RES[a] + (delta if (mod_pos is not None and i == mod_pos) else 0.0)
        for i, a in enumerate(seq)
    ]
    out = []
    n = len(seq)
    for i in range(1, n):
        b = sum(masses[:i]) + PROTON
        y = sum(masses[i:]) + H2O + PROTON
        out.append((b, f"b{i}", 1))
        out.append((y, f"y{n - i}", 1))
        out.append(((b + PROTON) / 2, f"b{i}", 2))
        out.append(((y + PROTON) / 2, f"y{n - i}", 2))
    return out, sum(masses) + H2O


def make_pair(seed=4):
    rng = np.random.default_rng(seed)
    seq = "AVLTSEGPSGFIYSEDQSR"  # arbitrary synthetic sequence (not a claim about any real peptide)
    k = 8  # modification on residue 9 (index 8), only in spectrum B
    fa, ma = fragments(seq, None, 0.0)
    fb, mb = fragments(seq, k, PHOSPHO)
    spectra = []
    for frags, M in ((fa, ma), (fb, mb)):
        mz, it = [], []
        for m, name, z in frags:
            if m < 120 or m > 1500:
                continue
            w = rng.lognormal(0, 0.8) * (1.0 if z == 1 else 0.35) * (1.0 if name[0] == "y" else 0.6)
            if rng.random() < 0.12:
                continue  # missing ion
            mz.append(m + rng.normal(0, 0.004))
            it.append(w * 1e4)
        nz = rng.uniform(150, 1450, 70)
        mz += list(nz)
        it += list(rng.lognormal(-1.2, 0.7, 70) * 1e4 * 0.5)
        mz = np.array(mz)
        it = np.array(it)
        o = np.argsort(mz)
        spectra.append((mz[o], it[o], M))
    return seq, k, spectra


def index(mz, it, precursor_mz, charge, ident):
    m_, i_ = ms._top_peak_arrays(mz, it, 60, precursor_mz=precursor_mz, precursor_exclusion_da=2.0)
    neutral = precursor_mz * charge - charge * PROTON
    return ms._IndexedSpectrum(
        ident,
        0.0,
        precursor_mz,
        charge,
        neutral,
        m_,
        i_,
        tuple(sorted({int(math.floor(v)) for v in m_})),
    )


def greedy(left, right, offset, tol=0.05):
    """Same greedy one-to-one matching as ms._matched_intensity_dot, returning index pairs."""
    i = j = 0
    pairs = []
    while i < left.size and j < right.size:
        d = float(right[j] - (left[i] + offset))
        if abs(d) <= tol:
            if (
                j + 1 < right.size
                and abs(float(right[j + 1] - (left[i] + offset))) < abs(d)
                and abs(float(right[j + 1] - (left[i] + offset))) <= tol
            ):
                j += 1
                continue
            if (
                i + 1 < left.size
                and abs(float(right[j] - (left[i + 1] + offset))) < abs(d)
                and abs(float(right[j] - (left[i + 1] + offset))) <= tol
            ):
                i += 1
                continue
            pairs.append((i, j))
            i += 1
            j += 1
        elif d < -tol:
            j += 1
        else:
            i += 1
    return pairs


def build():
    seq, k, spectra = make_pair()
    (mzA, itA, MA), (mzB, itB, MB) = spectra
    zA = zB = 2
    pA = (MA + 2 * PROTON) / 2
    pB = (MB + 2 * PROTON) / 2
    A = index(mzA, itA, pA, zA, 0)
    B = index(mzB, itB, pB, zB, 1)
    signed = B.neutral_mass - A.neutral_mass
    unchanged, total, sim = ms._relatedness(A, B, signed, fragment_match_da=0.05)
    res = {}
    for name, off in (("unchanged", 0.0), ("full", signed), ("half", signed / 2)):
        pairs = greedy(A.mz, B.mz, off)
        cnt, dot = ms._matched_intensity_dot(
            A.mz, A.intensity, B.mz, B.intensity, offset_da=off, tolerance_da=0.05
        )
        assert cnt == len(pairs), (name, cnt, len(pairs))
        res[name] = dict(offset=off, pairs=pairs, count=cnt, dot=dot)
    shared_bins = sorted(set(A.bins) & set(B.bins))
    return dict(
        seq=seq,
        A=A,
        B=B,
        rawA=(mzA, itA),
        rawB=(mzB, itB),
        signed=signed,
        unchanged=unchanged,
        total=total,
        sim=sim,
        res=res,
        shared=shared_bins,
        pA=pA,
        pB=pB,
    )


def family_data(sigma_ppm=2.5, seed=8):
    """Many accepted pairs -> clusters; real _cluster_observations and UniMod lookup."""
    rng = np.random.default_rng(seed)
    fams = [
        (1.00335, 40, 1500),
        (14.01565, 31, 1600),
        (15.99491, 58, 1700),
        (21.98194, 22, 1700),
        (28.0313, 27, 1600),
        (42.01057, 46, 1750),
        (57.02146, 36, 1750),
        (79.96633, 93, 1800),
    ]
    obs = []
    uid = 0
    for center, n, M in fams:
        sg = math.sqrt(2) * sigma_ppm * M / 1e6
        for _ in range(n):
            obs.append(
                ms._MassShiftObservation(
                    float(rng.normal(center, sg)), uid, uid + 1, M, float(rng.uniform(0.3, 0.9))
                )
            )
            uid += 2
    # background: unrelated accepted pairs scattered (never recur)
    for _ in range(140):
        obs.append(
            ms._MassShiftObservation(
                float(rng.uniform(1, 120)), uid, uid + 1, 1700.0, float(rng.uniform(0.25, 0.5))
            )
        )
        uid += 2
    M_rep = float(np.median([o.mean_neutral_mass for o in obs]))
    pair_sigma = math.sqrt(2) * sigma_ppm * M_rep / 1e6
    cluster_tol = max(0.01, 6.0 * pair_sigma)
    clusters = ms._cluster_observations(obs, cluster_tol, minimum_pairs=8, minimum_unique_spectra=6)
    return (
        obs,
        clusters,
        dict(M_rep=M_rep, pair_sigma=pair_sigma, cluster_tol=cluster_tol, sigma_ppm=sigma_ppm),
    )


if __name__ == "__main__":
    d = build()
    print(
        "signed delta",
        d["signed"],
        "unchanged",
        d["unchanged"],
        "total",
        d["total"],
        "sim",
        round(d["sim"], 3),
    )
    for k_, v in d["res"].items():
        print(k_, round(v["offset"], 4), v["count"], round(v["dot"], 3))
    print(
        "shared bins",
        len(d["shared"]),
        "A bins",
        len(d["A"].bins),
        "B bins",
        len(d["B"].bins),
        "kept A",
        d["A"].mz.size,
    )
    print("pA,pB", d["pA"], d["pB"])
    obs, cl, info = family_data()
    print(info)
    for c in cl:
        print(
            round(c.center_da, 4),
            round(c.sigma_da, 4),
            c.pair_support,
            c.unique_spectrum_support,
            round(c.median_similarity, 2),
        )
    recs = ms.load_openms_modifications()
    c = max(cl, key=lambda c: c.pair_support)
    m = ms._matching_modifications(c.center_da, 0.02, recs, maximum_candidates=24)
    print("cands", [(x["name"], round(x["residual_da"], 4), x["candidate_category"]) for x in m])
    near = [
        (r.delta_mass_da - 0, r.name, r.accession)
        for r in recs
        if abs(abs(r.delta_mass_da) - c.center_da) <= 0.12
    ]
    print(len(near))
    [print(" ", round(a, 4), b, cc) for a, b, cc in sorted(near)[:20]]
