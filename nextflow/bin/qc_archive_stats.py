#!/usr/bin/env python3
"""Archive-wide prideQC statistics for the PRIDE home page.

Reads an export runs.ndjson (one canonical run per line, sorted by project) and
writes one small JSON document (about 50 KB) with histograms and category
counts. Each distribution is given twice:

- runs:     every run counts once;
- projects: every project counts once, with the median (or the most common
            category) of its runs. One huge project (for example 237k MALDI
            runs) cannot then dominate the picture.

Default is a dry run: compute and write --output, no network call. --execute
also uploads the document to pride-mongo-ws POST /qc-stats (API key from the
environment variable named by --api-key-env; the key is never printed).
"""

from __future__ import annotations

import argparse
import bisect
import datetime as dt
import hashlib
import json
import math
import os
import statistics
import sys
import urllib.error
import urllib.request
from array import array
from collections import Counter
from collections.abc import Callable
from pathlib import Path
from typing import Any

UTC = dt.timezone.utc  # noqa: UP017 (Codon host runs Python 3.9)
SCHEMA_VERSION = 1
DEFAULT_URL = "https://pride-mongo-ws-hx-webadmin-11.wp-k8s.ebi.ac.uk/pride/mongo/ws"

PRECURSOR_TOL = "QCPRIDE:24117D188C618917D827"
FRAGMENT_TOL_PPM = "QCPRIDE:D2338E066BD0AD22664B"
FRAGMENT_TOL_DA = "QCPRIDE:3D7A0C71FC01C8616CDB"
FRAGMENT_ERR_DA = "QCPRIDE:61B6AF9C0CAB62CDBF63"
FRAGMENT_ERR_PPM = "QCPRIDE:7B04462C88F741A73E01"
ACQUISITION = "QCPRIDE:66216D97E9C2CB43C805"
MASS_SHIFTS = "QCPRIDE:DD1A0299E953F22C94FE"
RUN_DURATION = "MS:4000067"
MS1_SPECTRA = "MS:4000059"
MS2_SPECTRA = "MS:4000060"
MEDIAN_CHARGE = "MS:4000175"
MS1_CYCLE = "MS:4000192"

# Putative modification families that count for the "top modifications" chart:
# high support only, and not isotope/adduct/unknown artefacts.
PTM_CLASSES = {"putative-ptm", "putative-modification", "sample-prep-modification"}
PTM_TOP_N = 15


def _edges(start: float, stop: float, step: float) -> list[float]:
    n = int(round((stop - start) / step))
    return [round(start + i * step, 10) for i in range(n + 1)]


# Each histogram: title, unit, bin edges and whether the edges are log10 of the value.
# Bin i is [edges[i], edges[i+1]); values below edges[0] go to bin 0, values at
# or above edges[-2] go to the last bin, which is open ("≥ edges[-2]").
HISTOGRAMS: dict[str, dict[str, Any]] = {
    "precursor_tolerance_ppm": {
        "title": "Suggested precursor search tolerance", "unit": "ppm", "log10": False,
        "edges": _edges(0, 31, 1),
    },
    "fragment_tolerance_ppm": {
        "title": "Suggested fragment search tolerance (high-resolution MS2)", "unit": "ppm",
        "log10": False, "edges": _edges(0, 105, 5),
    },
    "fragment_tolerance_da": {
        "title": "Suggested fragment search tolerance (low-resolution MS2)", "unit": "Da",
        "log10": False, "edges": _edges(0, 1.05, 0.05),
    },
    "ms1_spectra": {
        "title": "MS1 spectra per run", "unit": "count", "log10": True, "edges": _edges(0, 6.25, 0.25),
    },
    "ms2_spectra": {
        "title": "MS2 spectra per run", "unit": "count", "log10": True, "edges": _edges(0, 6.25, 0.25),
    },
    "run_duration_min": {
        "title": "Chromatography (run) duration", "unit": "min", "log10": False, "edges": _edges(0, 250, 10),
    },
    "ms1_cycle_time_s": {
        "title": "MS1 median cycle time", "unit": "s", "log10": False, "edges": _edges(0, 10.5, 0.5),
    },
}

CATEGORIES = {
    "acquisition_method": "Acquisition method",
    "ms2_resolution": "MS2 resolution",
    "median_precursor_charge": "Median MS2 precursor charge",
}


def _metrics(run: dict) -> dict[str, Any]:
    return {m.get("accession"): m.get("value") for m in run.get("metrics") or []}


def _num(v: Any) -> float | None:
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return None
    return float(v) if v == v and v not in (float("inf"), float("-inf")) else None


def _tolerance(m: dict, acc: str) -> float | None:
    v = m.get(acc)
    return _num(v.get("suggested_tolerance")) if isinstance(v, dict) else None


def run_values(run: dict) -> tuple[dict[str, float], dict[str, str], set[tuple[str, str, float | None]]]:
    """Numeric values, categories and putative modification families of one run."""
    m = _metrics(run)
    values: dict[str, float] = {}

    def put(name: str, value: float | None, transform: Callable[[float], float] = lambda x: x) -> None:
        if value is not None and value >= 0:
            values[name] = transform(value)

    put("precursor_tolerance_ppm", _tolerance(m, PRECURSOR_TOL))
    put("fragment_tolerance_ppm", _tolerance(m, FRAGMENT_TOL_PPM))
    put("fragment_tolerance_da", _tolerance(m, FRAGMENT_TOL_DA))
    put("ms1_spectra", _num(m.get(MS1_SPECTRA)))
    put("ms2_spectra", _num(m.get(MS2_SPECTRA)))
    put("run_duration_min", _num(m.get(RUN_DURATION)), lambda s: s / 60.0)
    put("ms1_cycle_time_s", _num(m.get(MS1_CYCLE)))

    cats: dict[str, str] = {}
    acq = m.get(ACQUISITION)
    term = acq.get("cv_term") if isinstance(acq, dict) else None
    cats["acquisition_method"] = (
        str(term.get("name")) if isinstance(term, dict) and term.get("name") else "Not determined"
    )
    regime = ""
    for acc in (FRAGMENT_ERR_DA, FRAGMENT_ERR_PPM, FRAGMENT_TOL_PPM, FRAGMENT_TOL_DA):
        v = m.get(acc)
        if isinstance(v, dict) and v.get("resolution_regime"):
            regime = str(v["resolution_regime"])
            break
    cats["ms2_resolution"] = {"high-resolution": "High resolution", "low-resolution": "Low resolution"}.get(
        regime, "Not determined")
    charge = _num(m.get(MEDIAN_CHARGE))
    cats["median_precursor_charge"] = (
        "Not determined" if charge is None or charge <= 0
        else str(int(charge)) + "+" if charge < 5 else "5+ or more")

    ptms: set[tuple[str, str, float | None]] = set()
    shifts = m.get(MASS_SHIFTS)
    if isinstance(shifts, list):
        for item in shifts:
            if not isinstance(item, dict) or item.get("confidence") != "high-support":
                continue
            if item.get("classification") not in PTM_CLASSES:
                continue
            candidates = item.get("unimod_candidates") or []
            if candidates and isinstance(candidates[0], dict) and candidates[0].get("name"):
                delta = _num(item.get("delta_mass_da"))
                ptms.add((str(candidates[0]["name"]), str(candidates[0].get("unimod_accession") or ""),
                          round(delta, 4) if delta is not None else None))
    return values, cats, ptms


def _bin(edges: list[float], x: float) -> int:
    return min(max(bisect.bisect_right(edges, x) - 1, 0), len(edges) - 2)


def _quantiles(values: Any) -> dict[str, float] | None:
    data = sorted(values)
    if not data:
        return None

    def q(p: float) -> float:
        pos = (len(data) - 1) * p
        lo = int(pos)
        hi = min(lo + 1, len(data) - 1)
        return data[lo] + (data[hi] - data[lo]) * (pos - lo)

    return {"p05": q(0.05), "p25": q(0.25), "median": q(0.5), "p75": q(0.75), "p95": q(0.95)}


class Accumulator:
    """Streaming statistics; per-project state is kept only for the current project."""

    def __init__(self) -> None:
        self.run_values: dict[str, array] = {k: array("d") for k in HISTOGRAMS}
        self.project_values: dict[str, list[float]] = {k: [] for k in HISTOGRAMS}
        self.run_cats: dict[str, Counter] = {k: Counter() for k in CATEGORIES}
        self.project_cats: dict[str, Counter] = {k: Counter() for k in CATEGORIES}
        self.ptm_runs: Counter = Counter()
        self.ptm_projects: Counter = Counter()
        self.ptm_names: dict[str, str] = {}
        self.ptm_masses: dict[str, list[float]] = {}
        self.runs = 0
        self.projects = 0
        self.runs_with_mass_shift = 0
        self.projects_with_mass_shift = 0
        self.done: set[str] = set()
        self.current: str | None = None
        self._reset_project()

    def _reset_project(self) -> None:
        self.cur_values: dict[str, list[float]] = {k: [] for k in HISTOGRAMS}
        self.cur_cats: dict[str, Counter] = {k: Counter() for k in CATEGORIES}
        self.cur_ptms: set[str] = set()
        self.cur_mass_shift = False

    def _close_project(self) -> None:
        if self.current is None:
            return
        self.projects += 1
        for k, vals in self.cur_values.items():
            if vals:
                self.project_values[k].append(statistics.median(vals))
        for k, counter in self.cur_cats.items():
            if counter:
                # most common category of the project; ties broken by name for a stable result
                best = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
                self.project_cats[k][best] += 1
        for ptm in self.cur_ptms:
            self.ptm_projects[ptm] += 1
        if self.cur_mass_shift:
            self.projects_with_mass_shift += 1
        self.done.add(self.current)
        self._reset_project()

    def add(self, run: dict) -> None:
        accession = str(run.get("pxd_accession") or "")
        if accession != self.current:
            if accession in self.done:
                raise ValueError("runs.ndjson must be sorted by project; " + accession + " appears twice")
            self._close_project()
            self.current = accession
        values, cats, ptms = run_values(run)
        self.runs += 1
        for k, v in values.items():
            self.run_values[k].append(v)
            self.cur_values[k].append(v)
        for k, c in cats.items():
            self.run_cats[k][c] += 1
            self.cur_cats[k][c] += 1
        if any(m.get("accession") == MASS_SHIFTS for m in run.get("metrics") or []):
            self.runs_with_mass_shift += 1
            self.cur_mass_shift = True
        run_ptm_names: set[str] = set()
        for name, acc, delta in ptms:
            run_ptm_names.add(name)
            self.ptm_names[name] = acc
            if delta is not None:
                self.ptm_masses.setdefault(name, []).append(delta)
        for name in run_ptm_names:
            self.ptm_runs[name] += 1
        self.cur_ptms |= run_ptm_names

    def finish(self) -> dict:
        self._close_project()
        histograms = {}
        for k, spec in HISTOGRAMS.items():
            edges = spec["edges"]
            transform: Callable[[float], float] = math.log10 if spec["log10"] else (lambda x: x)
            out: dict[str, Any] = {key: spec[key] for key in ("title", "unit", "log10", "edges")}
            for who, vals in (("runs", self.run_values[k]), ("projects", self.project_values[k])):
                counts = [0] * (len(edges) - 1)
                zeros = 0
                for v in vals:
                    if spec["log10"] and v <= 0:
                        zeros += 1
                        continue
                    counts[_bin(edges, transform(v))] += 1
                out[who] = {"counts": counts, "n": len(vals), "quantiles": _quantiles(vals)}
                if spec["log10"]:
                    out[who]["zero"] = zeros
            histograms[k] = out
        categories = {}
        for k, title in CATEGORIES.items():
            keys = sorted(set(self.run_cats[k]) | set(self.project_cats[k]),
                          key=lambda c: (-self.project_cats[k][c], -self.run_cats[k][c], c))
            categories[k] = {"title": title, "items": [
                {"key": c, "runs": self.run_cats[k][c], "projects": self.project_cats[k][c]} for c in keys]}
        top = sorted(self.ptm_projects, key=lambda n: (-self.ptm_projects[n], -self.ptm_runs[n], n))[:PTM_TOP_N]
        return {
            "schemaVersion": SCHEMA_VERSION,
            "runs": self.runs,
            "projects": self.projects,
            "histograms": histograms,
            "categories": categories,
            "putativeModifications": {
                "title": "Most frequent putative modifications (high support)",
                "runsWithEvidence": self.runs_with_mass_shift,
                "projectsWithEvidence": self.projects_with_mass_shift,
                "items": [{"key": n, "unimod": self.ptm_names.get(n, ""), "runs": self.ptm_runs[n],
                           "projects": self.ptm_projects[n],
                           "deltaMassDa": round(statistics.median(self.ptm_masses[n]), 4) if self.ptm_masses.get(n) else None}
                          for n in top],
            },
        }


def compute(path: Path) -> dict:
    acc = Accumulator()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                acc.add(json.loads(line))
    return acc.finish()


def document(stats: dict, source: Path) -> dict:
    now = dt.datetime.now(UTC).replace(microsecond=0)
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return {
        "id": "prideqc-archive-stats-" + now.strftime("%Y%m%dT%H%M%SZ"),
        "generatedAt": now.isoformat().replace("+00:00", "Z"),
        "source": source.parent.name + "/" + source.name,
        "sourceSha256": digest.hexdigest(),
        "scientificApproval": False,
        **stats,
    }


def post(url: str, header: str, key: str, payload: dict) -> dict:
    req = urllib.request.Request(url, data=json.dumps(payload, allow_nan=False).encode(), method="POST",
                                 headers={"Content-Type": "application/json", header: key})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            return json.loads(resp.read() or b"{}")
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code} from {url}: {exc.read()[:300]!r}") from None


def main(argv: list | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs_ndjson", type=Path)
    ap.add_argument("--output", type=Path, required=True, help="where to write the stats JSON")
    ap.add_argument("--base-url", default=DEFAULT_URL)
    ap.add_argument("--api-key-header", default="mongo-ws-api-key")
    ap.add_argument("--api-key-env", default="PRIDE_MONGO_WS_API_KEY")
    ap.add_argument("--execute", action="store_true", help="also upload to pride-mongo-ws /qc-stats")
    args = ap.parse_args(argv)
    key = os.environ.get(args.api_key_env, "")
    if args.execute and not key:
        ap.error(f"--execute needs the API key in ${args.api_key_env}")
    doc = document(compute(args.runs_ndjson), args.runs_ndjson)
    text = json.dumps(doc, allow_nan=False, separators=(",", ":"))
    tmp = args.output.with_suffix(args.output.suffix + ".tmp")
    tmp.write_text(text + "\n", encoding="utf-8")
    tmp.replace(args.output)
    url = args.base_url.rstrip("/") + "/qc-stats"
    if args.execute:
        saved = post(url, args.api_key_header, key, doc)
        if saved.get("id") != doc["id"]:
            raise RuntimeError("server did not confirm the stats document")
    print(json.dumps({"event": "DONE", "mode": "execute" if args.execute else "dry_run", "id": doc["id"],
                      "runs": doc["runs"], "projects": doc["projects"], "bytes": len(text),
                      "output": str(args.output), "url": url}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
