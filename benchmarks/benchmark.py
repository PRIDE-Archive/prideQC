"""Report wall time and peak RSS for a synthetic stream or an actual mzML file."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

from prideqc import __version__
from prideqc.models import RunMetadata, Spectrum, SpectrumSink
from prideqc.pipeline import Analyzer


class SyntheticReader:
    engine_version = "synthetic-stream-no-IO"

    def __init__(self, scans: int, peaks: int) -> None:
        self.scans = scans
        self.peaks = peaks

    def read(self, path: Path, sink: SpectrumSink) -> RunMetadata:
        mz = np.linspace(100, 2000, self.peaks)
        intensity = np.linspace(1, 10000, self.peaks)
        for index in range(self.scans):
            sink.consume_spectrum(Spectrum(1 if index % 3 == 0 else 2,
                                           index / 10, mz, intensity))
        return RunMetadata()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file", type=Path, help="Real mzML input; excludes writing outputs")
    parser.add_argument("--scans", type=int, default=10000)
    parser.add_argument("--peaks", type=int, default=1000)
    args = parser.parse_args()
    if args.scans < 1 or args.peaks < 1:
        parser.error("scans and peaks must be positive")
    started = time.perf_counter()
    analyzer = Analyzer() if args.file else Analyzer(SyntheticReader(args.scans, args.peaks))
    result = analyzer.analyze(args.file or Path("synthetic.mzML"))
    elapsed = time.perf_counter() - started
    metrics = {metric.key: metric.value for metric in result.metrics}
    try:
        import resource

        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rss_bytes = rss if sys.platform == "darwin" else rss * 1024
    except ImportError:
        rss_bytes = None
    print(json.dumps({
        "mode": "real-mzML" if args.file else "synthetic-no-file-IO",
        "input": str(args.file) if args.file else None,
        "prideqc": __version__, "python": platform.python_version(),
        "numpy": np.__version__, "engine": result.engine_version,
        "platform": platform.platform(), "scans": metrics["NumberOfSpectra"],
        "peak_data_points": metrics["NumberOfSpectralPeaks"],
        "wall_seconds": elapsed, "scans_per_second": metrics["NumberOfSpectra"] / elapsed,
        "peak_rss_bytes": rss_bytes,
        "note": (
            "Synthetic mode excludes mzML decoding, conversion and output I/O; not an upstream "
            "speedup comparison."
        ),
    }, indent=2))


if __name__ == "__main__":
    main()
