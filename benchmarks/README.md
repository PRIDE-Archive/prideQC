# Performance measurement

Use a fresh process for each measurement because peak RSS is a lifetime maximum.

```bash
uv run python benchmarks/benchmark.py --scans 10000 --peaks 1000
uv run python benchmarks/benchmark.py --scans 10000 --peaks 10000
uv run python benchmarks/benchmark.py --scans 100000 --peaks 1000
uv run python benchmarks/benchmark.py --file /data/representative.mzML
```

The synthetic stream reuses input arrays and exercises the actual reduction and
statistics pipeline. It isolates memory retained by prideqc and arithmetic
cost. It does not simulate XML parsing, decompression, disk I/O, vendor
conversion or distribution of peak shapes. It cannot establish a speedup over
rawQC/TechSDRF. Recorded synthetic results are in `docs/benchmark-*.json`. These are historical
0.1 measurements; the new SDRF validation and PRIDE acquisition paths have not
been benchmarked. Measure those stages separately from the spectral analysis.

For a defensible comparison:

1. Use the same files and converter settings for both tools, including DDA,
   DIA, profile, centroid, MS3 and representative high-peak-count files.
2. Pin both environments and record platform, CPU, RAM and storage. Compare
   the overlapping metric set and disclose changed definitions.
3. Report separately conversion, QC+annotation, and output time. Use at least
   three repeats and distinguish cold/warm filesystem cache.
4. Record wall time and peak RSS at 1, 2 and 4 workers. Avoid simultaneous jobs
   during serial comparison; parallel workers can saturate storage.
5. Validate numerical outputs before comparing runtime. Retain fixtures and
   environment locks with the measurements.

The memory contract is O(scans + largest spectrum) per worker, with small
metadata counters and optional evidence state. Quantiles are exact. The most
expensive stages are peak reduction, optional diagnostic/peak-type algorithms,
and sorting of unsorted scalar distributions. Files with many distinct scan
energies also grow the corresponding metadata counters. There is no peak cache
and no all-pairs precursor-mass comparison in the core.
