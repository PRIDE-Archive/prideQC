# TMT reporter-ion evidence (experimental)

This opt-in, identification-free feature collects TMT reporter-ion evidence **during prideQC's existing spectrum pass**. It does not require peptide identifications or a new reader. It does not identify the TMT kit, infer used biological channels, map channels to samples, or populate SDRF `comment[label]`.

## Running one file

With a compatible pyOpenMS installation (the historical native Thermo reader was verified using pyOpenMS `3.6.0.dev20260910`):

```bash
uv run --active prideqc analyze \
  /path/to/known-TMT.raw \
  --output-dir /path/to/new-output \
  --tmt-reporter-evidence
```

For a more conservative local smoke test, begin with the 156 MiB PXD030342 MS3 run, already used in the separate exploratory pilot. Use a **new** output directory; existing prideQC outputs are not overwritten. The native Thermo RAW loader materializes spectra and can require substantially more RAM than the file's compressed size.

## Evidence contract

- Emits a `tmt_reporter_evidence` annotation to summary JSON and (as a locally defined `QCPRIDE:` table) mzQC; never an SDRF label.
- Reads exact reporter m/z values from *runtime pyOpenMS* TMT6, TMT10, TMT11, TMTpro16 and TMTpro18 quantitation classes using `getChannelInformation()`; missing classes are recorded explicitly.
- Counts eligible **centroid** MS2 and MS3 spectra independently. Profile and unknown spectrum types are counted as excluded, **not** as reporter-negative.
- In each scan, a channel hit is the strongest finite positive peak within ±10 ppm of its expected mass, provided its intensity is at least 0.1% of that scan's base peak. These default thresholds are a screen, **not** calibration or statistical FDR.
- For each method and channel, reports hit count, fraction of eligible centroid spectra, uncorrected total matched intensity and mean absolute ppm error. A method-level screen is positive when at least three channels match within a spectrum.
- Methods share channel masses; multiple positive method signatures are *not independent evidence* of multiple kits. No isobaric impurity correction, interference correction, occupancy classification or spectral profile centroiding is performed.
- The opt-in switch leaves legacy `--diagnostics` behavior and previous summary outputs unchanged.

## Real-data pilot and limitations

Seven downloaded runs were previously screened separately (PXD000001, PXD030259, PXD005180, PXD030342, PXD023865, PXD016491 and PXD024275). Five yielded informative centroid reporter evidence. PXD000001 and PXD005180 have relevant profile spectra and remain inconclusive. PXD020831 was excluded after two identical downloaded `.raw.gz` files failed `gzip -t`; do not retry as part of this collector pilot.

The current patch has **offline synthetic tests only**, not a completed native RAW integration run. Compare one new output against the prior saved `tmt-evidence-v1` report before broad deployment, and use independent SDRF/protocol labels to adjudicate intended channels. In particular, distinguish TMTpro15-like and TMTpro16-like *signal patterns* without asserting exact plex or sample assignment.
