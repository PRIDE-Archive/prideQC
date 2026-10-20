# PRIDEQC Nextflow DSL2 — archive-aware pilot

Entry point: `main.nf`; version 2 adapters under `bin/` and modules under
`modules/local/`. Requires a pinned SIF with prideqc + pridepy, Nextflow DSL2,
and a frozen **file-task** and **accession-report** manifest.

## October 19 safe first-flight lane (default)

`PRIDEQC_SOURCE_QC` runs download, extraction and analysis **inside the same
Slurm task scratch directory**. It accepts exactly **one logical RAW/mzML run
per source**, records a failure sentinel for a multi-run archive, and stages
out only per-run `result.json`, logs, `.mzQC`, `.obo`, and preparation status.
It does **not** route RAW directories via Nextflow's persistent work dir.

With the default Boolean `allow_complete_study=false`, the workflow writes an accession
`_audit/accession-status.json` with `status=partial` for an otherwise
successful pilot, and **never** publishes a study `<PXD>.mzQC`.
The workflow now **rejects** `allow_complete_study=true` outright until
an independently verified full-inventory attestation gate exists. Do not pass
`--allow_complete_study false` as a command-line string; use a typed JSON
params file containing `{ "allow_complete_study": false }`. The old separately
staged PREPARE/ANALYZE modules remain in source for historical reference but
are **not invoked** by the pilot workflow. They must not be re-enabled
without a RAW-free handoff design. Multi-run bundles and production
per-logical-run scatter require a separate architecture milestone.

## Pipeline

1. `PRIDEQC_MANIFEST`: fail-closed validation of frozen source-archive tasks.
2. Historical, inactive `PRIDEQC_PREPARE`: one pridepy download for each listed source file, then
   safely expand a ZIP or TAR.GZ to one *logical MS run* per unit directory.
   Supported inputs: `.raw`, `.raw.zip`, `.raw.gz`, `.mzML`, `.mzML.gz`,
   `.d.zip`, `.d`, `.tar.gz`, `.tgz` and generic bundles ending in `.zip`.
   ZIPs may contain several runs; `.d` folders are never split apart.
   Nested ZIP/TAR.GZ supported to bounded depth. Unsupported/corrupt archives
   produce failure sentinels and cannot lead to a complete study.
3. Historical, inactive `PRIDEQC_ANALYZE`: one isolated QC job per logical run. Every reached run
   emits `result.json` with success or failure, and a per-run mzQC only on success.
4. Active `PRIDEQC_SOURCE_QC` performs preparation and QC within one scratch task;
   it refuses source bundles with more than one logical run and stages out QC only.
5. `PRIDEQC_FINALIZE`: dynamically wait for the discovered number of unit
   receipts across the expected source-archive tasks. Write `_audit` JSON/TSV;
   **only** complete success permits a lossless study `_study/*.mzQC` merge.

`expected_file_tasks` counts original PRIDE archive entries, whereas
`expected_logical_runs` is calculated after validated extraction. These must
never be equated when a source ZIP contains multiple RAW files.

## Current limitations and production gates

* Requires `nextflow lint .`, local/remote profile resolution, Nextflow graph
  preview, a small real-archive end-to-end test, and Codon Slurm execution.
* Python offline tests do **not** prove FIRE availability, vendor `.d` reading,
  SIF completeness, resource limits, or Nextflow output/publish semantics.
* Supported v2 source types are metadata-only; SDRF-aligned per-member mapping,
  WIFF/Waters readers and arbitrary archive formats require separate work.
* Archive extraction is bounded to 500 GiB uncompressed per source by default;
  large submissions need a per-cluster capacity policy, not a silent override.
* A scheduler/node crash may prevent a per-unit receipt; Nextflow will report
  task failure and the study must be treated as incomplete, not as a clean run.
* Nextflow workflow completion must **not** be interpreted as accession QC
  success. Check `_audit/accession-status.json` (`status == complete`).
* MongoDB write is intentionally disabled. Export only validated complete mzQC
  through the separately reviewed `export_mongo_ndjson.py` artifact.
* PRIDE team-specific Codon instructions are maintained privately; do not
  publish internal filesystem paths, passwords or FIRE credentials.

## Mixed vendor RAW and converted mzML (reviewed representation manifest)

**One physical MS acquisition is one logical run**, even if PRIDE contains both
an original vendor file and a converted `.mzML`/`.mzML.gz` version. The formats
are alternate *representations*, not independent biological/technical replicates.
However, matching file stems alone **do not prove** equivalence: mzML conversion
may apply centroiding, metadata filters or other transformations.

Build a curated TSV with these mandatory columns:
`pxd_accession`, `file_name`, `file_category` (`RAW` or `PEAK`),
`size_bytes`, `logical_run_id`, `equivalence_evidence`,
`native_reader_verified` (`0` or `1`, for native `.raw` only).

Every row must have a reviewed acquisition ID. In a mixed vendor/mzML project,
all rows need explicit provenance strings starting with `MSRUN:`, `SDRF:`,
`CONVERSION:` or `CURATED:`. Different representations assigned to the same
logical run must share the same evidence string; the operator must validate
that the evidence is genuine. This convention is **not** automated verification
of an arbitrary claimed MSRun ID or annotation. A raw/converted pair can be
compared but must not be double-counted. Unsupported `.d`-folder archives,
unknown file mappings, unsafe basenames and size failures block selection.

```bash
python3 nextflow/bin/build_representation_manifest.py \
  --inventory /path/to/reviewed-all-format-files.tsv \
  --outdir /path/to/new-representation-manifest \
  --prefer native
```

Outputs: `file-manifest.tsv`, `accession-reports.tsv`,
`representation-decisions.tsv` (retains unused alternative file identities),
and `representation-provenance.json` (includes original TSV SHA256). A native
reader is considered eligible only after a separately verified reader test;
otherwise a *verified-equivalent* mzML is selected if available. The regular
`plan_manifest.py` rejects direct mixed vendor/mzML task lists lacking a
curated inventory hash and per-acquisition identity. Merely omitting the
converted files from a manifest **does not establish complete-study coverage**.
The selection tool always marks whole-study verification and merge authorization
as false.

## Pilot scientific and operations qualification

* Active Codon `PRIDEQC_SOURCE_QC` is limited to two concurrent 2-CPU,
  32-GB, 8-hour tasks by the checked-in profile. Retain scratch and bounded
  concurrency limits for the next multi-file pilot.
* Failed `ThermoRawFile` opens are classified only as
  `thermo_raw_open_error_unclassified`. This is not a corruption diagnosis.
  Inspect project MSRun/instrument metadata, upstream checksum, the complete
  file inventory (both RAW and PEAK categories), and a verified-equivalent
  open-format alternative before changing the reader or repeating a download.
* `nextflow/bin/inspect_pride_raw_metadata.py` inspects RAW-category file
  metadata without downloading data. It cannot rule out mzML in PEAK or other
  categories.
* Source-QC receipts now record requested mass-error and mass-shift estimator
  flags and optional curated representation identity. `audit_science_metrics.py`
  inventories mzQC metric groups but explicitly **does not validate** their
  numeric correctness, TMT chemistry, PTM localization, probability calibration,
  CV provenance or actual FIRE transfer.
* A full end-to-end flight of the patched source and typed JSON parameters,
  with explicit estimator enablement and a separate scientific review, remains
  necessary before any PRIDE-wide rollout.
