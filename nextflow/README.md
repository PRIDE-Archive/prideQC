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

The default `--allow_complete_study false` always writes an accession
`_audit/accession-status.json` with `status=partial` for an otherwise
successful pilot, and **never** publishes a study `<PXD>.mzQC`.
Only an externally verified full-inventory run may opt into
`--allow_complete_study true`; this option is an **operator assertion**, not
independent proof of PRIDE-wide file completeness. The old separately
staged PREPARE/ANALYZE modules remain in source for historical reference but
are **not invoked** by the pilot workflow. They must not be re-enabled
without a RAW-free handoff design. Multi-run bundles and production
per-logical-run scatter require a separate architecture milestone.

## Pipeline

1. `PRIDEQC_MANIFEST`: fail-closed validation of frozen source-archive tasks.
2. `PRIDEQC_PREPARE`: one pridepy download for each listed source file, then
   safely expand a ZIP or TAR.GZ to one *logical MS run* per unit directory.
   Supported inputs: `.raw`, `.raw.zip`, `.raw.gz`, `.mzML`, `.mzML.gz`,
   `.d.zip`, `.d`, `.tar.gz`, `.tgz` and generic bundles ending in `.zip`.
   ZIPs may contain several runs; `.d` folders are never split apart.
   Nested ZIP/TAR.GZ supported to bounded depth. Unsupported/corrupt archives
   produce failure sentinels and cannot lead to a complete study.
3. `PRIDEQC_ANALYZE`: one isolated QC job per logical run. Every reached run
   emits `result.json` with success or failure, and a per-run mzQC only on success.
4. `PRIDEQC_FINALIZE`: dynamically wait for the discovered number of unit
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
