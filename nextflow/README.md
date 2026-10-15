# PRIDEQC Nextflow DSL2 (pilot)

Entry point: `main.nf`; workflow: `workflows/prideqc.nf`; local modules:
`PRIDEQC_MANIFEST`, `PRIDEQC_ANALYZE`, `PRIDEQC_MERGE`.

This is a metadata-only, frozen-manifest **pilot**, not an archive-wide production
service. See `../docs/nextflow-codon.md` for safety, Codon prerequisites, protocol
selection, MongoDB staging and complete commands. Read that file before running.

No GitHub push is needed to develop this on a local feature branch. An updated,
checksum-pinned SIF must be built and staged independently before using on Codon.
