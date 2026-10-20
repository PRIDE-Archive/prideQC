process PRIDEQC_FINALIZE {
    tag "${accession}:${expected_units}-logical-runs"
    label 'process_low'
    container params.sif
    // QC-only result directories are staged as links by default. Python's
    // Path.rglob does not descend directory symlinks, yielding zero receipts.
    // Copy QC outputs (not RAW) to the finalizer's task work directory.
    stageInMode 'copy'
    // Publish BOTH declared audit outputs into the accession audit folder.
    // The old pattern 'final/*' missed nested _audit files in practice.
    publishDir { "${params.outdir}/${accession}/_audit" },
        mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.tokenize('/').last() }
    input:
    tuple val(accession), val(expected_units), path(run_dirs, stageAs: 'runs/*')
    path finalizer
    output:
    tuple val(accession), path('final/_audit/accession-status.json'),
        path('final/_audit/file-status.tsv'), emit: status
    script:
    """
    python3 ${finalizer} --accession ${accession} --expected-units ${expected_units} \
      --root runs --outdir final \
      ${params.allow_complete_study ? '--verified-complete-inventory' : ''}
    """
}
