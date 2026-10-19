process PRIDEQC_FINALIZE {
    tag "${accession}:${expected_units}-logical-runs"
    label 'process_low'
    container params.sif
    publishDir { "${params.outdir}/${accession}" },
        mode: 'copy', pattern: 'final/*', overwrite: true,
        saveAs: { filename -> filename.replaceFirst(/^final\//, '') }
    input:
    tuple val(accession), val(expected_units), path(run_dirs, stageAs: 'runs/*')
    path finalizer
    output:
    tuple val(accession), path('final/_audit/accession-status.json'), emit: status
    script:
    """
    python3 ${finalizer} --accession ${accession} --expected-units ${expected_units} \
      --root runs --outdir final \
      ${params.allow_complete_study ? '--verified-complete-inventory' : ''}
    """
}
