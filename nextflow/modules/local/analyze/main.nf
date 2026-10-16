process PRIDEQC_ANALYZE {
    tag "${accession}:${task_id}"
    label 'process_medium'
    container params.sif
    publishDir { "${params.outdir}/${accession}/repository-metadata/${task_id}_${slug}" },
        mode: 'copy', pattern: 'qc/*', overwrite: true,
        saveAs: { filename -> filename.replaceFirst(/^qc\//, '') }
    input:
    tuple val(accession), val(expected), val(task_id), val(slug), path(task_json)
    output:
    tuple val(accession), val(expected), path('qc/*.mzQC'), path('qc/*.obo'), emit: run
    script:
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    prideqc_run_file.py --task ${task_json} \
      --protocol ${params.protocol} \
      --output qc --download download \
      ${params.checksum_check ? '' : '--no-checksum-check'} \
      ${params.estimate_mass_error ? '--estimate-mass-error' : ''} \
      ${params.estimate_mass_shifts ? '--estimate-mass-shifts' : ''}
    """
}
