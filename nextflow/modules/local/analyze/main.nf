process PRIDEQC_ANALYZE {
    tag "${accession}:${unit_dir.name}"
    label 'process_medium'
    container params.sif
    scratch true
    publishDir { "${params.outdir}/${accession}/_run_status" },
        mode: 'copy', pattern: 'result/result.json', overwrite: true,
        saveAs: { filename -> "${unit_dir.name}.json" }
    publishDir { "${params.outdir}/${accession}/_run_logs" },
        mode: 'copy', pattern: 'result/*.log', overwrite: true,
        saveAs: { filename -> "${unit_dir.name}/${filename.tokenize('/').last()}" }
    input:
    tuple val(accession), path(unit_dir)
    path analyzer
    output:
    tuple val(accession), path('result'), emit: unit_result
    script:
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    python3 ${analyzer} --unit ${unit_dir} --output result \
      ${params.estimate_mass_error ? '--estimate-mass-error' : ''} \
      ${params.estimate_mass_shifts ? '--estimate-mass-shifts' : ''}
    """
}
