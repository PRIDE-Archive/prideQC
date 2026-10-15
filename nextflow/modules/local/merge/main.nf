process PRIDEQC_MERGE {
    tag "${accession}:expected-${expected}"
    label 'process_low'
    container params.sif
    publishDir { "${params.outdir}/${accession}/_study" },
        mode: 'copy', pattern: 'study/*', overwrite: true,
        saveAs: { filename -> filename.replaceFirst(/^study\//, '') }
    input:
    tuple val(accession), val(expected), path(mzqcs, stageAs: 'input/*'), path(obos, stageAs: 'input/*')
    output:
    tuple val(accession), path('study/*.mzQC'), path('study/*.obo'),
        path('study/*.manifest.json'), emit: study
    shell:
    '''
    mkdir -p study
    prideqc mzqc merge \
      --input-dir input \
      --output study/!{accession}.mzQC \
      --accession !{accession} \
      --expected-runs !{expected} \
      --cv-uri 'file://!{params.outdir}/!{accession}/_study/!{accession}.obo'
    '''
}
