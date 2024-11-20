process fastqc_raw {
    tag "${sample_id}"
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    publishDir "${params.outdir}/fastqc/raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.zip"), emit: fastqc_zip
    path("*.html"), emit: fastqc_html

    script:
    """
    fastqc --outdir . ${reads.join(' ')}
    """
}