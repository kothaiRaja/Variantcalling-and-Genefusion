process fastqc_trimmed {
    tag "${sample_id}"
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    publishDir "${params.outdir}/fastqc/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.zip"), path("*.html")

    script:
    """
    fastqc --outdir . ${reads.join(' ')}
    """
}
