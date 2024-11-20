process SAMTOOLS_FLAGSTAT {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/flagstat", mode: "copy"

    input:
    tuple val(replicateId), path(bam)

    output:
    path "${replicateId}_flagstat.txt"

    script:
    """
    samtools flagstat ${bam} > ${replicateId}_flagstat.txt
    """
}