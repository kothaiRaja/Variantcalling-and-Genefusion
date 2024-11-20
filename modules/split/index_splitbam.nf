// Process to run samtools indexing
process SAMTOOLS_INDEX_SPLIT_BAM {
    container "https://depot.galaxyproject.org/singularity/samtools%3A0.1.18--h50ea8bc_13"
    publishDir "${params.outdir}/SNG/splitbam", mode: 'copy'

    input:
    tuple val(replicateId), path(bam)

    output:
    tuple val(replicateId), path("split.bam"), path("split.bam.bai")

    script:
    """
    samtools index $bam
    """
}
