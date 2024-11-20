process SAMTOOLS_INDEX_BAM {
    tag "$replicateId"
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.19.1--h50ea8bc_0"
    publishDir "${params.outdir}/recalibrate/bam", mode: 'copy'

    input:
    tuple val(replicateId), path(bam)

    output:
    tuple val(replicateId), path("${bam.baseName}.bam"), path("${bam.baseName}.bam.bai")

    script:
    """
    if [ -f $bam ]; then
        samtools index -b $bam ${bam.baseName}.bam.bai
    else
        echo "Error: BAM file not found" >&2
        exit 1
    fi
    """
}