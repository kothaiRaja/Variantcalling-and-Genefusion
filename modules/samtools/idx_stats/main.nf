process SAMTOOLS_IDXSTATS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_idxstats.*"

    input:
    tuple val(sample_id), path(bam), path(bai), val(strandedness)

    output:
    tuple val(sample_id), path("${sample_id}_idxstats.txt"), val(strandedness)

    script:
    """
    samtools idxstats ${bam} > ${sample_id}_idxstats.txt
    """
}