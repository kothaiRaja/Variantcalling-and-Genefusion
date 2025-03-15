process BCFTOOLS_STATS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.3--h577a1d6_9"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_bcftools_stats.*"

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}_bcftools_stats.txt")

    script:
    """
    bcftools stats ${vcf} > ${sample_id}_bcftools_stats.txt
    """
}