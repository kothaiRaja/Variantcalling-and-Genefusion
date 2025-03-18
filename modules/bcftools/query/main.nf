process BCFTOOLS_QUERY {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.3--h577a1d6_9"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_query.*"

    input:
    tuple val(sample_id), path(filtered_vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}_variant_summary.txt")

    script:
    """
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n' ${filtered_vcf} > ${sample_id}_variant_summary.txt
    """
}
