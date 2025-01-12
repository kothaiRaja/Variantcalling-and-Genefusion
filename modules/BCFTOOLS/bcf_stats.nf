process BCFTOOLS_STATS {
    tag "bcftools_stats"

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.21--h8b25389_0" 
    publishDir "${params.outdir}/bcftools_stats", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)   

    output:
    tuple val(sample_id), path("haplotypecaller_stats_${sample_id}.txt")   

    script:
    """
    # Generate stats
    bcftools stats output_${sample_id}.vcf.gz > haplotypecaller_stats_${sample_id}.txt

    """
}