process BCFTOOLS_STATS {
    tag "bcftools_stats"

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.21--h8b25389_0" 
    publishDir "${params.outdir}/bcftools_stats/beforefilteration", mode: "copy"

    input:
    path(vcf_file)    // Input VCF file
    path(vcf_index)   // Input VCF index (.tbi) file
	path(tsv_file)

    output:
    path("stats.txt")     // Output stats file

    script:
    """
    # Generate stats
    bcftools stats ${vcf_file} > stats.txt

    """
}