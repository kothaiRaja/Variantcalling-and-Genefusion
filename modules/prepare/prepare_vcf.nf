//Filtered and recorded variants

process PREPARE_VCF_FILE {
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.outdir}/genome/vcf", mode: "copy"
	 
    input: 
        path variantsFile
        path denylisted

    output:
        tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"),
              path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")
  
    script:  
    """
    # Exclude regions from the blacklist using bcftools
    bcftools view -T ^${denylisted} ${variantsFile} -Oz -o ${variantsFile.baseName}.filtered.recode.vcf.gz

    # Index the filtered VCF
    tabix -p vcf ${variantsFile.baseName}.filtered.recode.vcf.gz
    """
}

