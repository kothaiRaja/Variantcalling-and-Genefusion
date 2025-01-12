process PREPARE_VCF_FILE {
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.actual_data_dir}", mode: 'copy'

    input: 
    path snpsFile
    path indelsFile
    path denylisted

    output:
    tuple path("merged.filtered.recode.vcf.gz"),
          path("merged.filtered.recode.vcf.gz.tbi")

    script:  
    """
    # Filter SNPs file
    bcftools view -T ^${denylisted} ${snpsFile} -Oz -o ${snpsFile.baseName}.filtered.recode.vcf.gz
    tabix -p vcf ${snpsFile.baseName}.filtered.recode.vcf.gz

    # Filter INDELs file
    bcftools view -T ^${denylisted} ${indelsFile} -Oz -o ${indelsFile.baseName}.filtered.recode.vcf.gz
    tabix -p vcf ${indelsFile.baseName}.filtered.recode.vcf.gz

    # Merge the filtered SNPs and INDELs into a single VCF file
    bcftools merge ${snpsFile.baseName}.filtered.recode.vcf.gz ${indelsFile.baseName}.filtered.recode.vcf.gz \
        -Oz -o merged.filtered.recode.vcf.gz

    # Create a tabix index for the merged VCF
    tabix -p vcf merged.filtered.recode.vcf.gz
	
	# Check and fix contig prefixes in the merged VCF file
    if zcat merged.filtered.recode.vcf.gz | grep -q "^##contig=<ID=chr"; then
        echo "Renaming contig prefixes..."
        zcat merged.filtered.recode.vcf.gz | sed 's/^##contig=<ID=chr/##contig=<ID=/' | bgzip > fixed_merged.filtered.recode.vcf.gz
        tabix -p vcf fixed_merged.filtered.recode.vcf.gz
        mv fixed_merged.filtered.recode.vcf.gz merged.filtered.recode.vcf.gz
        mv fixed_merged.filtered.recode.vcf.gz.tbi merged.filtered.recode.vcf.gz.tbi
    else
        echo "Contig prefixes are already correct."
    fi
    """
}
