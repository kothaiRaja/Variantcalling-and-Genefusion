process FILTER_AND_MERGE_VCF {
    tag "Filter and Merge VCF Files"
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    

    input: 
    path variants_snp
    path variants_snp_index 
    path variants_indels
    path variants_indels_index   
    path denylist

    output:
    path "merged.filtered.recode.vcf.gz", emit: merged_vcf
    path "merged.filtered.recode.vcf.gz.tbi", emit: merged_vcf_tbi

    script:
    """
    # Number of threads for bcftools
    THREADS=${task.cpus}

    # Filter SNPs with multi-threading
    bcftools view -T ^${denylist} ${variants_snp} -Oz -o filtered_snps.vcf.gz -@ \$THREADS
    tabix -p vcf filtered_snps.vcf.gz  

    # Filter INDELs with multi-threading
    bcftools view -T ^${denylist} ${variants_indels} -Oz -o filtered_indels.vcf.gz -@ \$THREADS
    tabix -p vcf filtered_indels.vcf.gz  

    # Merge the filtered SNPs and INDELs into a single VCF file with multi-threading
    bcftools merge filtered_snps.vcf.gz filtered_indels.vcf.gz -Oz -o merged.filtered.recode.vcf.gz -@ \$THREADS
    tabix -p vcf merged.filtered.recode.vcf.gz  

    # Check and fix contig prefixes before finalizing
    if zcat merged.filtered.recode.vcf.gz | grep -q "^##contig=<ID=chr"; then
        echo "Renaming contig prefixes..."
        zcat merged.filtered.recode.vcf.gz | sed 's/^##contig=<ID=chr/##contig=<ID=/' | bgzip -@ \$THREADS > fixed_merged.filtered.recode.vcf.gz
        tabix -p vcf fixed_merged.filtered.recode.vcf.gz
        mv fixed_merged.filtered.recode.vcf.gz merged.filtered.recode.vcf.gz
        mv fixed_merged.filtered.recode.vcf.gz.tbi merged.filtered.recode.vcf.gz.tbi
    else
        echo "Contig prefixes are already correct."
    fi
    """
}

