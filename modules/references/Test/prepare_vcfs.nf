process FILTER_AND_MERGE_VCF {
    tag "Filter and Merge VCF Files"
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'
	
    cpus params.get('filter_merge_vcf_cpus', 12)
    memory params.get('filter_merge_vcf_memory', '32 GB')
    time params.get('filter_merge_vcf_time', '6h')

    input: 
    path variants_snp
    path variants_snp_index 
    path variants_indels
    path variants_indels_index   
    path (denylist)  

    output:
    path "merged.filtered.recode.vcf.gz", emit: merged_vcf
    path "merged.filtered.recode.vcf.gz.tbi", emit: merged_vcf_tbi

    script:
    """
    # Number of threads for bcftools
    THREADS=${task.cpus}

    # Extract the denylist file path from the tuple
    DENYLIST_FILE=${denylist}

    echo "Filtering SNP variants..."
    bcftools view -T ^\$DENYLIST_FILE ${variants_snp} -Oz -o filtered_snps.vcf.gz --threads \$THREADS
    tabix -p vcf filtered_snps.vcf.gz  

    echo "Filtering INDEL variants..."
    bcftools view -T ^\$DENYLIST_FILE ${variants_indels} -Oz -o filtered_indels.vcf.gz --threads \$THREADS
    tabix -p vcf filtered_indels.vcf.gz  

    echo "Merging filtered SNPs and INDELs..."
    bcftools merge filtered_snps.vcf.gz filtered_indels.vcf.gz -Oz -o merged.filtered.recode.vcf.gz --threads \$THREADS
    tabix -p vcf merged.filtered.recode.vcf.gz  

    echo "Fixing chromosome naming..."
    # Generate a temporary chromosome rename file
    zcat merged.filtered.recode.vcf.gz | grep -v "^#" | cut -f1 | sort -u | grep "^chr" | awk '{print \$1"\t"substr(\$1,4)}' > rename_chr.txt


    # Apply renaming
    bcftools annotate --rename-chrs rename_chr.txt -Oz -o fixed_merged.filtered.recode.vcf.gz merged.filtered.recode.vcf.gz --threads \$THREADS

    tabix -p vcf fixed_merged.filtered.recode.vcf.gz

    echo "Replacing original merged VCF with corrected version..."
    mv fixed_merged.filtered.recode.vcf.gz merged.filtered.recode.vcf.gz
    mv fixed_merged.filtered.recode.vcf.gz.tbi merged.filtered.recode.vcf.gz.tbi

    # Cleanup temporary rename file
    rm rename_chr.txt
    """
}