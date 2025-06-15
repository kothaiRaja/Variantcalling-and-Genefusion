process FILTER_AND_MERGE_VCF {
    tag "Filter and Merge VCF Files"
	label 'process_high'
	
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.main_data_dir}/reference", mode: 'copy'
	
   

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
    THREADS=${task.cpus}
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

    echo "Checking Chromosome Naming in merged VCF..."
    zgrep -v '^#' merged.filtered.recode.vcf.gz | cut -f1 | sort -u | head -1 | awk '{ print "  → Detected contig in merged VCF: " \$1 }'

    """
}