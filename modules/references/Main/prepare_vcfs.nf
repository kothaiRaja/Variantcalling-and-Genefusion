process FILTER_AND_MERGE_VCF {
    tag "Filter and Merge VCF Files"
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'
	
    cpus params.get('filter_merge_vcf_cpus', 12)
    memory params.get('filter_merge_vcf_memory', '32 GB')
    time params.get('filter_merge_vcf_time', '6h')

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

    # Filter SNPs
    bcftools view -T ^${denylist} ${variants_snp} -Oz -o filtered_snps.vcf.gz --threads \$THREADS
    tabix -p vcf filtered_snps.vcf.gz  

    # Filter INDELs
    bcftools view -T ^${denylist} ${variants_indels} -Oz -o filtered_indels.vcf.gz --threads \$THREADS
    tabix -p vcf filtered_indels.vcf.gz  

    # Merge SNPs and INDELs
    bcftools merge filtered_snps.vcf.gz filtered_indels.vcf.gz -Oz -o merged.filtered.recode.vcf.gz --threads \$THREADS
    tabix -p vcf merged.filtered.recode.vcf.gz  

    # Fix chromosome naming in both header and variant lines
    zcat merged.filtered.recode.vcf.gz | awk '
        BEGIN {OFS="\t"} 
        /^##contig=<ID=chr/ {sub(/^##contig=<ID=chr/, "##contig=<ID=")} 
        !/^#/ && \$1 ~ /^chr/ {sub(/^chr/, "", \$1)} 
        {print}
    ' | bgzip --threads \$THREADS > fixed_merged.filtered.recode.vcf.gz

    tabix -p vcf fixed_merged.filtered.recode.vcf.gz
    mv fixed_merged.filtered.recode.vcf.gz merged.filtered.recode.vcf.gz
    mv fixed_merged.filtered.recode.vcf.gz.tbi merged.filtered.recode.vcf.gz.tbi
    """
} 
