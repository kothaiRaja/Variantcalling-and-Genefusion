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
    # Number of threads for bcftools
    THREADS=${task.cpus}

    # Extract denylist path
	DENYLIST_FILE="${denylist}"

	# Ensure denylist exists
	if [[ ! -s "\$DENYLIST_FILE" ]]; then
		echo "ERROR: Denylist file not found or empty at \$DENYLIST_FILE" >&2
		exit 1
	fi



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
	echo "Fixing chromosome naming..."

	# Use a fixed name for the rename file
	RENAME_FILE="rename_chr.tmp.txt"

	zcat merged.filtered.recode.vcf.gz | grep -v "^#" | cut -f1 | sort -u | grep "^chr" | awk '{print \$1"\t"substr(\$1,4)}' > "\$RENAME_FILE"

	bcftools annotate --rename-chrs "\$RENAME_FILE" \
    -Oz -o fixed_merged.filtered.recode.vcf.gz merged.filtered.recode.vcf.gz --threads \$THREADS

	tabix -p vcf fixed_merged.filtered.recode.vcf.gz

	echo "Replacing original merged VCF with corrected version..."
	mv fixed_merged.filtered.recode.vcf.gz merged.filtered.recode.vcf.gz
	mv fixed_merged.filtered.recode.vcf.gz.tbi merged.filtered.recode.vcf.gz.tbi

	# Cleanup
	rm "\$RENAME_FILE"

    """
}