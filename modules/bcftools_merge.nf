//====================================================Merging all the vcf files============================================//

process BCFTOOLS_MERGE {
	container "https://depot.galaxyproject.org/singularity/bcftools%3A1.19--h8b25389_1"
	publishDir "${params.outdir}/BCFTOOLS_MERGE", mode: "copy"
	
    input:
    path vcfs

    output:
    path "merged_output.vcf.gz"
    path "merged_output.vcf.gz.tbi"
	path "variants_per_sample.tsv" 


    script:
    """
    # Ensure absolute paths for all VCFs
    vcfs_absolute=\$(realpath ${vcfs} | tr '\\n' ' ')
    
    echo "Merging VCFs: \$vcfs_absolute"

    # Run BCFtools merge
    bcftools merge \\
        \$vcfs_absolute \\
        -O z -o merged_output.vcf.gz

    # Index the merged VCF
    tabix -p vcf merged_output.vcf.gz

    # Generate tabular summary of variants per sample
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT\\t]\\n' merged_output.vcf.gz > variants_per_sample.tsv
    """
}