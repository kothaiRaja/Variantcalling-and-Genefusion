// Process to download the actual variants VCF
process DOWNLOAD_VARIANTS_SNP {
    tag "Download actual variants VCF"
	container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "variants_snp.vcf"
	
	when:
    !file("${params.actual_data_dir}/reference/variants_snp.vcf").exists()

    script:
    """
    wget -q -O variants_snp.vcf ${params.actual_data_dbsnp}
    """
}