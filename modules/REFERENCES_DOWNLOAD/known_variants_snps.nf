// Process to download the test variants VCF
process DOWNLOAD_TEST_VARIANTS_SNP {
    tag "Download test variants VCF"
	container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_snp.vcf"

    script:
    """
    wget -q -O variants_snp.vcf ${params.test_data_dbsnp}
    """
}