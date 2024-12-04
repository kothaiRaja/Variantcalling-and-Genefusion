// Process to download the test variants VCF
process DOWNLOAD_TEST_VARIANTS {
    tag "Download test variants VCF"
	container null
    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "variants.vcf"

    script:
    """
    wget -q -O variants.vcf ${params.test_data_url}/subset_chr22.vcf.gz
    """
}