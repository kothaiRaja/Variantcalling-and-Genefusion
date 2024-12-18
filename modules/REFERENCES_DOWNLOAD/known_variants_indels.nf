// Process to download the test variants VCF
process DOWNLOAD_TEST_VARIANTS_INDELS {
    tag "Download test variants VCF"
	container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf"

    script:
    """
    wget -q -O variants_indels.vcf ${params.test_data_known_indels }
    """
}

