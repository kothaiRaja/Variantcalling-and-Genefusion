// Process to download the test variants VCF
process DOWNLOAD_VARIANTS_INDELS {
    tag "Download test variants VCF"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf"

    when:
    !file("${params.test_data_dir}/reference/variants_indels.vcf").exists()

    script:
    """
    wget -q -O variants_indels.vcf ${params.test_data_known_indels}
    """
}
