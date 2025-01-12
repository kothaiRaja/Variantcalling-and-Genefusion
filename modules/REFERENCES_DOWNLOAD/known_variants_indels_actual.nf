// Process to download the variants VCF
process DOWNLOAD_VARIANTS_INDELS {
    tag "Download variants Indels VCF"
	container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf"

    script:
    """
    wget -q -O variants_indels.vcf ${params.actual_data_known_indels }
    """
}