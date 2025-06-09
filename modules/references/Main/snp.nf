// ========================== Download SNP VCF ========================== //
process CHECK_OR_DOWNLOAD_VARIANTS_SNP {
    tag "Download SNP VCF"
    container null
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    output:
    path "variants_snp.vcf.gz", emit: variants_snp

    when:
    !file("${params.main_data_dir}/reference/variants_snp.vcf.gz").exists()

    script:
    """
    echo "Downloading SNP VCF from: ${params.variants_snp_download_url}"
    wget -q -O variants_snp.vcf.gz ${params.variants_snp_download_url}
    """
}

// ========================== Download SNP VCF Index ========================== //
process CHECK_OR_DOWNLOAD_VARIANTS_SNP_INDEX {
    tag "Download SNP Index"
    container null
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    output:
    path "variants_snp.vcf.gz.tbi", emit: snp_index

    when:
    !file("${params.main_data_dir}/reference/variants_snp.vcf.gz.tbi").exists()

    script:
    """
    echo "Downloading index for SNP VCF..."
    wget -q -O variants_snp.vcf.gz.tbi ${params.variants_snp_index_download_url}
    """
}

// ========================== Index SNP VCF (if missing) ========================== //
process INDEX_SNP_VCF {
    tag "Index SNP VCF"
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    input:
    path(vcf_file)

    output:
    path("variants_snp.vcf.gz.tbi"), emit: snp_index

    when:
    !file("${params.main_data_dir}/reference/variants_snp.vcf.gz.tbi").exists()

    script:
    """
    echo "Indexing SNP VCF file..."
    tabix -p vcf ${vcf_file}
    """
}
