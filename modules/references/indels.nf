// ========================== Download Indels Variants VCF ========================== //
process CHECK_OR_DOWNLOAD_VARIANTS_INDELS {
    tag "Download INDEL VCF"
    container null
    publishDir "${params.ref_base}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf.gz", emit: variants_indels


    script:
    """
    echo "Downloading INDEL VCF from: ${params.variants_indels_download_url}"
    wget -q -O variants_indels.vcf.gz ${params.variants_indels_download_url}
    """
}

// ========================== Download INDEL Index ========================== //
process CHECK_OR_DOWNLOAD_VARIANTS_INDELS_INDEX {
    tag "Download INDEL Index"
    container null
    publishDir "${params.ref_base}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf.gz.tbi", emit: vcf_indel_index

    

    script:
    """
    echo "Downloading INDEL index from: ${params.variants_indels_index_download_url}"
    wget -q -O variants_indels.vcf.gz.tbi ${params.variants_indels_index_download_url}
    """
}

// ========================== Index INDEL VCF (if needed) ========================== //
process INDEX_INDEL_VCF {
    tag "Index INDEL VCF"
    label 'process_low'
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.ref_base}/reference", mode: 'copy'

    input:
    path vcf_file

    output:
    path "*.tbi", emit: indels_index


    script:
    """
    echo "Indexing INDEL VCF: ${vcf_file}"
    tabix -p vcf ${vcf_file}
    """
}
