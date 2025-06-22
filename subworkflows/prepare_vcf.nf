include { CHECK_OR_DOWNLOAD_DENYLIST } from '../modules/references/denylist.nf'
include { CHECK_OR_DOWNLOAD_VARIANTS_SNP } from '../modules/references/snp.nf'
include { CHECK_OR_DOWNLOAD_VARIANTS_SNP_INDEX } from '../modules/references/snp.nf'
include { INDEX_SNP_VCF } from '../modules/references/snp.nf'
include { CHECK_OR_DOWNLOAD_VARIANTS_INDELS } from '../modules/references/indels.nf'
include { CHECK_OR_DOWNLOAD_VARIANTS_INDELS_INDEX } from '../modules/references/indels.nf'
include { INDEX_INDEL_VCF } from '../modules/references/indels.nf'
include { FILTER_AND_MERGE_VCF } from '../modules/references/prepare_vcfs.nf'



workflow DOWNLOAD_AND_PREPARE_VARIANT_VCFS {

   
	
	main:

    // --------- DENYLIST ---------
    denylist_ch = Channel.empty()
    if (params.denylist_bed) {
        println " Using provided denylist BED: ${params.denylist_bed}"
        denylist_ch = Channel
            .fromPath(params.denylist_bed, checkIfExists: true)
            .map { [it] }
            .collect()
    } else {
        println " Downloading denylist BED..."
        denylist_ch = CHECK_OR_DOWNLOAD_DENYLIST()
    }

    // --------- SNP VCF ---------
    snp_vcf_ch = Channel.empty()
    if (params.variants_snp) {
        println " Using provided SNP VCF: ${params.variants_snp}"
        snp_vcf_ch = Channel
            .fromPath(params.variants_snp, checkIfExists: true)
            .map { [it] }
            .collect()
    } else {
        println " Downloading SNP VCF..."
        snp_vcf_ch = CHECK_OR_DOWNLOAD_VARIANTS_SNP()
    }

    
    // --------- SNP VCF Index ---------
    snp_index_ch = Channel.empty()
    if (params.variants_snp_index) {
        println " Using provided SNP VCF index: ${params.variants_snp_index}"
        snp_index_ch = Channel
            .fromPath(params.variants_snp_index, checkIfExists: true)
            .map { [it] }
            .collect()
    } else if (params.variants_snp_index_download_url) {
        println " Downloading SNP VCF index from URL..."
        snp_index_ch = CHECK_OR_DOWNLOAD_VARIANTS_SNP_INDEX()
    } else {
        println " Generating SNP VCF index from VCF..."
        snp_index_ch = INDEX_SNP_VCF(snp_vcf_ch)
    }

     

    // --------- INDEL VCF ---------
    indel_vcf_ch = Channel.empty()
    if (params.variants_indels) {
        println " Using provided INDEL VCF: ${params.variants_indels}"
        indel_vcf_ch = Channel
            .fromPath(params.variants_indels, checkIfExists: true)
            .map { [it] }
            .collect()
    } else {
        println " Downloading INDEL VCF..."
        indel_vcf_ch = CHECK_OR_DOWNLOAD_VARIANTS_INDELS()
    }

    // --------- INDEL VCF Index ---------
    indel_index_ch = Channel.empty()
    if (params.variants_indels_index) {
        println " Using provided INDEL VCF index: ${params.variants_indels_index}"
        indel_index_ch = Channel
            .fromPath(params.variants_indels_index, checkIfExists: true)
            .map { [it] }
            .collect()
    } else if (params.variants_indels_index_download_url) {
        println " Downloading INDEL VCF index from URL..."
        indel_index_ch = CHECK_OR_DOWNLOAD_VARIANTS_INDELS_INDEX()
    } else {
        println " Generating INDEL VCF index from VCF..."
        indel_index_ch = INDEX_INDEL_VCF(indel_vcf_ch)
    }
	
	// --------- FILTER + MERGE ---------

if (params.merged_filtered_vcf && file(params.merged_filtered_vcf).exists() &&
    params.merged_filtered_vcf_index && file(params.merged_filtered_vcf_index).exists()) {

    println " Using provided merged and filtered VCF + index:"
    println "  ${params.merged_filtered_vcf}"
    println "  ${params.merged_filtered_vcf_index}"

    ch_merged_vcf   = Channel.value(file(params.merged_filtered_vcf))
    ch_merged_index = Channel.value(file(params.merged_filtered_vcf_index))

} else {
    println " Running FILTER_AND_MERGE_VCF to create merged and filtered VCF..."

    merge_results = FILTER_AND_MERGE_VCF(
        snp_vcf_ch,
        snp_index_ch,
        indel_vcf_ch,
        indel_index_ch,
        denylist_ch
    )

    ch_merged_vcf   = merge_results.merged_vcf
    ch_merged_index = merge_results.merged_vcf_tbi
}

    
	
	emit:
    merged_vcf       = ch_merged_vcf
    merged_vcf_index = ch_merged_index
    

}
