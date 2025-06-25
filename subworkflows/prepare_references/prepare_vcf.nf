include { CHECK_OR_DOWNLOAD_ALL_VARIANTS } from '../../modules/references/variants.nf'
include { SPLIT_AND_INDEX_VCF }            from '../../modules/references/prepare_vcf.nf'
include { INDEX_VCF }                      from '../../modules/references/index_vcf.nf'

workflow SPLIT_AND_INDEX_ALL_VARIANTS {

    take:
    all_variants_vcf_optional

    main:
    snps_vcf_ch = Channel.empty()
    snps_index_ch = Channel.empty()
    indels_vcf_ch = Channel.empty()
    indels_index_ch = Channel.empty()

    if (params.all_variants_vcf) {
        println "Using provided GCF VCF file: ${params.all_variants_vcf}"
        all_variants_vcf_ch = Channel.fromPath(params.all_variants_vcf, checkIfExists: true)
        split_results = SPLIT_AND_INDEX_VCF(all_variants_vcf_ch)

        snps_vcf_ch     = split_results.out.snps_vcf
        snps_index_ch   = split_results.out.snps_index
        indels_vcf_ch   = split_results.out.indels_vcf
        indels_index_ch = split_results.out.indels_index

    } else if (params.known_snps_vcf && params.known_indels_vcf) {
        println "Using separately provided known SNPs and INDELs VCF files."

        snps_vcf_ch   = Channel.fromPath(params.known_snps_vcf, checkIfExists: true)
        indels_vcf_ch = Channel.fromPath(params.known_indels_vcf, checkIfExists: true)

        snps_index_path   = "${params.known_snps_vcf}.tbi"
        indels_index_path = "${params.known_indels_vcf}.tbi"

        if (!new File(snps_index_path).exists()) {
            snps_index_result = INDEX_VCF(snps_vcf_ch)
            snps_index_ch     = snps_index_result.index
        } else {
            snps_index_ch = Channel.fromPath(snps_index_path, checkIfExists: true)
        }

        if (!new File(indels_index_path).exists()) {
            indels_index_result = INDEX_VCF(indels_vcf_ch)
            indels_index_ch     = indels_index_result.index
        } else {
            indels_index_ch = Channel.fromPath(indels_index_path, checkIfExists: true)
        }
    }
	else {
    println "No known variant VCFs provided. Attempting to download 00-All VCF..."
    all_variants_vcf_ch = CHECK_OR_DOWNLOAD_ALL_VARIANTS().out.all_variants_vcf
    split_results = SPLIT_AND_INDEX_VCF(all_variants_vcf_ch)

    snps_vcf     = split_results.out.snps_vcf
    snps_index   = split_results.out.snps_index
    indels_vcf   = split_results.out.indels_vcf
    indels_index = split_results.out.indels_index
}


    emit:
    snps_vcf     = snps_vcf_ch
    snps_index   = snps_index_ch
    indels_vcf   = indels_vcf_ch
    indels_index = indels_index_ch
}
