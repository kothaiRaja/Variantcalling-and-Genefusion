nextflow.enable.dsl = 2

// Include required processes
include { VCF2MAF } from '../modules/maftools/vcf2maf/main.nf'
include { MAF_VISUALIZATION } from '../modules/maftools/visualisation/main.nf'

workflow MAF_ANALYSIS {

    take:
    vcf_tuples
    fasta
    vep_cache
    rscript

    main:

    log.info " Starting maftools analysis workflow..."
	
	ch_versions = Channel.empty()

    // Step 1: Convert VCF to MAF
    ch_vcf2maf = VCF2MAF(
        vcf_tuples,
        fasta,
        vep_cache
    )

    ch_maf = VCF2MAF.out.maf
    ch_versions = ch_versions.mix(VCF2MAF.out.versions)

    // Step 2: Visualize MAF
    ch_visualisation = MAF_VISUALIZATION(
        ch_maf,
        rscript
    )

    ch_pdf = MAF_VISUALIZATION.out.maf_plots
    ch_versions = ch_versions.mix(MAF_VISUALIZATION.out.versions)

    emit:
    maf         = ch_maf
    maf_plots   = ch_pdf
    versions    = ch_versions
}
