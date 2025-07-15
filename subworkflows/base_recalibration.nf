nextflow.enable.dsl = 2

include { GATK_BASERECALIBRATOR } from '../modules/gatk/baserecalibration/main.nf'
include { GATK_APPLYBQSR }       from '../modules/gatk/applybsqr/main.nf'

workflow BASE_RECALIBRATION {

    take:
    bam_input     
	recalib_ch
    intervalsch                 
    reference_genome
    reference_genome_index
    reference_genome_dict
    known_snps_vcf
    known_snps_index
    known_indels_vcf
    known_indels_index

    main:

    ch_versions = Channel.empty()
//    log.info " Starting Base Recalibration Workflow (per interval)..."



    // Combine known sites (VCFs)
    // Merge VCFs
ch_known_sites_vcf = known_snps_vcf
    .concat(known_indels_vcf)
    .collect()

// Merge VCF indexes
ch_known_sites_index = known_snps_index
    .concat(known_indels_index)
    .collect()

   
    baserecalibrator_results = GATK_BASERECALIBRATOR(
        recalib_ch,
        reference_genome,
        reference_genome_index,
        reference_genome_dict,
        ch_known_sites_vcf,
        ch_known_sites_index
    )

    recal_tables_ch = baserecalibrator_results.recal_table
    ch_versions     = ch_versions.mix(baserecalibrator_results.versions)

    //  Join recalibration tables with BAMs by meta
    applybqsr_joined_ch = bam_input.join(recal_tables_ch)

    //  Combine with intervals again
   ch_bqsr_apply_input = applybqsr_joined_ch
    .combine(intervalsch)
    .map { it ->
    def bam_meta         = it[0]
    def bam              = it[1]
    def bai              = it[2]
    def table            = it[3]
    def interval_meta    = it[4]
    def interval_file    = it[5]
    def meta             = [id: bam_meta.id, strandedness: bam_meta.strandedness]
    tuple(meta, bam, bai, table, interval_file)
}



    // ️ Run ApplyBQSR
    GATK_APPLYBQSR(
        ch_bqsr_apply_input,
        reference_genome,
        reference_genome_index,
        reference_genome_dict
    )

    recalibrated_bams_ch = GATK_APPLYBQSR.out.recalibrated_bam
    ch_versions           = ch_versions.mix(GATK_APPLYBQSR.out.versions)

    emit:
    recalibrated_bams = recalibrated_bams_ch
    versions          = ch_versions
}
