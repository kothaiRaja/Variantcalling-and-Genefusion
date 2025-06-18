nextflow.enable.dsl = 2

include { GATK_BASERECALIBRATOR } from '../modules/gatk/baserecalibration/main.nf'
include { GATK_APPLYBQSR }       from '../modules/gatk/applybsqr/main.nf'

workflow BASE_RECALIBRATION {
    take:
    bam_input_ch           // (sample_id, strandedness, bam, bai)
    intervals_ch           // One interval file per entry (scattered intervals)
    reference_genome
    reference_genome_index
    reference_genome_dict
    merged_vcf
    merged_vcf_index

    main:
    ch_versions = Channel.empty()

    log.info " Starting Base Recalibration Workflow (per interval)..."

    // STEP 1: Run GATK BaseRecalibrator once per sample (not per interval)
    ch_recal_input = bam_input_ch.map { sample_id, strandedness, bam, bai -> 
        tuple(sample_id, strandedness, bam, bai)
    }

    GATK_BASERECALIBRATOR_RESULT = GATK_BASERECALIBRATOR(
        ch_recal_input,
        reference_genome,
        reference_genome_index,
        reference_genome_dict,
        merged_vcf,
        merged_vcf_index
    )

    ch_recal_tables = GATK_BASERECALIBRATOR.out.recal_table
    ch_versions     = ch_versions.mix(GATK_BASERECALIBRATOR.out.versions)

    ch_recal_tables.view { " Recalibration Table: $it" }

    // STEP 2: Combine BAM input with intervals
    ch_bam_interval = bam_input_ch
        .combine(intervals_ch)
        .map { sample_id, strandedness, bam, bai, interval ->
            tuple([sample_id], sample_id, strandedness, bam, bai, interval)
        }

    // STEP 3: Prepare recal table keyed by [sample_id]
    ch_recal_keyed = ch_recal_tables
        .map { sample_id, strandedness, recal_table ->
            tuple([sample_id], recal_table)
        }

    // STEP 4: Join BAMs+intervals with recal tables
    ch_apply_input = ch_bam_interval
        .join(ch_recal_keyed)
        .map { key, sample_id, strandedness, bam, bai, interval, recal_table ->
            tuple(sample_id, strandedness, bam, bai, recal_table, interval)
        }

    ch_apply_input.view { " ApplyBQSR Input: $it" }

    // STEP 5: Run ApplyBQSR for each BAM+interval+recal_table
    GATK_APPLYBQSR_RESULT = GATK_APPLYBQSR(
        ch_apply_input,
        reference_genome,
        reference_genome_index,
        reference_genome_dict
    )

    bams_base_recalibrated_ch = GATK_APPLYBQSR.out.recalibrated_bam
    ch_versions               = ch_versions.mix(GATK_APPLYBQSR.out.versions)

    bams_base_recalibrated_ch.view { " Final BQSR Bams: $it" }

    emit:
    recalibrated_bams = bams_base_recalibrated_ch
    versions          = ch_versions
}
