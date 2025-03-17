nextflow.enable.dsl = 2

// Import necessary processes
include { SAMTOOLS_FILTER_ORPHANS } from '../modules/samtools/filter_orphans/main.nf'
include { GATK_MARK_DUPLICATES } from '../modules/gatk/mark_duplicates/main.nf'

workflow MARK_DUPLICATES {
    take:
    bam_input_ch   

    main:
    log.info "Starting MarkDuplicates Workflow..."

    // **Step 1: Determine Input BAMs**
    filtered_bams_ch = Channel.empty()

    if (params.input_bam) {
        log.info "Using user-provided BAM files for MarkDuplicates..."
        bam_channel = Channel.fromPath(params.input_bam)
        filtered_bams_ch = SAMTOOLS_FILTER_ORPHANS(bam_channel)
    } else {
        log.info "Using filtered BAMs from STAR alignment..."
        filtered_bams_ch = bam_input_ch 
    }

    // **Step 2: Mark Duplicates**
    dedup_bams_ch = GATK_MARK_DUPLICATES(filtered_bams_ch)

    // **Emit Outputs**
    emit:
    marked_bams_bai = dedup_bams_ch.marked_bams_bai
	marked_bams_bai_metrics = dedup_bams_ch.marked_bams_bai_metrics
}
