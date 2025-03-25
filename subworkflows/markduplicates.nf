nextflow.enable.dsl = 2



include { GATK_MARK_DUPLICATES } from '../modules/gatk/mark_duplicates/main.nf'



workflow MARK_DUPLICATES {
    take:
    bam_input_ch   

    main:
	
	ch_versions = Channel.empty()
    log.info "Starting MarkDuplicates Workflow..."

    

    // **Step 1: Mark Duplicates**
    dedup_bams = GATK_MARK_DUPLICATES(bam_input_ch)
	
	dedup_bams_ch = GATK_MARK_DUPLICATES.out.marked_bams_bai
	dedup_metrics_ch = GATK_MARK_DUPLICATES.out.marked_bams_bai_metrics
		.collect { it[2] }
		.ifEmpty([])
	ch_versions = ch_versions.mix(GATK_MARK_DUPLICATES.out.versions.first())
	
	

    // **Emit Outputs**
    emit:
    marked_bams_bai = dedup_bams_ch
	marked_bams_bai_metrics = dedup_metrics_ch
	versions  = ch_versions
}
