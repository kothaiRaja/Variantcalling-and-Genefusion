nextflow.enable.dsl = 2

// Import subworkflows
include { PREPROCESSING } from './subworkflows/preprocessing.nf'
include { STAR_ALIGN } from './subworkflows/star_align.nf'
include { INTERVAL_PROCESSING } from './subworkflows/interval_processing.nf'
include { MARK_DUPLICATES } from './subworkflows/markduplicates.nf'
include { SPLIT_MERGE_BAMS } from './subworkflows/split_merge_bams.nf'
include { BASE_RECALIBRATION } from './subworkflows/base_recalibration.nf'
include { VARIANT_CALLING } from './subworkflows/variant_calling.nf'
include { ANNOTATE } from './subworkflows/variant_annotations.nf'








workflow {

    // Ensure only one of the specialized workflows is selected
    if (params.only_qc && params.skip_star) { 
        error "You cannot enable both only_qc and skip_star at the same time. Set one to false."
    }

    //============================== PREPROCESSING / TRIMMING ===============================//

    trimmed_reads_ch = Channel.empty() // Initialize empty channel

    if (params.input_reads) {
        log.info "Using provided input reads for STAR Alignment..."
        trimmed_reads_ch = Channel.fromPath(params.input_reads)
    } else {
        log.info "No input reads provided. Running Preprocessing..."
        PREPROCESSING(params.samplesheet)
        
        trimmed_reads_ch  = PREPROCESSING.out.trimmed_reads
        fastp_reports_ch  = PREPROCESSING.out.fastp_reports
        qc_results_ch     = PREPROCESSING.out.qc_reports
        multiqc_quality   = PREPROCESSING.out.multiqc
    }

    // If user wants QC-only mode, exit after preprocessing
    if (params.only_qc) {
        log.info("QC completed. Exiting pipeline...")
        return
    }

    //==================================== STAR ALIGNMENT ====================================//

    if (!params.skip_star) {
        log.info "Running STAR Alignment..."
        STAR_ALIGN(trimmed_reads_ch, params.star_genome_index, params.gtf_annotation)
		// Capture Outputs from STAR Alignment
        star_bam_ch         = STAR_ALIGN.out.bam_sorted
        chimeric_reads_ch   = STAR_ALIGN.out.chimeric_reads
        flagstats_ch        = STAR_ALIGN.out.flagstats
        align_stats_ch      = STAR_ALIGN.out.align_stats
        star_logs_ch        = STAR_ALIGN.out.star_logs
		filtered_bams_ch 	= STAR_ALIGN.out.filtered_bams
    } else {
        log.info "Skipping STAR Alignment as per user request."
    }
	
	//=====================================Intervals processing============================//
	
	// **Step 1: Convert BED file to Channel**
    ch_exon_bed = Channel
        .fromPath(params.exons_BED)
        .map { file_path -> tuple(file_path.baseName, file_path) }
        .view { meta, file -> println "  Exons BED Tuple: ${meta}, File: ${file}" }
		
	log.info "Starting Interval Processing..."
    INTERVAL_PROCESSING(ch_exon_bed, params.reference_genome, params.reference_genome_dict)

    // Capture Outputs from Interval Processing
    intervals_ch = INTERVAL_PROCESSING.out.intervals

    // **View the intervals to confirm output**
    intervals_ch.view { file -> "Generated Interval: ${file}" }
	
	//===================================Markduplicates==============================//
	
	log.info "Starting MarkDuplicates Processing..."

    MARK_DUPLICATES(filtered_bams_ch)

    // Capture Outputs
    dedup_bam_ch = MARK_DUPLICATES.out.marked_bams_bai
	
	// ==================== SPLIT & MERGE BAMs ==================== //

	log.info "Running Split & Merge BAMs..."

	SPLIT_MERGE_BAMS(
		dedup_bam_ch,        // BAMs after MarkDuplicates
		intervals_ch,		// Scattered intervals from Interval Processing
		params.reference_genome,
		params.reference_genome_index,
		params.reference_genome_dict
		
	)

	// Set the final BAMs channel
	final_bams_ch = SPLIT_MERGE_BAMS.out.merged_calmd_bams
	
	//======================BASE_RECALIBRATION=========================//
	
	if (!params.skip_base_recalibration) {
    log.info "Running Base Recalibration..."

    BASE_RECALIBRATION(
        final_bams_ch,    // BAMs after merging & CALMD
        intervals_ch,		// Scattered intervals from Interval Processing
    	params.reference_genome,
		params.reference_genome_index,
		params.reference_genome_dict,
		params.merged_vcf,
		params.merged_vcf_index)

    recalibrated_bams_ch = BASE_RECALIBRATION.out.recalibrated_bams
} else {
    log.info "Skipping Base Recalibration. Using previous step's BAMs..."
    recalibrated_bams_ch = final_bams_ch
}
	
	// =========== Step 6: Variant Calling =========== //
    log.info "🔹 Running Variant Calling..."
    VARIANT_CALLING(
        recalibrated_bams_ch,            // BAMs after splitting & merging
        params.reference_genome,
        params.reference_genome_index,
        params.reference_genome_dict,
        params.merged_vcf,
        params.merged_vcf_index
    )

    // Capture Variant Calling Outputs
    filtered_vcf_ch = VARIANT_CALLING.out.final_variants
    selected_snps_ch = VARIANT_CALLING.out.selected_snps
    selected_indels_ch = VARIANT_CALLING.out.selected_indels
	selected_variants_ch = VARIANT_CALLING.out.selected_variants

    log.info "✅ Variant Calling Completed!"
	
	// ===================== Step 7: Variant Annotation ===================== //
	log.info "🔹 Running Variant Annotation..."

	annotation_results = ANNOTATE(
        filtered_vcf_ch,
        params.annotation_tools,
		params.snpeff_jar,
        params.snpeff_db,
		params.snpeff_config,
        params.genomedb,
        params.genome_assembly,
        params.species,
        params.cache_version,
        params.vep_cache_dir,
		params.clinvar,
		params.clinvartbi
    )

    log.info "Variant annotation complete!"

	final_annotated_vcf = ANNOTATE.out.final_vcf_annotated
	html_report 		= ANNOTATE.out.reports_html
	
}
   
	
	
    
