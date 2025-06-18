nextflow.enable.dsl = 2

// Import subworkflows
include { PREPROCESSING } from '../subworkflows/preprocessing.nf'
include { STAR_ALIGN } from '../subworkflows/star_align.nf'
include { INTERVAL_PROCESSING } from '../subworkflows/interval_processing.nf'
include { MARK_DUPLICATES } from '../subworkflows/markduplicates.nf'
include { SPLIT_MERGE_BAMS } from '../subworkflows/split_merge_bams.nf'
include { BASE_RECALIBRATION } from '../subworkflows/base_recalibration.nf'
include { VARIANT_CALLING } from '../subworkflows/variant_calling.nf'
include { ANNOTATE } from '../subworkflows/variant_annotations.nf'
include { GENE_FUSION } from '../subworkflows/gene_fusion.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nfcore/software_versions/main.nf'
include { MultiQC } from '../modules/multiqc_quality/main.nf'
include { GATK_VCF_TO_TABLE } from '../modules/gatk/vcf2table/main.nf'
include { MAF_ANALYSIS } from '../subworkflows/maf_analysis.nf'












workflow RNA_VARIANT_CALLING_GENE_FUSION {

    ch_versions = Channel.empty()
	reports_ch = Channel.empty()


    // ============================ VALIDATION ============================
    if (params.only_qc && params.only_star) {
        error " Cannot run both 'only_qc' and 'only_star'. Please set only one to true."
    }

    // ============================ ONLY QC MODE ============================
    if (params.only_qc) {
        log.info " Running QC-only mode..."

        PREPROCESSING(params.samplesheet, params.dump_script)

        reports_ch       = PREPROCESSING.out.reports
        multiqc_quality  = PREPROCESSING.out.multiqc
        ch_versions      = ch_versions.mix(PREPROCESSING.out.versions)

        log.info(" QC completed. Exiting pipeline.")
        return
    }

    // ============================ ONLY STAR MODE ============================
    if (params.only_star) {
        log.info " Running ONLY STAR alignment..."

        PREPROCESSING(params.samplesheet, params.dump_script)

        STAR_ALIGN(
            PREPROCESSING.out.trimmed_reads,
            params.star_genome_index,
            params.gtf_annotation,
            null,
            null
        )

        STAR_ALIGN.out.bam_sorted.view { " STAR sorted BAM: $it" }

        return
    }

    // ============================ FULL PIPELINE MODE ============================
    log.info " Starting Preprocessing using sample sheet..."

    PREPROCESSING(params.samplesheet, params.dump_script)

    reports_ch 		= reports_ch.mix(PREPROCESSING.out.qc_results.collect { it[1] }.ifEmpty([]))
	reports_ch    = reports_ch.mix(PREPROCESSING.out.fastp_reports.collect { it[1] }.ifEmpty([]))
	trimmed_reads_ch = PREPROCESSING.out.trimmed_reads
    QC_reports_ch       = PREPROCESSING.out.reports
    multiqc_quality    = PREPROCESSING.out.multiqc
    ch_versions        = ch_versions.mix(PREPROCESSING.out.versions)

    // ========================== STAR ALIGNMENT LOGIC ==========================
    if (!params.skip_star) {
        log.info " Running STAR Alignment..."

        def aligned_bam_samplesheet = params.aligned_bam_samplesheet ? file(params.aligned_bam_samplesheet) : null
        def aligned_bam_folder      = params.aligned_bam_folder      ? file(params.aligned_bam_folder) : null

        STAR_ALIGN(
            trimmed_reads_ch, 
            params.star_genome_index, 
            params.gtf_annotation, 
            aligned_bam_samplesheet, 
            aligned_bam_folder
        )

        star_bam_ch        = STAR_ALIGN.out.bam_sorted
        chimeric_reads_ch  = STAR_ALIGN.out.chimeric_reads
        filtered_bams_ch   = STAR_ALIGN.out.filtered_bams
		flagstats_ch 	   = STAR_ALIGN.out.flagstats
		align_stats_ch     = STAR_ALIGN.out.align_stats
		star_logs_ch       = STAR_ALIGN.out.star_logs
        ch_versions        = ch_versions.mix(STAR_ALIGN.out.versions)
		reports_ch = reports_ch.mix(flagstats_ch.collect { it[2] }.ifEmpty([]))
							.mix(align_stats_ch.collect { it[2] }.ifEmpty([]))
							.mix(star_logs_ch.collect { it[2] }.ifEmpty([]))

    } else {
        log.info " Skipping STAR alignment as per user request."
		
		def aligned_bam_folder = params.aligned_bam_folder ? file(params.aligned_bam_folder) : null

    if (!aligned_bam_folder) {
        error "You must set 'aligned_bam_folder' when skipping STAR alignment."
    }

    filtered_bams_ch = Channel
        .fromPath("${aligned_bam_folder}/*.bam")
        .map { bam_file ->
            def sample_id = bam_file.baseName.replaceAll(/\.bam$/, "")
            tuple(sample_id, bam_file)
        }

        star_bam_ch        = Channel.empty()
        chimeric_reads_ch  = Channel.empty()
        flagstats_ch       = Channel.empty()
        align_stats_ch     = Channel.empty()
        star_logs_ch       = Channel.empty()
        
    }



	
	//=====================================Intervals processing============================//
	
	// Step 1: Convert BED file to Channel
    ch_exon_bed = Channel
        .fromPath(params.exons_BED)
        .map { file_path -> tuple(file_path.baseName, file_path) }
        .view { meta, file -> println "  Exons BED Tuple: ${meta}, File: ${file}" }
		
	log.info "Starting Interval Processing..."
    INTERVAL_PROCESSING(ch_exon_bed, params.reference_genome, params.reference_genome_index, params.reference_genome_dict)

    // Capture Outputs from Interval Processing
    intervals_ch = INTERVAL_PROCESSING.out.intervals
	ch_versions        = ch_versions.mix(INTERVAL_PROCESSING.out.versions)

    // View the intervals to confirm output
    intervals_ch.view { file -> "Generated Interval: ${file}" }
	
	//===================================Markduplicates==============================//
	
	log.info "Starting MarkDuplicates Processing..."

    MARK_DUPLICATES(filtered_bams_ch)

    // Capture Outputs
    dedup_bam_ch = MARK_DUPLICATES.out.marked_bams_bai
	reports_ch = reports_ch.mix(MARK_DUPLICATES.out.marked_bams_bai_metrics.ifEmpty([]))

	
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
	final_bams_ch 	   = SPLIT_MERGE_BAMS.out.merged_calmd_bams
	ch_versions        = ch_versions.mix(SPLIT_MERGE_BAMS.out.versions)
	
	
	//======================BASE_RECALIBRATION=========================//
	
	
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
	ch_versions        = ch_versions.mix(BASE_RECALIBRATION.out.versions)

	
	// =========== Step 6: Variant Calling =========== //
    log.info " Running Variant Calling..."
    VARIANT_CALLING(
        recalibrated_bams_ch,            // BAMs after splitting & merging
		intervals_ch,
        params.reference_genome,
        params.reference_genome_index,
        params.reference_genome_dict,
        params.merged_vcf,
        params.merged_vcf_index
    )

    // Capture Variant Calling Outputs
    filtered_vcf_ch = VARIANT_CALLING.out.final_variants
	ch_versions        = ch_versions.mix(VARIANT_CALLING.out.versions)
	reports_ch = reports_ch
    .mix(
        (VARIANT_CALLING.out.bcftools_stats
            .collect { it[1] }          
        ).ifEmpty([])                   
    )
    .mix(
        (VARIANT_CALLING.out.bcftools_query
            .collect { it[1] }
        ).ifEmpty([])
    )


    log.info " Variant Calling Completed!"
	
	// ===================== Step 7: Variant Annotation ===================== //
	log.info " Running Variant Annotation..."

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
		params.vep_plugins_dir
		
    )

    log.info "Variant annotation complete!"

	final_annotated_vcf = ANNOTATE.out.final_vcf_annotated
	uncompressed_annotated_vcf = ANNOTATE.out.uncompressed_vcf_annotated
	report_ch   		= reports_ch.mix(ANNOTATE.out.reports_html.ifEmpty([]))
	ch_versions        = ch_versions.mix(ANNOTATE.out.versions)
	
	//=====================Create Table============================//
	
	GATK_VCF_TO_TABLE(final_annotated_vcf)
	
	table_ch 			= GATK_VCF_TO_TABLE.out.vcf_table
	ch_versions        = ch_versions.mix(GATK_VCF_TO_TABLE.out.versions)
	
	//=====================Maftools Visualisation======================//
	
	if (params.maftools) {
        log.info "Running MAF_ANALYSIS Subworkflow..."

        MAF_ANALYSIS(
            uncompressed_annotated_vcf,
            params.reference_genome,
            params.vep_cache_dir,
            params.rscript
        )

        // Capture outputs
        maf_reports_ch = MAF_ANALYSIS.out.maf_plots
        ch_versions = ch_versions.mix(MAF_ANALYSIS.out.versions)
    } else {
        log.info "Skipping MAF_ANALYSIS Subworkflow (params.maftools = false)"
        maf_reports_ch = Channel.empty()
    }
	
	
	
	//================== Step 8: Run Gene Fusion Analysis on STAR Chimeric Reads ================//

if (params.run_fusion && !params.skip_star) {
    log.info " Running Gene Fusion Analysis..."

    fusions = GENE_FUSION(
        star_bam_ch,  
        params.reference_genome,
        params.gtf_annotation,
        params.arriba_blacklist,
        params.arriba_known_fusions,
		params.cytobands,
		params.protein_domains
		
    )

    ARRIBA_fusion_ch = GENE_FUSION.out.fusion_results
    reports_ch = report_ch.mix(GENE_FUSION.out.fusion_visualizations.collect { it[1] }.ifEmpty([]))
    ch_versions = ch_versions.mix(GENE_FUSION.out.versions)

} else if (params.run_fusion && params.skip_star) {
    log.warn " STAR alignment was skipped. Gene fusion analysis requires STAR chimeric reads."
} else {
    log.info " Skipping gene fusion detection as per configuration (params.run_fusion = false)"
}


 //=======================================Version collection==============================//
 
 ch_version_yaml = Channel.empty()
 CUSTOM_DUMPSOFTWAREVERSIONS(
		ch_versions.unique().collectFile(name: "software_versions_input.yml"),
		params.dump_script
	)
report_ch = reports_ch.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.ifEmpty([]))

 
 // ======================= Final Report Collection for MultiQC ======================= //

// Collect everything
final_reports_ch = Channel
    .empty()
    .mix(reports_ch)
    .mix(ANNOTATE.out.reports_html.ifEmpty([]))
    .unique { it.name }   
    .collect()

// Run MultiQC
multiqc_quality = MultiQC(final_reports_ch)




	
}
   
	
	
    
