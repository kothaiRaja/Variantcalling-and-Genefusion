nextflow.enable.dsl = 2
// Import subworkflows for reference files 

include { DOWNLOAD_REFERENCE_GENOME } from '../subworkflows/download_genome.nf'
include { DOWNLOAD_GTF_ANNOTATION } from '../subworkflows/download_gtf.nf'
include { DOWNLOAD_AND_PREPARE_VARIANT_VCFS } from '../subworkflows/prepare_vcf.nf'
include { BUILD_STAR_INDEX } from '../subworkflows/star_index.nf'
include { SNPEFF_SETUP } from '../subworkflows/snpeff_setup.nf'
include { ARRIBA_SETUP} from '../subworkflows/arriba.nf'
include { VEP_SETUP } from '../subworkflows/VEP_setup.nf'









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
	
	// ============================ DOWNLOAD OR USE GENOME + GTF ============================
   
    genome_refs = DOWNLOAD_REFERENCE_GENOME()
    reference_genome_ch       = genome_refs.genome
    reference_genome_index_ch = genome_refs.genome_index
    reference_genome_dict_ch  = genome_refs.genome_dict
	
	// ============================ DOWNLOAD OR USE GTF + EXONS BED ============================

	gtf_outputs = DOWNLOAD_GTF_ANNOTATION()
	gtf_ch        = gtf_outputs.gtf
	exons_bed_ch  = gtf_outputs.exons_bed
	
	//=========================PREPARING KNOWN_VCFs=====================
	
	prepare_variants_output = DOWNLOAD_AND_PREPARE_VARIANT_VCFS()
	known_variants_ch = prepare_variants_output.merged_vcf 
	known_variants_index = prepare_variants_output.merged_vcf_index
	
	//========================Creating STAR index=========================
	BUILD_STAR_INDEX(reference_genome_ch, gtf_ch)
	
	star_index_ch = BUILD_STAR_INDEX.out.star_index
	
	//======================SNPEFF TOOL===============================
	
	SNPEFF_SETUP(params.genomedb)

	snpeff_jar_ch    = SNPEFF_SETUP.out.snpeff_jar
	snpeff_config_ch = SNPEFF_SETUP.out.snpeff_config
	snpeff_db_dir_ch = SNPEFF_SETUP.out.snpeff_db_dir
	
	//====================== Run the arriba setup
	ARRIBA_SETUP()


	arriba_dir_ch = ARRIBA_SETUP.out.arriba_dir
	
	//==================VEP_SETUP======================================
	vep_outputs = VEP_SETUP()  

	vep_cache_ch   = vep_outputs.vep_cache
	vep_plugins_ch = vep_outputs.vep_plugins
	
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
        log.info " Running STAR Alignment..."

        
        STAR_ALIGN(
            trimmed_reads_ch, 
            star_index_ch, 
            gtf_ch
        )

        star_bam_ch        = STAR_ALIGN.out.bam_sorted
        chimeric_reads_ch  = STAR_ALIGN.out.chimeric_reads
		chimeric_junction_ch = STAR_ALIGN.out.chimeric_junction
        filtered_bams_ch   = STAR_ALIGN.out.filtered_bams
		flagstats_ch 	   = STAR_ALIGN.out.flagstats
		align_stats_ch     = STAR_ALIGN.out.align_stats
		star_logs_ch       = STAR_ALIGN.out.star_logs
        ch_versions        = ch_versions.mix(STAR_ALIGN.out.versions)
		reports_ch = reports_ch.mix(flagstats_ch.collect { it[2] }.ifEmpty([]))
							.mix(align_stats_ch.collect { it[2] }.ifEmpty([]))
							.mix(star_logs_ch.collect { it[2] }.ifEmpty([]))
							
		
							
							
	// ===================== Intervals Processing ===================== //
		log.info "Starting Interval Processing..."

	// Use exons_bed_ch that was generated earlier
	ch_exon_bed = exons_bed_ch.map { file -> tuple(file.baseName, file) }

	INTERVAL_PROCESSING(
		ch_exon_bed,
		reference_genome_ch,
		reference_genome_index_ch,
		reference_genome_dict_ch
	)

	// Capture Outputs from Interval Processing
	intervals_ch = INTERVAL_PROCESSING.out.intervals
	ch_versions  = ch_versions.mix(INTERVAL_PROCESSING.out.versions)

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
		dedup_bam_ch,        
		intervals_ch,		
		reference_genome_ch,
		reference_genome_index_ch,
		reference_genome_dict_ch
		
	)

	// Set the final BAMs channel
	final_bams_ch 	   = SPLIT_MERGE_BAMS.out.merged_calmd_bams
	ch_versions        = ch_versions.mix(SPLIT_MERGE_BAMS.out.versions)
	
	//======================BASE_RECALIBRATION=========================//
	
	
    log.info "Running Base Recalibration..."

    BASE_RECALIBRATION(
        final_bams_ch,    // BAMs after merging & CALMD
        intervals_ch,		// Scattered intervals from Interval Processing
    	reference_genome_ch,
		reference_genome_index_ch,
		reference_genome_dict_ch,
		known_variants_ch,
		known_variants_index)

    recalibrated_bams_ch = BASE_RECALIBRATION.out.recalibrated_bams
	ch_versions        = ch_versions.mix(BASE_RECALIBRATION.out.versions)
	
	// =========== Step 6: Variant Calling =========== //
    log.info " Running Variant Calling..."
    VARIANT_CALLING(
        recalibrated_bams_ch,            // BAMs after splitting & merging
		intervals_ch,
        reference_genome_ch,
		reference_genome_index_ch,
		reference_genome_dict_ch,
		known_variants_ch,
		known_variants_index
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
		snpeff_jar_ch,
		snpeff_db_dir_ch,
		snpeff_config_ch,
		params.genomedb,
		params.genome_assembly,
		params.species,
		params.cache_version,
		vep_cache_ch,
		vep_plugins_ch
		
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
            reference_genome_ch,
            vep_cache_ch,
            params.rscript
        )

        // Capture outputs
        maf_reports_ch = MAF_ANALYSIS.out.maf_plots
        ch_versions = ch_versions.mix(MAF_ANALYSIS.out.versions)
    } else {
        log.info "Skipping MAF_ANALYSIS Subworkflow (params.maftools = false)"
        maf_reports_ch = Channel.empty()
    }
	
	
	//================== Step 8: Run Gene Fusion Analysis Independently ================//

if (params.run_fusion) {
    log.info " Running Gene Fusion Analysis..."

    fusions = GENE_FUSION(
        star_bam_ch,  
        reference_genome_ch,
        gtf_ch,
        params.arriba_blacklist,
        params.arriba_known_fusions
		
		
    )

    ARRIBA_fusion_ch = GENE_FUSION.out.fusion_results
    reports_ch = report_ch.mix(GENE_FUSION.out.fusion_visualizations.collect { it[1] }.ifEmpty([]))
    ch_versions = ch_versions.mix(GENE_FUSION.out.versions)


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






	

	

   
	
	
    
