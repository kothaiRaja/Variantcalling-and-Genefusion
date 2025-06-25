nextflow.enable.dsl = 2

// Convert params to channels

ch_genomedb                = params.genomedb                ? Channel.value(params.genomedb) : Channel.empty()
ch_annotation_tools        = params.annotation_tools        ? Channel.value(params.annotation_tools) : Channel.empty()
ch_rscript                 = params.rscript                 ? Channel.fromPath(params.rscript) : Channel.empty()
ch_ver_script 			   = params.dump_script				? Channel.fromPath(params.dump_script) : Channel.empty()

ch_arriba_blacklist        = params.arriba_blacklist        ? Channel.fromPath(params.arriba_blacklist) : Channel.empty()
ch_arriba_known_fusions    = params.arriba_known_fusions    ? Channel.fromPath(params.arriba_known_fusions) : Channel.empty()

ch_genome_assembly         = params.genome_assembly         ? Channel.value(params.genome_assembly) : Channel.empty()
ch_species                 = params.species                 ? Channel.value(params.species) : Channel.empty()
ch_cache_version           = params.cache_version           ? Channel.value(params.cache_version) : Channel.empty()

ch_known_snps_vcf     	   = params.known_snps_vcf     ? Channel.fromPath(params.known_snps_vcf, checkIfExists: true) : Channel.empty()
ch_known_snps_index   	   = params.known_snps_vcf_index     ? Channel.fromPath("${params.known_snps_vcf_index}", checkIfExists: true) : Channel.empty()

ch_known_indels_vcf  	   = params.known_indels_vcf   ? Channel.fromPath(params.known_indels_vcf, checkIfExists: true) : Channel.empty()
ch_known_indels_index 	   = params.known_indels_vcf_index   ? Channel.fromPath("${params.known_indels_vcf_index}", checkIfExists: true) : Channel.empty()



// Import subworkflows for reference files 
include { BUILD_REFERENCES } from '../subworkflows/build_references.nf'







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
include { MULTIQC_WRAPPER } from '../subworkflows/multiqc.nf'
include { GATK_VCF_TO_TABLE } from '../modules/gatk/vcf2table/main.nf'
include { MAF_ANALYSIS } from '../subworkflows/maf_analysis.nf'












workflow RNA_VARIANT_CALLING_GENE_FUSION {

    ch_versions = Channel.empty()
	reports_ch = Channel.empty()
	
	//===========================Building References======================
	
	BUILD_REFERENCES(params.samplesheet, ch_genomedb)
	
	reference_genome_ch        = BUILD_REFERENCES.out.reference_genome
	reference_genome_index_ch  = BUILD_REFERENCES.out.reference_genome_index
	reference_genome_dict_ch   = BUILD_REFERENCES.out.reference_genome_dict

	gtf_ch        = BUILD_REFERENCES.out.gtf_annotation
	exons_bed_ch  = BUILD_REFERENCES.out.exons_BED

	star_index_ch = BUILD_REFERENCES.out.star_genome_index

	snpeff_jar_ch    = BUILD_REFERENCES.out.snpeff_jar
	snpeff_config_ch = BUILD_REFERENCES.out.snpeff_config
	snpeff_db_dir_ch = BUILD_REFERENCES.out.snpeff_db_dir

	arriba_dir_ch    = BUILD_REFERENCES.out.arriba_dir

	vep_cache_ch     = BUILD_REFERENCES.out.vep_cache
	vep_plugins_ch   = BUILD_REFERENCES.out.vep_plugins


	
	validated_reads_ch = BUILD_REFERENCES.out.validated_reads
	
	


	
	// ============================ FULL PIPELINE MODE ============================
   log.info " Starting Preprocessing using validated reads from build references..."

	PREPROCESSING(validated_reads_ch )

	reports_ch        = reports_ch.mix(PREPROCESSING.out.qc_results.collect { it[1] }.ifEmpty([]))
	reports_ch        = reports_ch.mix(PREPROCESSING.out.fastp_reports.collect { it[1] }.ifEmpty([]))
	trimmed_reads_ch  = PREPROCESSING.out.trimmed_reads
	ch_versions       = ch_versions.mix(PREPROCESSING.out.versions)

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
		ch_known_snps_vcf,
		ch_known_snps_index,
		ch_known_indels_vcf,
		ch_known_indels_index
		 )

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
		ch_known_snps_vcf,
		ch_known_snps_index
		
		
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
		ch_genomedb,
		ch_genome_assembly,
		ch_species,
		ch_cache_version,
		vep_cache_ch,
		vep_plugins_ch
		
    )

    log.info "Variant annotation complete!"

	final_annotated_vcf = ANNOTATE.out.final_vcf_annotated
	uncompressed_annotated_vcf = ANNOTATE.out.uncompressed_vcf_annotated
	reports_ch   		= reports_ch.mix(ANNOTATE.out.reports_html.ifEmpty([]))
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
            ch_rscript
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
        ch_arriba_blacklist,
        ch_arriba_known_fusions
		
		
    )

    ARRIBA_fusion_ch = GENE_FUSION.out.fusion_results
    reports_ch = reports_ch.mix(GENE_FUSION.out.fusion_visualizations.collect { it[1] }.ifEmpty([]))
    ch_versions = ch_versions.mix(GENE_FUSION.out.versions)


} else {
    log.info " Skipping gene fusion detection as per configuration (params.run_fusion = false)"
}

//=======================================Version collection==============================//
 
 ch_version_yaml = Channel.empty()
 CUSTOM_DUMPSOFTWAREVERSIONS(
		ch_versions.unique().collectFile(name: "software_versions_input.yml"),
		ch_ver_script
	)
reports_ch = reports_ch.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.ifEmpty([]))

 
 // ======================= Final Report Collection for MultiQC ======================= //

// Collect everything



// Run MultiQC wrapper subworkflow
multiqc_result = MULTIQC_WRAPPER(
    PREPROCESSING.out.qc_results.collect { it[1] }.ifEmpty([]),
    PREPROCESSING.out.fastp_reports.collect { it[1] }.ifEmpty([]),
    STAR_ALIGN.out.star_logs.collect { it[2] }.ifEmpty([]),
    STAR_ALIGN.out.flagstats.collect { it[2] }.ifEmpty([]),
    MARK_DUPLICATES.out.marked_bams_bai_metrics.ifEmpty([]),
    VARIANT_CALLING.out.bcftools_stats.collect { it[1] }.ifEmpty([]),
    VARIANT_CALLING.out.bcftools_query.collect { it[1] }.ifEmpty([]),
    ANNOTATE.out.reports_html.ifEmpty([]),
    GENE_FUSION.out.fusion_visualizations.collect { it[1] }.ifEmpty([]),
    maf_reports_ch.ifEmpty([]),
    CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.ifEmpty([])
)





	

    
}






	

	

   
	
	
    
