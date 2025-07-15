nextflow.enable.dsl = 2

// Import subworkflows for reference files 
include { BUILD_REFERENCES } from '../subworkflows/build_references.nf'

// Import subworkflows
include { INPUT_PAIRED_READS } from '../subworkflows/input.nf'
include { PREPROCESSING } from '../subworkflows/preprocessing.nf'
include { STAR_ALIGN } from '../subworkflows/star_align.nf'
include { BED_TO_INTERVAL_LIST } from '../modules/Intervals/bed_to_intervals/main.nf'
include { SCATTER_INTERVAL_LIST } from '../modules/Intervals/scattered_intervals/main.nf'
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

	// Convert params to channels
	ch_input = file(params.samplesheet)
	ch_genomedb = params.genomedb ? Channel.value(params.genomedb) : Channel.empty()
	ch_annotation_tools = params.annotation_tools ? Channel.value(params.annotation_tools) : Channel.empty()
	ch_rscript = params.rscript ? Channel.fromPath(params.rscript) : Channel.empty()
	ch_ver_script = params.dump_script ? Channel.fromPath(params.dump_script) : Channel.empty()
	
	ch_arriba_blacklist = params.arriba_blacklist ? Channel.fromPath(params.arriba_blacklist) : Channel.empty()
	ch_arriba_known_fusions = params.arriba_known_fusions ? Channel.fromPath(params.arriba_known_fusions) : Channel.empty()

	ch_genome_assembly = params.genome_assembly ? Channel.value(params.genome_assembly) : Channel.empty()
	ch_species = params.species ? Channel.value(params.species) : Channel.empty()
	ch_cache_version = params.cache_version ? Channel.value(params.cache_version) : Channel.empty()

	ch_known_snps_vcf = params.known_snps_vcf ? Channel.fromPath(params.known_snps_vcf, checkIfExists: true) : Channel.empty()
	ch_known_snps_index = params.known_snps_vcf_index ? Channel.fromPath("${params.known_snps_vcf_index}", checkIfExists: true) : Channel.empty()

	ch_known_indels_vcf = params.known_indels_vcf ? Channel.fromPath(params.known_indels_vcf, checkIfExists: true) : Channel.empty()
	ch_known_indels_index = params.known_indels_vcf_index ? Channel.fromPath("${params.known_indels_vcf_index}", checkIfExists: true) : Channel.empty()

	ch_multiqc_config = Channel.fromPath(file("$projectDir/assets/multiqc_config.yml", checkIfExists: true))

	ch_versions = Channel.empty()
	reports_ch = Channel.empty()
	
	//===========================Building References======================
	BUILD_REFERENCES(ch_genomedb)

	reference_genome_ch = BUILD_REFERENCES.out.reference_genome.collect()
	reference_genome_index_ch = BUILD_REFERENCES.out.reference_genome_index.collect()
	reference_genome_dict_ch = BUILD_REFERENCES.out.reference_genome_dict.collect()

	gtf_ch = BUILD_REFERENCES.out.gtf_annotation.collect()
	exons_bed_ch = BUILD_REFERENCES.out.exons_BED
	ch_exons_bed_tuple = exons_bed_ch.map { bed -> tuple(bed.baseName, bed) }


	star_index_ch = BUILD_REFERENCES.out.star_genome_index.collect()

	snpeff_jar_ch = BUILD_REFERENCES.out.snpeff_jar
	snpeff_config_ch = BUILD_REFERENCES.out.snpeff_config
	snpeff_db_dir_ch = BUILD_REFERENCES.out.snpeff_db_dir

	arriba_dir_ch = BUILD_REFERENCES.out.arriba_dir

	vep_cache_ch = BUILD_REFERENCES.out.vep_cache
	vep_plugins_ch = BUILD_REFERENCES.out.vep_plugins

	//========================Input channel ==================================
	INPUT_PAIRED_READS(ch_input)
	reads_ch = INPUT_PAIRED_READS.out.paired_reads
//	reads_ch.view { "READS_CH: $it" }

	// ============================ PREPROCESSING ============================
	PREPROCESSING(reads_ch)

	reports_ch = reports_ch.mix(PREPROCESSING.out.reports.ifEmpty([]))
	trimmed_reads_ch = PREPROCESSING.out.trimmed_reads
	ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

	// ========================== STAR ALIGNMENT LOGIC ==========================
//	log.info " Running STAR Alignment..."

	STAR_ALIGN(
		trimmed_reads_ch, 
		star_index_ch, 
		gtf_ch
	)

	star_bam_ch = STAR_ALIGN.out.bam_sorted
	chimeric_reads_ch = STAR_ALIGN.out.chimeric_reads
	chimeric_junction_ch = STAR_ALIGN.out.chimeric_junction
	filtered_bams_ch = STAR_ALIGN.out.filtered_bams
	flagstats_ch 	   = STAR_ALIGN.out.flagstats
		align_stats_ch     = STAR_ALIGN.out.align_stats
		star_logs_ch       = STAR_ALIGN.out.star_logs
        ch_versions        = ch_versions.mix(STAR_ALIGN.out.versions)
		reports_ch = reports_ch.mix(flagstats_ch.map { it[1] }.ifEmpty([]))
							.mix(align_stats_ch.map { it[1] }.ifEmpty([]))
							.mix(star_logs_ch.map { it[1] }.ifEmpty([]))
							

		
	// ===================== Intervals Processing ===================== //
//		log.info "Starting Interval Processing..."
		
		// **Step 1: Convert BED to Interval List**
		BED_TO_INTERVAL_LIST(ch_exons_bed_tuple, reference_genome_ch, reference_genome_dict_ch)
		def interval_list = BED_TO_INTERVAL_LIST.out.interval_list
		
		

		
		

	// Step 2: Scatter if enabled, else use the full list
    def scattered_intervals_ch = Channel.empty()

    if (params.scatterintervals) {
        log.info " Scattering intervals for parallel execution..." 
		SCATTER_INTERVAL_LIST(interval_list, reference_genome_ch, reference_genome_index_ch, reference_genome_dict_ch)
		scattered_intervals_ch = SCATTER_INTERVAL_LIST.out.scattered_intervals
								.map{ _meta, bed -> [bed] }.collect()
        
        scattered_intervals_ch.view { file -> " Scattered interval: $file" }
        ch_versions = ch_versions.mix(SCATTER_INTERVAL_LIST.out.versions)
    } else {
        log.info " Using full interval list without scattering..."
        scattered_intervals_ch = interval_list.map { _meta, bed -> bed }  

    }
	
	//===================================Markduplicates==============================//
	
	

    MARK_DUPLICATES(filtered_bams_ch)

    // Capture Outputs
    dedup_bam_ch = MARK_DUPLICATES.out.marked_bams_bai
	reports_ch = reports_ch.mix(MARK_DUPLICATES.out.marked_bams_bai_metrics.ifEmpty([]))
	
	// ==================== SPLIT & MERGE BAMs ==================== //

	

	SPLIT_MERGE_BAMS(
		dedup_bam_ch,        
		scattered_intervals_ch,		
		reference_genome_ch,
		reference_genome_index_ch,
		reference_genome_dict_ch
		
	)

	// Set the final BAMs channel
	final_bams_ch 	   = SPLIT_MERGE_BAMS.out.merged_calmd_bams
	ch_versions        = ch_versions.mix(SPLIT_MERGE_BAMS.out.versions)
	
	//==========================BASE_RECALIBRATION======================//
	
	def intervals_ch_recalib = interval_list.map { meta_id, interval_path -> tuple(meta_id, interval_path) }
	ch_recalib_input = final_bams_ch
    .combine(intervals_ch_recalib)
    .map { bam_meta, bam, bai, interval_meta, interval_file ->
        def sample_id    = bam_meta.id
        def strandedness = bam_meta.strandedness
        def meta         = [id: sample_id, strandedness: strandedness]
        tuple(meta, bam, bai, interval_file)
    }
   

	

//	  ch_recalib_input.view {"recalib_bams : $it "}



	BASE_RECALIBRATION(
        final_bams_ch, 
		ch_recalib_input,		
        intervals_ch_recalib,		
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
    
	
	VARIANT_CALLING(
        recalibrated_bams_ch,            
		scattered_intervals_ch,
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
            .map { it[1] }          
        ).ifEmpty([])                   
    )
    .mix(
        (VARIANT_CALLING.out.bcftools_query
            .map { it[1] }
        ).ifEmpty([])
    )
	
	// ===================== Step 7: Variant Annotation ===================== //
	

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

//    log.info "Variant annotation complete!"

	final_annotated_vcf = ANNOTATE.out.final_vcf_annotated
	uncompressed_annotated_vcf = ANNOTATE.out.uncompressed_vcf_annotated
	reports_ch   		= reports_ch.mix(ANNOTATE.out.reports_html.ifEmpty([]))
	ch_versions        = ch_versions.mix(ANNOTATE.out.versions)
	
	//=====================Maftools Visualisation======================//
	
	if (params.maftools) {
//        log.info "Running MAF_ANALYSIS Subworkflow..."

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
//    log.info " Running Gene Fusion Analysis..."

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

// Combine all report sources into one channel
all_reports_ch = Channel
    .empty()
    .mix(reports_ch)
    .mix(ch_multiqc_config.ifEmpty([]))
    .filter { it instanceof Path && it.exists() }
    .unique { it.name }
    .view { it -> " MultiQC input file: ${it?.getClass()?.getSimpleName()} - ${it}" }
    .collect()
    
// Run MultiQC
multiqc_quality = MultiQC(all_reports_ch)






	
		
} 
