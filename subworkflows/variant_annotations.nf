nextflow.enable.dsl = 2

// Import annotation modules

include { VARIANT_ANNOTATION } from '../subworkflows/snpeff_annotations.nf'
include { VEP_ANNOTATION_WORKFLOW as COMBINED_ANNOTATE } from '../subworkflows/ensemblvep.nf'
include { VEP_ANNOTATION_WORKFLOW } from '../subworkflows/ensemblvep.nf'
include { COLLECT_VARIANT_CALLING_METRICS } from '../modules/gatk/variant_metrics/main.nf'


workflow ANNOTATE {
    take:
    vcf          
    tools
	snpeff_jar
    snpeff_db
    snpeff_config
	genomedb
    genome_assembly
    species
    vep_version
    vep_cache
	vep_plugins
	dbsnp_vcf
    dbsnp_tbi
    ref_fasta
    ref_dict
	
	
	main:

    // Initialize empty channels
    annotated_vcfs = Channel.empty()
    annotation_reports = Channel.empty()
	ch_versions = Channel.empty()
	uncompressed_vcf = Channel.empty()
    
	
//	log.info "Starting variant annotation workflow..."

    // Case 1: Run SnpEff (if 'snpeff' or 'combine' is selected)
    if (tools.contains('snpeff') || tools.contains('combine')) {
//        log.info "Running SnpEff annotation..."
		
		VARIANT_ANNOTATION(
		vcf,  
		snpeff_jar,
		snpeff_config,
		snpeff_db,
		genomedb
	)
	
	annotated_vcfs = annotated_vcfs.mix(VARIANT_ANNOTATION.out.final_annotated_variants)
	uncompressed_vcf = uncompressed_vcf.mix(VARIANT_ANNOTATION.out.annotated_variants)
    annotation_reports = annotation_reports.mix(VARIANT_ANNOTATION.out.annotated_summary)
	ch_versions = ch_versions.mix(VARIANT_ANNOTATION.out.versions)
	
	}
	
	if (tools.contains('combine')) {
        vcf_ann_for_merge = VARIANT_ANNOTATION.out.final_annotated_variants.map{ sample_id, vcf, tbi -> [sample_id, vcf, tbi] }
        COMBINED_ANNOTATE (
            vcf_ann_for_merge,
            vep_cache,
			genome_assembly,
			vep_version,
			species,
			vep_plugins
        )
        annotated_vcfs  = annotated_vcfs.mix( COMBINED_ANNOTATE.out.final_vep_annotated_vcf)
		uncompressed_vcf = uncompressed_vcf.mix(COMBINED_ANNOTATE.out.annotated_variants)
        annotation_reports  = annotation_reports.mix(COMBINED_ANNOTATE.out.annotated_summary)
		ch_versions = ch_versions.mix(COMBINED_ANNOTATE.out.versions)
		
		
    }
	
	if (tools.contains('vep')) {
        VEP_ANNOTATION_WORKFLOW(
        vcf,
		vep_cache,
        genome_assembly,
		vep_version,
		species,
		vep_plugins
       
    )
        annotated_vcfs  = annotated_vcfs.mix(VEP_ANNOTATION_WORKFLOW.out.final_vep_annotated_vcf)
		uncompressed_vcf = uncompressed_vcf.mix(VEP_ANNOTATION_WORKFLOW.out.annotated_variants)
        annotation_reports  = annotation_reports.mix(VEP_ANNOTATION_WORKFLOW.out.annotated_summary)
		ch_versions = ch_versions.mix(VEP_ANNOTATION_WORKFLOW.out.versions)
		
        
    }
	
	COLLECT_VARIANT_CALLING_METRICS(
        annotated_vcfs,   // tuple: val(meta), path(vcf.gz), path(.tbi)
        dbsnp_vcf, dbsnp_tbi,
        ref_fasta, ref_dict
    )
    metrics_tables = COLLECT_VARIANT_CALLING_METRICS.out.metrics
    ch_versions    = ch_versions.mix(COLLECT_VARIANT_CALLING_METRICS.out.versions)
	
	emit: 
		final_vcf_annotated     = annotated_vcfs
		uncompressed_vcf_annotated = uncompressed_vcf
		metrics                    = metrics_tables
        reports_html     		= annotation_reports
		versions 				= ch_versions
		
    }