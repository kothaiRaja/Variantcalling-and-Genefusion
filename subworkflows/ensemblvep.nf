nextflow.enable.dsl=2

// Include modules
include { ANNOTATEVARIANTS_VEP } from '../modules/ensemblvep/main.nf'
include { BGZIP_TABIX_ANNOTATIONS } from '../modules/tabix/bziptabix/main.nf'

workflow VEP_ANNOTATION_WORKFLOW {

    take:
    vcf_files_ch 
	vep_cache
	vep_plugin
	genome_assembly
	vep_version
	species
	
    main:
    log.info "Starting Ensembl VEP annotation..."
	
	ch_versions = Channel.empty()

    // Step 1: Annotate VCF using Ensembl VEP
    annotated_vcf = ANNOTATEVARIANTS_VEP(
        vcf_files_ch,
		vep_cache,
		vep_plugin,
        genome_assembly,
		vep_version,
		species
       
    )
	
	annotated_vcf_ch = ANNOTATEVARIANTS_VEP.out.annotated_vcf
	annotated_summary_ch = ANNOTATEVARIANTS_VEP.out.summary
	ch_versions = ch_versions.mix(ANNOTATEVARIANTS_VEP.out.versions) 

    log.info "Compressing and indexing annotated VCFs..."

    // Step 2: Compress and index annotated VCFs
    compressed_vcf = BGZIP_TABIX_ANNOTATIONS(annotated_vcf_ch)
	
	compressed_vcf_ch = BGZIP_TABIX_ANNOTATIONS.out.compressed_indexed
	ch_versions = ch_versions.mix(BGZIP_TABIX_ANNOTATIONS.out.versions)
	
    emit:
	annotated_variants = annotated_vcf_ch
	annotated_summary = annotated_summary_ch
    final_vep_annotated_vcf = compressed_vcf_ch
	versions = ch_versions
}
