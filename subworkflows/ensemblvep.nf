nextflow.enable.dsl=2

// Include modules
include { ANNOTATEVARIANTS_VEP } from '../modules/ensemblvep/main.nf'
include { BGZIP_TABIX_ANNOTATIONS } from '../modules/tabix/bziptabix/main.nf'

workflow VEP_ANNOTATION_WORKFLOW {

    take:
    vcf_files_ch 
	vep_cache
	clinvar_vcf
	clinvar_vcf_tbi
	genome_assembly
	vep_version
	species
	
    main:
    log.info "Starting Ensembl VEP annotation..."

    // Step 1: Annotate VCF using Ensembl VEP
    annotated_vcf_ch = ANNOTATEVARIANTS_VEP(
        vcf_files_ch,
		vep_cache,
        clinvar_vcf,
        clinvar_vcf_tbi,
        genome_assembly,
		vep_version,
		species
       
    )

    log.info "Compressing and indexing annotated VCFs..."

    // Step 2: Compress and index annotated VCFs
    compressed_vcf_ch = BGZIP_TABIX_ANNOTATIONS(annotated_vcf_ch[0])

    emit:
	annotated_variants = annotated_vcf_ch.annotated_vcf
	annotated_html = annotated_vcf_ch.summary_html
    final_vep_annotated_vcf = compressed_vcf_ch.compressed_indexed
}
