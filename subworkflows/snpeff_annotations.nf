nextflow.enable.dsl = 2

// Import annotation modules
include { ANNOTATE_VARIANTS } from '../modules/snpeff/main.nf'
include { BGZIP_TABIX_ANNOTATIONS } from '../modules/tabix/bziptabix/main.nf'



workflow VARIANT_ANNOTATION {

    take:
    filtered_variants_ch  // From VARIANT_CALLING
    snpEffJar
    snpEffConfig
    snpEffDbDir
    genomedb

    main:
	ch_versions = Channel.empty()
	
//    log.info " Starting Variant Annotation Workflow..."

   annotated_variants = ANNOTATE_VARIANTS(filtered_variants_ch, snpEffJar, snpEffConfig, snpEffDbDir, genomedb)
   annotated_variants_ch = ANNOTATE_VARIANTS.out.annotated_vcf
   annotated_summary_ch = ANNOTATE_VARIANTS.out.summary
   annotated_html_ch = ANNOTATE_VARIANTS.out.summary_html
   ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions) 
   
   compressed_variants = BGZIP_TABIX_ANNOTATIONS(annotated_variants_ch)
   compressed_variants_ch = BGZIP_TABIX_ANNOTATIONS.out.compressed_indexed
   ch_versions = ch_versions.mix(BGZIP_TABIX_ANNOTATIONS.out.versions)
   

   
    // Emit Final Annotated Variants
    emit:
    annotated_variants = annotated_variants_ch
	annotated_summary = annotated_summary_ch
	annotated_html = annotated_html_ch
	final_annotated_variants = compressed_variants_ch
	versions = ch_versions
}
