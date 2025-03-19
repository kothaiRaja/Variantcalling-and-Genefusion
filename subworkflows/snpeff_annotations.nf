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
    log.info " Starting Variant Annotation Workflow..."

   annotated_variants_ch = ANNOTATE_VARIANTS(filtered_variants_ch, snpEffJar, snpEffConfig, snpEffDbDir, genomedb)
   
   compressed_variants_ch = BGZIP_TABIX_ANNOTATIONS(annotated_variants_ch[0])

   
    // Emit Final Annotated Variants
    emit:
    annotated_variants = annotated_variants_ch.annotated_vcf
	annotated_html = annotated_variants_ch.summary_html
	annotated_csv = annotated_variants_ch.annotation_csv
	final_annotated_variants = compressed_variants_ch.compressed_indexed
}
