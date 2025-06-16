#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ========================== Custom Config Check ========================== //
def customConfigPath = file("custom.config")

if (customConfigPath.exists()) {
    println "[INFO]  'custom.config' found and may be in use (if passed via -c)."
} else {
    println "[WARN]   'custom.config' not found in the current directory."

    if (!params.reference_genome) {
        println """
[ERROR]  Missing required parameter: 'reference_genome'.

It seems 'custom.config' was not passed or does not contain this required value.

 To fix this, either:
   • Create a config file:    cp custom.config.template custom.config
   • Or run the pipeline with:
     nextflow run main.nf -c nextflow_main.config -c custom.config -profile singularity
        """.stripIndent()
        exit 1
    }
}




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Custom RNA-seq Variant Calling & Gene Fusion Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INCLUDE WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUILD_REFERENCE }             from './workflows/build_reference_main.nf'
include { BUILD_REFERENCE_TEST }         from './workflows/build_reference_test.nf'
include { RNA_VARIANT_CALLING_GENE_FUSION } from './workflows/varfuse.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONDITIONAL MAIN WORKFLOW CONTROLLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MASTER_PIPELINE {

    main:

    if (params.build_references) {
        BUILD_REFERENCE()
    } else if (params.build_references_test) {
        BUILD_REFERENCE_TEST()
    } else {
        RNA_VARIANT_CALLING_GENE_FUSION()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FINAL ENTRYPOINT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    MASTER_PIPELINE()
}
