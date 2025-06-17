#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Custom RNA-seq Variant Calling & Gene Fusion Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INCLUDE WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUILD_REFERENCE }                 from './workflows/build_reference.nf'
include { RNA_VARIANT_CALLING_GENE_FUSION } from './workflows/varfuse.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONDITIONAL MAIN WORKFLOW WRAPPER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MASTER_PIPELINE {
    main:
    if (params.build_references) {
        BUILD_REFERENCE()
    } else {
        RNA_VARIANT_CALLING_GENE_FUSION()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FINAL Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    def base = System.getProperty('user.dir')

    // ======= If running build_references only =======
    if (params.build_references) {
        // === Reference-building mode ===
        if (params.ref_base == "${base}/reference") {
            println "️  You did not provide --ref_base."
            println "    Reference files will be saved to: ${base}/reference"
            println "    It's better to provide a folder with --ref_base"
        }

        println """
Reminder: When running --build_references, you must also provide:
A custom config file with tool database requirements and reference_file paths using:
-c custom.config
"""
        MASTER_PIPELINE()
    }

    // ======= If running full pipeline =======
    else {
        def missing = []
        if (params.ref_base == "${base}/reference")    missing << "--ref_base"
        if (params.cache_dir == "${base}/cache")       missing << "--cache_dir"
        if (params.resultsdir == "${base}/results")    missing << "--resultsdir"
        if (!params.samplesheet)                       missing << "--samplesheet"

        if (missing) {
            println """
	Note: The following parameters were not provided:
    ${missing.join(', ')}

    Default locations under '${base}' will be used.
    It's recommended to set these manually when running the full pipeline.
"""
        }
        MASTER_PIPELINE()
    }
}
