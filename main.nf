nextflow.enable.dsl = 2

// Import subworkflows
include { PREPROCESSING } from './subworkflows/preprocessing.nf'
include { STAR_ALIGN } from './subworkflows/star_align.nf'

workflow {

    // Ensure only one of the specialized workflows is selected
    if (params.only_qc && params.skip_star) { 
        error "You cannot enable both only_qc and skip_star at the same time. Set one to false."
    }

    //============================== PREPROCESSING / TRIMMING ===============================//

    trimmed_reads_ch = Channel.empty() // Initialize empty channel

    if (params.input_reads) {
        log.info "Using provided input reads for STAR Alignment..."
        trimmed_reads_ch = Channel.fromPath(params.input_reads)
    } else {
        log.info "No input reads provided. Running Preprocessing..."
        PREPROCESSING(params.samplesheet)
        
        trimmed_reads_ch  = PREPROCESSING.out.trimmed_reads
        fastp_reports_ch  = PREPROCESSING.out.fastp_reports
        qc_results_ch     = PREPROCESSING.out.qc_reports
        multiqc_quality   = PREPROCESSING.out.multiqc
    }

    // If user wants QC-only mode, exit after preprocessing
    if (params.only_qc) {
        log.info("QC completed. Exiting pipeline...")
        return
    }

    //==================================== STAR ALIGNMENT ====================================//

    if (!params.skip_star) {
        log.info "Running STAR Alignment..."
        STAR_ALIGN(trimmed_reads_ch, params.star_genome_index, params.gtf_annotation)
    } else {
        log.info "Skipping STAR Alignment as per user request."
    }
}
   
	
	
    
