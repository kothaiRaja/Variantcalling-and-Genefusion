nextflow.enable.dsl = 2

// Include required processes
include { CONCAT_FASTQ } from '../modules/cat_fastq/main.nf'
include { FASTQC_RAW   } from '../modules/fastqc/main.nf'
include { TRIM_READS   } from '../modules/fastp/main.nf'

workflow PREPROCESSING {

    take:
    validated_reads

    main:
//    log.info " Starting Preprocessing Steps..."

    ch_versions       = Channel.empty()
    qc_results_ch     = Channel.empty()
    trimmed_reads_ch  = Channel.empty()
    reports_ch        = Channel.empty()

    // Step 1: Concatenate FASTQs (optional)
    concatenated_reads_ch = params.concatenate ?
        CONCAT_FASTQ(validated_reads) :
        validated_reads.map { meta, reads ->
            tuple(meta, reads[0], reads[1])
        }

    // Step 2: Run FastQC on raw reads
	qc_results = FASTQC_RAW(concatenated_reads_ch)
	qc_results_ch = qc_results_ch.mix(FASTQC_RAW.out.qc_results)
	reports_ch = reports_ch
		.mix(FASTQC_RAW.out.qc_results.map { meta, r1_zip, r1_html, r2_zip, r2_html -> r1_zip })
		.mix(FASTQC_RAW.out.qc_results.map { meta, r1_zip, r1_html, r2_zip, r2_html -> r1_html })
		.mix(FASTQC_RAW.out.qc_results.map { meta, r1_zip, r1_html, r2_zip, r2_html -> r2_zip })
		.mix(FASTQC_RAW.out.qc_results.map { meta, r1_zip, r1_html, r2_zip, r2_html -> r2_html })

	// Version tracking remains unchanged
	ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())


    // Step 3: Trim reads using Fastp
    trimmed_reads = TRIM_READS(concatenated_reads_ch)

    trimmed_reads_ch  = trimmed_reads_ch.mix(TRIM_READS.out.trimmed_reads)
    fastp_reports_ch  = TRIM_READS.out.fastp_reports
    reports_ch = reports_ch.mix(TRIM_READS.out.fastp_reports.map { it[1] })


    ch_versions       = ch_versions.mix(TRIM_READS.out.versions.first())

//    trimmed_reads_ch.view { " TRIMMED READ: $it" }
//	reports_ch.view { "  REPORT FOR MULTIQC: $it" }


//    log.info " Preprocessing Completed."

    emit:
        qc_results     = qc_results_ch
        fastp_reports  = fastp_reports_ch
        trimmed_reads  = trimmed_reads_ch
        reports        = reports_ch
        versions       = ch_versions
}
