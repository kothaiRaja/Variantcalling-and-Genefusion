nextflow.enable.dsl = 2

// Include required processes
include { CONCAT_FASTQ } from '../modules/cat_fastq/main.nf'
include { FASTQC_RAW } from '../modules/fastqc/main.nf'
include { TRIM_READS } from '../modules/fastp/main.nf'
include { MultiQC_quality } from '../modules/multiqc_quality/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nfcore/software_versions/main.nf'

workflow PREPROCESSING {
    take:
    samplesheet
	dump_script

    main:

    log.info "Starting Preprocessing Steps..."
	
	ch_versions = Channel.empty()
	qc_results_ch = Channel.empty()
	trimmed_reads_ch = Channel.empty()
	reports_ch = Channel.empty()
	combined_channel = Channel.empty()



    // Step 1: Read samplesheet
    samples_ch = Channel.fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)], row.strandedness) }

    // Step 2: Concatenate FASTQ files if required
    concatenated_reads_ch = params.concatenate ? 
        CONCAT_FASTQ(samples_ch) : 
        samples_ch.map { sample_id, reads, strandedness -> tuple(sample_id, reads[0], reads[1], strandedness) }

    // Step 3: Perform FastQC on raw reads
    qc_results = FASTQC_RAW(concatenated_reads_ch)
	
	qc_results_ch = qc_results_ch.mix(FASTQC_RAW.out.qc_results)
	reports_ch = reports_ch.mix(FASTQC_RAW.out.qc_results.map { it[1] })
	ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())


    // Step 4: Trim reads using Fastp
    trimmed_reads = TRIM_READS(concatenated_reads_ch)
	
	trimmed_reads_ch =trimmed_reads_ch.mix(TRIM_READS.out.trimmed_reads)  
	reports_ch = reports_ch.mix(TRIM_READS.out.fastp_reports.map { it[1] })
	reports_ch.view { " QC report: $it" }
	
	ch_versions = ch_versions.mix(TRIM_READS.out.versions.first())
	
	// Step 5: Version Dump
	

	CUSTOM_DUMPSOFTWAREVERSIONS(
		ch_versions.unique().collectFile(name: "software_versions_input.yml"),
		dump_script
	)
	
	reports_ch = reports_ch.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)

    // Step 6: Collect QC reports for MultiQC
    
    
    // FLATTEN the tuples into a flat list of files
	collected_reports_ch = reports_ch.collect()

	// Optional: Debug what goes into MultiQC
	collected_reports_ch.view { " Report passed to MultiQC: $it" }

	// Pass to MultiQC
	multiqc_quality = MultiQC_quality(collected_reports_ch)
	ch_versions = ch_versions.mix(MultiQC_quality.out.versions.first())

    log.info " Preprocessing Completed."

    emit:
    trimmed_reads = trimmed_reads_ch
    reports        = reports_ch
    multiqc        = multiqc_quality.report
    versions       = ch_versions

}
