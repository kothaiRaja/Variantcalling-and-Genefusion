process MULTIQC_REPORT {
    tag "MultiQC"
	container "https://depot.galaxyproject.org/singularity/multiqc%3A1.26--pyhdfd78af_0"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_report_files 

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc ${qc_report_files.join(' ')} -o .
    """
}