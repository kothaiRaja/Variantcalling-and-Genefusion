process MULTIQC_REPORT {
    tag "Generate MultiQC report"

    container "https://depot.galaxyproject.org/singularity/multiqc%3A1.20--pyhdfd78af_0"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files_dir

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc ${qc_files_dir} -o .
    """
}
