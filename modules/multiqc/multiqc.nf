process multiqc {
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
	container "https://depot.galaxyproject.org/singularity/multiqc%3A1.24.1--pyhdfd78af_0"
    input:
    path results_dir

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc ${results_dir} --outdir .
    """
}
