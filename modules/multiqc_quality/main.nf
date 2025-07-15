process MultiQC {

    tag { "MultiQC" }
    label 'process_low'

    container params.multiqc_container
    publishDir params.multiqc_outdir, mode: 'copy'

    input:
    path(report_files)

    output:
    path "multiqc_report.html", emit: multiqc_report
    path "multiqc_data",        optional: true, emit: data
    path "multiqc_plots",       optional: true, emit: plots
    path "versions.yml",        emit: versions

    script:
    """
    echo "Running MultiQC on directly passed files..."

    multiqc ${report_files.join(' ')} -c multiqc_config.yml -o .

    # Capture MultiQC version
    multiqc_version=\$(multiqc --version | grep -oP '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo "not_detected")
    cat <<EOF > versions.yml
"MultiQC":
  multiqc: "\$multiqc_version"
EOF
    """
}
