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
    echo " Running MultiQC on collected files..."
    mkdir multiqc_inputs

    for file in ${report_files}; do
        cp \$file multiqc_inputs/
    done

    ls -lh multiqc_inputs

    multiqc multiqc_inputs -c multiqc_config.yml -o .

    # Capture MultiQC version
    multiqc_version=\$(multiqc --version | grep -oP '[0-9]+\\.[0-9]+(\\.[0-9]+)?')
    cat <<EOF > versions.yml
"MultiQC":
  multiqc: "\$multiqc_version"
EOF
    """
}
