process MultiQC {
  tag "MultiQC"
  label 'process_low'

  publishDir params.multiqc_quality_outdir, mode: "copy"
  container params.multiqc_quality_container

  input:
    path report_dir

  output:
    path "multiqc_report.html", emit: report
    path "versions.yml", emit: versions

  script:
  """
  echo "Running MultiQC in ${report_dir}"
  ls -lh ${report_dir}

  multiqc ${report_dir} -o .

  # Capture MultiQC version
multiqc_version=\$(multiqc --version | grep -oP '[0-9]+\\.[0-9]+(\\.[0-9]+)?')
cat <<EOF > versions.yml
 "${task.process}":
multiqc: "\${multiqc_version}"
EOF
  """
}
