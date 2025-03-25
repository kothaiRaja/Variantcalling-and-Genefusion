process MultiQC {
  tag "MultiQC_quality"
  label 'process_low'

  publishDir params.multiqc_qualtiy_outdir, mode: "copy"
  container params.multiqc_quality_container

  input:
    path report_files

  output:
    path "multiqc_report.html", emit: report
    path "versions.yml", emit: versions

  script:
  """
  echo "Running MultiQC on quality control reports..."

  multiqc ${report_files.join(' ')} -o .

  # Capture MultiQC version
  multiqc_version=\$(multiqc --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+(\\.[0-9]+)?')
  cat <<EOF > versions.yml
  "${task.process}":
    multiqc: "\${multiqc_version}"
  EOF
  """
}

