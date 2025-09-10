process ARRIBA_VISUALIZATION {

  tag { "${meta.id}" }
  label 'process_medium'
  container params.arriba_container
  publishDir params.arriba_outdir, mode: 'copy'

  input:
  tuple val(meta), path(fusions_file), path(bam_file), path(bai_file)
  path gtf

  output:
  path("*.pdf"), emit: fusion_plot
  path("versions.yml"), emit: versions

  script:
  """
  set -euo pipefail

  echo "Running fusion visualization for sample: ${meta.id} (${meta.strandedness})"
  echo "Fusions file: \$(realpath ${fusions_file})"
  echo "GTF file: \$(realpath ${gtf})"

  # Is there at least one data row (beyond header)?
  has_data=\$(awk 'NR>1{print 1; exit} END{print 0}' "${fusions_file}")

  if [ ! -s "${fusions_file}" ] || [ "\$has_data" -eq 0 ]; then
    echo "No fusions detected for ${meta.id}. Creating empty PDF."
    printf "%s\n" "%PDF-1.1" "%EOF" > "${meta.id}_${meta.strandedness}.fusion_plot.pdf" || touch "${meta.id}_${meta.strandedness}.fusion_plot.pdf"
  else
    echo "Fusions detected, generating plot..."
    draw_fusions.R \
      --fusions="${fusions_file}" \
      --alignments="${bam_file}" \
      --output="${meta.id}_${meta.strandedness}.fusion_plot.pdf" \
      --annotation="${gtf}" \
    || { echo "draw_fusions.R failed — creating empty PDF"; printf "%s\n" "%PDF-1.1" "%EOF" > "${meta.id}_${meta.strandedness}.fusion_plot.pdf" || touch "${meta.id}_${meta.strandedness}.fusion_plot.pdf"; }
  fi

  cat <<-END_VERSIONS > versions.yml
"ARRIBA_VISUALIZATION":
  arriba: "\$(arriba -h 2>&1 | grep 'Version:' | sed 's/Version: //')"
END_VERSIONS
  """
}
