process COLLECT_VARIANT_CALLING_METRICS {

  tag { "${meta.id}_${task.process}" }
  label 'process_medium'
  
  container params.gatk_container
  publishDir params.variant_metrics_outdir, mode: 'copy'

  input:
  tuple val(meta), path(vcf), path(vcf_tbi)            
  path dbsnp_vcf
  path dbsnp_tbi                      
  path ref_fasta
  path ref_dict                       
                      

  output:
  tuple val(meta),
        path("${meta.id}.picard_vc.variant_calling_summary_metrics"),
        path("${meta.id}.picard_vc.variant_calling_detail_metrics"),
        emit: metrics
  path "versions.yml", emit: versions

  script:
  def mem = task.memory?.giga ?: 4

  """
  set -euo pipefail

  gatk CollectVariantCallingMetrics \\
    -I "${vcf}" \\
    -O "${meta.id}.picard_vc" \\
    -DBSNP "${dbsnp_vcf}" \\
    -R "${ref_fasta}" \\
    -SD "${ref_dict}" \\
    --VALIDATION_STRINGENCY LENIENT

  
  # Versions
  gatk_version=\$(gatk --version | head -n 1)
  cat > versions.yml <<EOF
"${task.process}":
  gatk: "\${gatk_version}"
EOF
  """
}
