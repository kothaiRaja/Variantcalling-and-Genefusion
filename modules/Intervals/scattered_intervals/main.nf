process SCATTER_INTERVAL_LIST {
    tag { "${task.process}" }
    label 'process_medium'

    container params.gatk_container
    publishDir params.scatter_intervals_outdir, mode: 'copy'

    input:
    tuple val(meta), path(interval_list)
    path genome_fasta 
    path genome_index
    path genome_dict

    output:
    tuple val(meta), path("*.interval_list"), emit: scattered_intervals
    path("versions.yml"), emit: versions

    script:
    """
    gatk SplitIntervals \
      -R ${genome_fasta} \
      -L ${interval_list} \
      --scatter-count ${params.scatter_count} \
      --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
      -O ./

    # Log version
    gatk_version=\$(gatk --version | head -n 1)

cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
