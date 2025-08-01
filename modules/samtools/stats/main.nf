process SAMTOOLS_STATS {
    tag { "${meta.id}_${task.process}" }
    label 'process_low'

    container params.samtools_container
    publishDir params.samtools_stats_outdir, mode: "copy", pattern: "*_stats_report.*"

    input:
    tuple val(meta), path(sorted_bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_stats_report.txt"), emit: stats
    path("versions.yml"), emit: versions

    script:
    def sample_id = meta.id

    """
    # Run samtools stats on the sorted BAM file
    samtools stats ${sorted_bam} > ${sample_id}_stats_report.txt

    # Capture samtools version
    samtools_version=\$(samtools --version | head -n 1 | awk '{print \$2}')
    cat <<EOF > versions.yml
"${task.process}":
  samtools: "\${samtools_version}"
EOF
    """
}
