process SAMTOOLS_IDXSTATS {
    tag { "${sample_id}_${task.process}" }
    label 'process_low'

    container params.samtools_container
    publishDir params.samtools_idx_outdir, mode: "copy", pattern: "*_idxstats.*"

    input:
    tuple val(sample_id), path(bam), path(bai), val(strandedness)

    output:
    tuple val(sample_id), path("${sample_id}_idxstats.txt"), val(strandedness), emit: idxstats
    path("versions.yml"), emit: versions

    script:
    """
    echo "Running samtools idxstats for sample: ${sample_id}"

    samtools idxstats "${bam}" > "${sample_id}_idxstats.txt"

    # Capture samtools version
    samtools_version=\$(samtools --version | head -n 1 | awk '{print \$2}')
cat <<EOF > versions.yml
"${task.process}":
  samtools: "\${samtools_version}"
EOF

    echo "Finished samtools idxstats
	"""
}