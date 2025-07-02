process SAMTOOLS_SORT_INDEX {
    tag { "${meta.id}_${task.process}" }
    label 'process_medium'

    container params.samtools_container
    publishDir params.samtools_outdir, mode: "copy"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_sorted.bam"), path("${meta.id}_sorted.bam.bai"), emit: bam_sorted
    path("versions.yml"), emit: versions

    script:
    def sample_id = meta.id

    """
    echo "Sorting and indexing BAM for sample: ${sample_id}"

    # Sort the BAM file
    samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam "${bam}"

    # Index the sorted BAM file
    samtools index ${sample_id}_sorted.bam

    # Capture the version
    samtools_version=\$(samtools --version | head -n 1 | awk '{print \$2}')
    cat <<EOF > versions.yml
"${task.process}":
  SAMTOOLS: "\${samtools_version}"
EOF
    """
}
