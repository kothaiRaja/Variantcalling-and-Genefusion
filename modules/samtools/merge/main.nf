process MERGE_BAMS {
    tag { "${meta.id}_${task.process}" }
    label 'process_medium'

    container params.samtools_container
    publishDir params.merge_bam_outdir, mode: "copy"

    input:
    tuple val(meta), path(bam_list), path(bai_list)

    output:
    tuple val(meta), 
          path("${meta.id}_merged.bam"), 
          path("${meta.id}_merged.bam.bai"), emit: merged_bams
    path("versions.yml"), emit: versions

    script:
    def sample_id = meta.id
    def bam_files = bam_list.join(' ')

    return """
    echo "Merging BAM files for sample: ${sample_id}"

    # Merge BAMs directly using list expansion
    samtools merge -@ ${task.cpus} -o "${sample_id}_merged.bam" ${bam_files}

    # Index the merged BAM
    samtools index "${sample_id}_merged.bam"

    # Capture Samtools version
    samtools_version=\$(samtools --version | head -n 1 | awk '{print \$2}')

cat <<EOF > versions.yml
"${task.process}":
  samtools: "\${samtools_version}"
EOF

    echo "Merge complete for sample: ${sample_id}"
    """
}
