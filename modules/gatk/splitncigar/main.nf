process SPLIT_NCIGAR_READS {
    tag { "${meta.id}_${interval.baseName}_${task.process}" }

    label 'process_high'

    container params.gatk_container
    publishDir params.split_ncigar_outdir, mode: "copy"

    input:
    tuple val(meta), path(bam), path(bai), path(interval)
    path genome_fasta
    path index
    path genome_dict

    output:
    tuple val(meta), 
          path("${meta.id}_split_${interval.baseName}.bam"), 
          path("${meta.id}_split_${interval.baseName}.bai"), emit: split_interval_bams
    path("versions.yml"), emit: versions

    script:
    def sample_id = meta.id
    def avail_mem = task.memory ? task.memory.giga : 3
    def interval_command = interval ? "--intervals ${interval}" : ""

    return """
    echo "Running SplitNCigarReads for sample: ${sample_id} on interval: ${interval.baseName}"

    gatk --java-options "-Xmx${avail_mem}g" SplitNCigarReads \\
        -R "${genome_fasta}" \\
        -I "${bam}" \\
        -O "${sample_id}_split_${interval.baseName}.bam" \\
        --skip-mapping-quality-transform false \\
        --max-mismatches-in-overhang 1 \\
        --max-bases-in-overhang 50 \\
        --create-output-bam-index true \\
        --process-secondary-alignments true \\
        ${interval_command}

    # Capture version
    gatk_version=\$(gatk --version | head -n 1)

cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF

   
    """
}
