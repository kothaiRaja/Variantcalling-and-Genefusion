process SPLIT_NCIGAR_READS {
    tag { "${sample_id}_${interval.baseName}" } 
    label 'process_high'

    container params.gatk_container
    publishDir params.split_ncigar_outdir, mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai), path(interval)
    path genome_fasta
    path index
    path genome_dict

    output:
    tuple val(sample_id), val(strandedness), 
          path("${sample_id}_split_${interval.baseName}.bam"), 
          path("${sample_id}_split_${interval.baseName}.bai"), emit: split_interval_bams
    path("versions.yml"), emit: versions

    script:
    def avail_mem = task.memory ? task.memory.giga : 3
    def interval_command = interval ? "--intervals ${interval}" : ""

    """
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
    gatk_version=\$(gatk --version | grep -Eo '[0-9.]+' | head -n 1)
    cat <<EOF > versions.yml
    "${task.process}":
      gatk: "\${gatk_version}"
    EOF

    echo "SplitNCigarReads completed for sample: ${sample_id} on interval: ${interval.baseName}"
    """
}
