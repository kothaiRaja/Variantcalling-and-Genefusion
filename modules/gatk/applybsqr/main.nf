process GATK_APPLYBQSR {
    tag { "${meta.id}_${task.process}" }

    label 'process_high'

    container params.gatk_container
    publishDir params.recalibrated_bams_outdir, mode: "copy"

    input:
    tuple val(meta), path(bam), path(bai), path(recal_table), path(intervals)
    path genome_fasta
    path index
    path dict

    output:
    tuple val(meta),
          path("${meta.id}_recalibrated_${intervals.getBaseName()}.bam"),
          path("${meta.id}_recalibrated_${intervals.getBaseName()}.bai"),
          path(intervals), emit: recalibrated_bam
    path("versions.yml"), emit: versions

    script:
    def avail_mem = task.memory ? task.memory.giga : 3
    def interval_option = intervals ? "-L ${intervals}" : ""

    """
    THREADS=${task.cpus}

    echo "Applying BQSR for sample: ${meta.id}"

    # Step 1: ApplyBQSR
    gatk --java-options "-Xmx${avail_mem}g" ApplyBQSR \\
        -R "${genome_fasta}" \\
        -I "${bam}" \\
        --bqsr-recal-file "${recal_table}" \\
        ${interval_option} \\
        -O "${meta.id}_recalibrated_${intervals.getBaseName()}.bam"

    # Verify output BAM exists
    if [ ! -s "${meta.id}_recalibrated_${intervals.getBaseName()}.bam" ]; then
        echo "Error: Recalibrated BAM not generated for ${meta.id}" >&2
        exit 1
    fi

    # Step 2: Index BAM
    gatk BuildBamIndex \\
         -I "${meta.id}_recalibrated_${intervals.getBaseName()}.bam"

    # Capture version
    gatk_version=\$(gatk --version | head -n 1)

    cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
