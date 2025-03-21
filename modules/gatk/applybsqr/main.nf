process GATK_APPLYBQSR {
    tag { sample_id }
    label 'process_high'

    container params.gatk_container
    publishDir params.recalibrated_bams_outdir, mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai), path(recal_table), path(interval)
    path genome_fasta
    path index
    path dict

    output:
    tuple val(sample_id), val(strandedness),
          path("${sample_id}_${interval.baseName}_recalibrated.bam"),
          path("${sample_id}_${interval.baseName}_recalibrated.bai"), emit: recalibrated_bam
    path("versions.yml"), emit: versions

    script:
    def interval_command = interval ? "--intervals ${interval}" : ""

    """
    THREADS=\${task.cpus}

    echo "Applying BQSR for sample: ${sample_id}, interval: ${interval.baseName}"

    # Step 1: ApplyBQSR
    gatk ApplyBQSR \\
        -R "${genome_fasta}" \\
        -I "${bam}" \\
        --bqsr-recal-file "${recal_table}" \\
        -O "${sample_id}_${interval.baseName}_recalibrated.bam" \\
        ${interval_command}

    # Verify output BAM exists
    if [ ! -s "${sample_id}_${interval.baseName}_recalibrated.bam" ]; then
        echo "Error: Recalibrated BAM not generated for ${sample_id}" >&2
        exit 1
    fi

    # Step 2: Index BAM
    gatk BuildBamIndex \\
        -I "${sample_id}_${interval.baseName}_recalibrated.bam"

    # Capture version
    gatk_version=\$(gatk --version | grep -Eo '[0-9.]+' | head -n 1)
    cat <<EOF > versions.yml
    "${task.process}":
      gatk: "\${gatk_version}"
    EOF
    """
}
