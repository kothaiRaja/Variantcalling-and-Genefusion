process GATK_BASERECALIBRATOR {
    tag { "${sample_id}_${task.process}" }

    label 'process_high'

    container params.gatk_container
    publishDir params.recalibration_table_outdir, mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai)
    path genome_fasta
    path index
    path dict
    path known_variants
    path known_variants_index

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_recal_data.table"), emit: recal_table
    path("versions.yml"), emit: versions

    script:
    def avail_mem = 3
    if (task.memory) {
        avail_mem = task.memory.giga
    } else {
        log.info '[GATK BaseRecalibrator] No memory set — defaulting to 3GB.'
    }

    """
    THREADS=${task.cpus}

    echo "Running GATK BaseRecalibrator for sample: ${sample_id}"

    gatk --java-options "-Xmx${avail_mem}g" BaseRecalibrator \\
        -R "${genome_fasta}" \\
        -I "${bam}" \\
        --known-sites "${known_variants}" \\
        -O "${sample_id}_recal_data.table"

    # Check if recalibration table was created
    if [ ! -s "${sample_id}_recal_data.table" ]; then
        echo "Error: Recalibration table not generated for ${sample_id}" >&2
        exit 1
    fi

    # Capture GATK version
    gatk_version=\$(gatk --version | head -n 1)

    cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF

    echo "BaseRecalibrator completed for ${sample_id}"
    """
}
