process GATK_BASERECALIBRATOR {
    tag { "${meta.id}_${task.process}" }

    label 'process_high'

    container params.gatk_container
    publishDir params.recalibration_table_outdir, mode: "copy"

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)
    path genome_fasta
    path index
    path dict
    path known_variants
    path known_variants_index

    output:
    tuple val(meta), path("${meta.id}_recal_data.table"), emit: recal_table
    path("versions.yml"), emit: versions

    script:
    def avail_mem = task.memory ? task.memory.giga : 3
    def interval_option = intervals ? "-L ${intervals}" : ""
    def known_sites_command = known_variants.collect { "--known-sites ${it}" }.join(' ')

    """
    THREADS=${task.cpus}

    echo "Running GATK BaseRecalibrator for sample: ${meta.id}"

    gatk --java-options "-Xmx${avail_mem}g" BaseRecalibrator \\
        -R "${genome_fasta}" \\
        -I "${bam}" \\
        ${known_sites_command} \\
        ${interval_option} \\
        -O "${meta.id}_recal_data.table"

    if [ ! -s "${meta.id}_recal_data.table" ]; then
        echo "Error: Recalibration table not generated for ${meta.id}" >&2
        exit 1
    fi

    gatk_version=\$(gatk --version | head -n 1)

    cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF

    echo "BaseRecalibrator completed for ${meta.id}"
    """
}
