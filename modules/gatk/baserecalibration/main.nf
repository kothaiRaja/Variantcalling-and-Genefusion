process GATK_BASERECALIBRATOR {
    tag { sample_id }
    label 'process_high'

    container params.gatk_container
    publishDir params.recalibration_table_outdir, mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai), path(interval)
    path genome_fasta
    path index
    path dict
    path known_variants
    path known_variants_index

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_recal_data.table"), emit: recal_table
    path("versions.yml"), emit: versions

    script:
    def interval_command = interval ? "--intervals ${interval}" : ""

    """
    THREADS=${task.cpus}

    echo "Running GATK BaseRecalibrator for sample: ${sample_id}"

    gatk BaseRecalibrator \\
        -R "${genome_fasta}" \\
        -I "${bam}" \\
        --known-sites "${known_variants}" \\
        -O "${sample_id}_recal_data.table" \\
        ${interval_command}

    # Check if recalibration table was created
    if [ ! -s "${sample_id}_recal_data.table" ]; then
        echo "Error: Recalibration table not generated for ${sample_id}" >&2
        exit 1
    fi

    # Capture GATK version
    gatk_version=\$(gatk --version | grep -Eo '[0-9.]+' | head -n 1)
    cat <<EOF > versions.yml
    "${task.process}":
      gatk: "\${gatk_version}"
    EOF

    echo "BaseRecalibrator completed for ${sample_id}"
    """
}
