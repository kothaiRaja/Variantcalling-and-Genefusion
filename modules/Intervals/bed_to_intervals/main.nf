process BED_TO_INTERVAL_LIST {
    tag { "${meta_id}_${bed_file.baseName}_${task.process}" }

    label 'process_medium'

    container params.gatk_container
    publishDir params.bed_to_interval_outdir, mode: "copy"

    input:
    tuple val(meta_id), path(bed_file)
    path(genome_fasta)
    path(genome_dict)

    output:
    tuple val(meta_id), path("${bed_file.baseName}.interval_list"), emit: interval_list
	path("versions.yml"),emit:versions

    script:
    """
    echo "Converting BED to interval list for: ${meta_id}"

    gatk BedToIntervalList \\
        -I "${bed_file}" \\
        -O "${bed_file.baseName}.interval_list" \\
        -SD "${genome_dict}"

    # Capture GATK version
    gatk_version=\$(gatk --version | head -n 1)

cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
