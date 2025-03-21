process SCATTER_INTERVAL_LIST {
    tag "Scatter interval list"
    label 'process_medium'

    container params.gatk_container
    publishDir params.scatter_intervals_outdir, mode: "copy"

    input:
    tuple val(meta), path(interval_list)
    path(genome_dict)

    output:
    tuple val(meta), path("*.interval_list"), emit: scattered_intervals
    path("versions.yml"), emit: versions

    script:
    """
    mkdir -p scattered_intervals

    gatk IntervalListTools \\
        --INPUT "${interval_list}" \\
        --OUTPUT scattered_intervals \\
        --SCATTER_COUNT ${params.scatter_count} \\
        --UNIQUE true \\
        --SEQUENCE_DICTIONARY "${genome_dict}"

    # Rename scatter outputs for uniqueness
    for f in scattered_intervals/*/*; do
        dir_name=\$(basename \$(dirname "\$f"))
        file_name=\$(basename "\$f")
        mv "\$f" "scattered_intervals/\${dir_name}_\${file_name}.interval_list"
    done

    # Move final scattered interval files to working directory
    mv scattered_intervals/*.interval_list .
    rm -r scattered_intervals

    # Capture GATK version
    gatk_version=\$(gatk --version | grep -Eo '[0-9.]+' | head -n 1)
    cat <<EOF > versions.yml
    "${task.process}":
      gatk: "\${gatk_version}"
    EOF
    """
}
