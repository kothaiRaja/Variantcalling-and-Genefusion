process SCATTER_INTERVAL_LIST {
    tag { "${task.process}" }
	label 'process_medium'

    container params.gatk_container
    publishDir params.scatter_intervals_outdir, mode: "copy"

    input:
    tuple val(meta), path(interval_list)
	path(genome_dict)

    output:
    tuple val(meta), path("*.interval_list"), emit: scattered_intervals
	path ("versions.yml"), emit: versions

    script:
    """
    mkdir -p scattered_intervals
    gatk IntervalListTools \
        --INPUT ${interval_list} \
        --OUTPUT scattered_intervals \
        --SCATTER_COUNT ${params.scatter_count} \
        --UNIQUE true

    # Move and rename the output files to the working directory with unique names
    for f in scattered_intervals/*/*; do
        dir_name=\$(basename \$(dirname "\$f"))  # Get subdirectory name
        file_name=\$(basename "\$f")            # Get original file name
        mv "\$f" "scattered_intervals/\${dir_name}_\${file_name}.interval_list"
    done

    # Move the uniquely renamed files to the working directory
    mv scattered_intervals/*.interval_list .
    rm -r scattered_intervals
	
	# Log the generated intervals
    echo "Scattered intervals:"
    cat *.interval_list
	
	# Capture GATK version
    gatk_version=\$(gatk --version | head -n 1)

cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
