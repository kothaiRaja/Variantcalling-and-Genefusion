process ARRIBA_VISUALIZATION {
    tag { sample_id }

    cpus params.get('arriba_visualization_cpus', 4)
    memory params.get('arriba_visualization_memory', '16 GB')
    time params.get('arriba_visualization_time', '2h')

    container params.r_base_container
    publishDir params.visualisation_outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(fusions_tsv), val(strandedness)
    path r_script
    path fasta
    path gtf

    output:
    path "*.fusion_plot.pdf", emit: fusion_plot
    path "versions.yml", emit: versions

    script:
    """
    PREFIX=\$(basename ${fusions_tsv} .fusions.tsv)

    echo "Running fusion visualization for sample: ${sample_id}"
    echo "Resolved fusions file: \$(realpath ${fusions_tsv})"
    echo "Resolved GTF file: \$(realpath ${gtf})"

    if [ ! -f "\$(realpath ${fusions_tsv})" ]; then
        echo "Error: Fusions file not found." >&2
        exit 1
    fi
    if [ ! -f "\$(realpath ${gtf})" ]; then
        echo "Error: GTF file not found." >&2
        exit 1
    fi

    # Skip if no fusions
    if [ -z "\$(awk 'NR > 1 {print; exit}' ${fusions_tsv})" ]; then
        echo "No fusions found for ${sample_id}. Creating empty plot file."
        touch \${PREFIX}.fusion_plot.pdf
    else
        Rscript ${r_script} \\
            --fusions "\$(realpath ${fusions_tsv})" \\
            --annotation "\$(realpath ${gtf})" \\
            --output "\${PREFIX}.fusion_plot.pdf"
    fi

    # Version tracking
r_version=\$(Rscript --version 2>&1 | awk '{print \$NF}')

cat <<EOF > versions.yml
"${task.process}":
  R: "\${r_version}"
EOF
    """
}
