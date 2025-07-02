process ARRIBA_VISUALIZATION {

    tag { "${meta.id}" }

    input:
    tuple val(meta), path(fusions_file), path(bam_file), path(bai_file)
    path gtf

    output:
    path("*.pdf"), emit: fusion_plot
    path("versions.yml"), emit: versions

    script:
    """
    set -euo pipefail

    echo "Running fusion visualization for sample: ${meta.id} (${meta.strandedness})"
    echo "Fusions file: \$(realpath ${fusions_file})"
    echo "GTF file: \$(realpath ${gtf})"

    if [ ! -s "${fusions_file}" ] || ! awk 'NR > 1 { exit 1 }' "${fusions_file}"; then
        echo "No fusions detected for ${meta.id}. Creating placeholder plot."
        touch "${meta.id}_${meta.strandedness}.fusion_plot.pdf"
    else
        echo "Fusions detected, generating plot..."
        draw_fusions.R \
            --fusions="${fusions_file}" \
            --alignments="${bam_file}" \
            --output="${meta.id}_${meta.strandedness}.fusion_plot.pdf" \
            --annotation="${gtf}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "ARRIBA_VISUALIZATION":
      arriba: "\$(arriba -h 2>&1 | grep 'Version:' | sed 's/Version: //')"
    END_VERSIONS
    """
}
