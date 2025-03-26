process ARRIBA_VISUALIZATION {
    tag { sample_id }

    label 'process_medium'

    container params.arriba_container
    publishDir params.visualisation_outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(fusions_tsv)
    path gtf

    output:
    tuple val(sample_id), path("*.pdf"), emit: fusion_plot
    path "versions.yml", emit: versions

    script:
    def prefix = sample_id
    def args = task.ext.args ?: ''
    def cytobands_arg = params.cytobands ? "--cytobands=${params.cytobands}" : ""
    def protein_domains_arg = params.protein_domains ? "--proteinDomains=${params.protein_domains}" : ""

    """
    set -euo pipefail

    echo " Running fusion visualization for sample: ${sample_id}"
    echo " Fusions file: \$(realpath ${fusions_tsv})"
    echo " GTF file: \$(realpath ${gtf})"

    if [ ! -s "${fusions_tsv}" ] || ! awk 'NR > 1 { exit 1 }' ${fusions_tsv}; then
        echo " No fusions detected for ${sample_id}. Creating placeholder plot."
        echo "No fusions found for ${sample_id}." | convert -background white -fill black -font Helvetica -pointsize 20 text:- ${prefix}.fusion_plot.pdf
    else
        echo " Fusions detected, generating plot..."
        draw_fusions.R \\
            --fusions=${fusions_tsv} \\
            --alignments=${bam} \\
            --output=${prefix}.fusion_plot.pdf \\
            --annotation=${gtf} \\
            ${cytobands_arg} \\
            ${protein_domains_arg} \\
            ${args}
    fi

    arriba_version=\$(arriba -h | grep 'Version:' 2>&1 | sed 's/Version:\\s//')

cat <<EOF > versions.yml
"${task.process}":
  arriba: "\${arriba_version}"
EOF
    """
}
