process BCFTOOLS_STATS {

    tag { "${meta.id}_${task.process}" }
    label 'process_low'

    container params.bcftools_container
    publishDir params.bcftools_stats_outdir, mode: "copy", pattern: "*_bcftools_stats.*"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}_bcftools_stats.txt"), emit: stats_output
    path("versions.yml"), emit: versions

    script:
    """
    echo "Running bcftools stats for sample: ${meta.id}"

    bcftools stats "${vcf}" > "${meta.id}_bcftools_stats.txt"

    bcftools_version=\$(bcftools --version | head -n 1 | awk '{print \$2}')
    cat <<EOF > versions.yml
"${task.process}":
  bcftools: "\${bcftools_version}"
EOF
    """
}
