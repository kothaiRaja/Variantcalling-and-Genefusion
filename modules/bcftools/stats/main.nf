process BCFTOOLS_STATS {
    tag { sample_id }
    label 'process_low'

    container params.bcftools_container
    publishDir params.bcftools_stats_outdir, mode: "copy", pattern: "*_bcftools_stats.*"

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}_bcftools_stats.txt"), emit: stats_output
    path("versions.yml"), emit: versions

    script:
    """
    echo "Running bcftools stats for sample: ${sample_id}"

    bcftools stats "${vcf}" > "${sample_id}_bcftools_stats.txt"

    # Capture bcftools version
    bcftools_version=\$(bcftools --version | head -n 1 | awk '{print \$2}')
    cat <<EOF > versions.yml
    "${task.process}":
      bcftools: "\${bcftools_version}"
    EOF
    """
}
