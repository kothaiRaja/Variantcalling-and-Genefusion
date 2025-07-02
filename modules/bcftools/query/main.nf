process BCFTOOLS_QUERY {
    tag { "${meta}_${task.process}" }
    label 'process_low'

    container params.bcftools_container
    publishDir params.bcftools_query_outdir, mode: "copy", pattern: "*_variant_summary.txt"

    input:
    tuple val(meta), path(filtered_vcf), path(tbi)

    output:
    tuple val(meta), path("${meta}_variant_summary.txt"), emit: query_output
    path("versions.yml"), emit: versions

    script:
    """
    echo "Querying variant fields for sample: ${meta}"

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/DP\\t%INFO/AF\\n' "${filtered_vcf}" > "${meta}_variant_summary.txt"

    # Capture bcftools version
    bcftools_version=\$(bcftools --version | head -n 1 | awk '{print \$2}')
cat <<EOF > versions.yml
"${task.process}":
  bcftools: "\${bcftools_version}"
EOF
    """
}
