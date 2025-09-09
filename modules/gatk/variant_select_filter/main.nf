process GATK_VARIANT_SELECT_FILTER {
    tag { "${meta.id}_${task.process}" }

    label 'process_high'
    publishDir params.variant_filter_outdir, mode: 'copy'
    container params.gatk_container

    input:
    tuple val(meta), path(vcf_file), path(vcf_index)
    path genome
    path genome_index
    path genome_dict

    output:
    tuple val(meta), path("${meta.id}_filtered.vcf"), emit: filtered_vcf
    path("versions.yml"), emit: versions

    script:
    def sample_id = meta.id
    def avail_mem = 3
    if (task.memory) {
        avail_mem = task.memory.giga
    } else {
        log.info '[GATK VariantFiltration] No memory set — defaulting to 3GB.'
    }

    """
    echo "Selecting and filtering variants for sample: ${sample_id}"

    # 1. Select SNPs
    gatk --java-options "-Xmx${avail_mem}g" SelectVariants \\
        -R ${genome} \\
        -V ${vcf_file} \\
        --select-type-to-include SNP \\
        -O ${sample_id}_snps.vcf

    # 2. Select INDELs
    gatk --java-options "-Xmx${avail_mem}g" SelectVariants \\
        -R ${genome} \\
        -V ${vcf_file} \\
        --select-type-to-include INDEL \\
        -O ${sample_id}_indels.vcf

    # 3. Filter SNPs
    gatk --java-options "-Xmx${avail_mem}g" VariantFiltration \\
        -R ${genome} \\
        -V ${sample_id}_snps.vcf \\
        -O ${sample_id}_snps_filtered.vcf \\
        --filter-name "FAIL_SNP" \\
        --filter-expression "QD < 2.0 || FS > 30.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

    # 4. Filter INDELs
    gatk --java-options "-Xmx${avail_mem}g" VariantFiltration \\
        -R ${genome} \\
        -V ${sample_id}_indels.vcf \\
        -O ${sample_id}_indels_filtered.vcf \\
        --filter-name "FAIL_INDEL" \\
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

    # 5. Merge both filtered VCFs
    gatk --java-options "-Xmx${avail_mem}g" MergeVcfs \\
        -I ${sample_id}_snps_filtered.vcf \\
        -I ${sample_id}_indels_filtered.vcf \\
        -O ${sample_id}_filtered.vcf

    # 6. Capture GATK version
    gatk_version=\$(gatk --version | head -n 1)
    cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
