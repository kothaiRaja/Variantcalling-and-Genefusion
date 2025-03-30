process GATK_VARIANT_FILTER {
    tag { "${sample_id}_${task.process}" }

    label 'process_high'

    container params.gatk_container
    publishDir params.variant_filter_outdir, mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path genome
    path genome_index
    path genome_dict

    output:
    tuple val(sample_id), 
          path("${sample_id}_filtered.vcf.gz"), 
          path("${sample_id}_filtered.vcf.gz.tbi"), emit: filtered_vcf
    path("versions.yml"), emit: versions

    script:
	
	def avail_mem = 3
if (task.memory) {
    avail_mem = task.memory.giga
} else {
    log.info '[GATK VariantFiltration] No memory set — defaulting to 3GB.'
}

	
    """
    THREADS=${task.cpus}

    echo "Running GATK VariantFiltration for sample: ${sample_id}"

    gatk --java-options "-Xmx${avail_mem}g" VariantFiltration \\
        -R "${genome}" \\
        -V "${vcf_file}" \\
        --cluster-window-size ${params.gatk_vf_window_size} \\
        --cluster-size ${params.gatk_vf_cluster_size} \\
        --filter-name "LowQual"        --filter-expression "QUAL < ${params.gatk_vf_qual_filter}" \\
        --filter-name "LowQD"          --filter-expression "QD < ${params.gatk_vf_qd_filter}" \\
        --filter-name "HighFS"         --filter-expression "FS > ${params.gatk_vf_fs_filter}" \\
        --filter-name "LowMQ"          --filter-expression "MQ < ${params.gatk_vf_mq_filter}" \\
        --filter-name "HighSOR"        --filter-expression "SOR > ${params.gatk_vf_sor_filter}" \\
        --filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum < ${params.gatk_vf_read_pos_filter}" \\
        --filter-name "LowBaseQRankSum"  --filter-expression "BaseQRankSum < ${params.gatk_vf_baseq_filter}" \\
        -O "${sample_id}_filtered.vcf.gz"

    gatk --java-options "-Xmx${avail_mem}g" IndexFeatureFile -I "${sample_id}_filtered.vcf.gz"


    if [ ! -s "${sample_id}_filtered.vcf.gz" ] || [ ! -s "${sample_id}_filtered.vcf.gz.tbi" ]; then
        echo "Error: Filtered VCF or index is empty for ${sample_id}" >&2
        exit 1
    fi

    # Capture version
    gatk_version=\$(gatk --version | head -n 1)
cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
