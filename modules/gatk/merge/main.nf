process GATK_MERGEVCFS {
    tag { "${sample_id}_${task.process}" }

    label 'process_medium'

    container params.gatk_container
    publishDir params.merged_vcf_outdir, mode: "copy"

    input:
    tuple val(sample_id), path(vcf_list), path(tbi_list)

    output:
    tuple val(sample_id), 
          path("merged_${sample_id}.vcf.gz"), 
          path("merged_${sample_id}.vcf.gz.tbi"), emit: merged_vcf
    path("versions.yml"), emit: versions

    script:
    def vcf_inputs = vcf_list.collect { "-I \"${it}\"" }.join(" \\\n")
	def avail_mem = 3
if (task.memory) {
    avail_mem = task.memory.giga
} else {
    log.info '[GATK MergeVcfs] No memory set — defaulting to 3GB.'
}


    """
    echo "Merging VCFs for sample: ${sample_id}"

    gatk --java-options "-Xmx${avail_mem}g" MergeVcfs \\
    ${vcf_inputs} \\
    -O "merged_${sample_id}.vcf.gz"

    gatk --java-options "-Xmx${avail_mem}g" IndexFeatureFile -I "merged_${sample_id}.vcf.gz"


    # Capture GATK version
    gatk_version=\$(gatk --version | head -n 1)

cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}

