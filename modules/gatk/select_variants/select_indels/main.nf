process SELECT_INDELs {
    tag { "${sample_id}_${task.process}" }

    label 'process_medium'

    container params.gatk_container
    publishDir params.indels_select_outdir, mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path genome
    path index
    path dict

    output:
    tuple val(sample_id), path("${sample_id}_indels.vcf.gz"), path("${sample_id}_indels.vcf.gz.tbi"), emit: selected_indels
    path("versions.yml"), emit: versions

    script:
	def avail_mem = 3
if (task.memory) {
    avail_mem = task.memory.giga
} else {
    log.info '[GATK SelectVariants] No memory set — defaulting to 3GB.'
}

	
    """
    echo "Selecting INDELs for sample: ${sample_id}"

    gatk --java-options "-Xmx${avail_mem}g" SelectVariants \\
        -R "${genome}" \\
        -V "${vcf_file}" \\
        --select-type-to-include INDEL \\
        -O "${sample_id}_indels.vcf.gz"

    gatk --java-options "-Xmx${avail_mem}g" IndexFeatureFile -I "${sample_id}_indels.vcf.gz"


    # Capture GATK version
    gatk_version=\$(gatk --version | head -n 1)
cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
