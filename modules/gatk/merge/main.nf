process GATK_MERGEVCFS {
    tag { sample_id }
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

    """
    echo "Merging VCFs for sample: ${sample_id}"

    gatk MergeVcfs \\
    ${vcf_inputs} \\
    -O "merged_${sample_id}.vcf.gz"

    gatk IndexFeatureFile -I "merged_${sample_id}.vcf.gz"

    # Capture GATK version
    gatk_version=\$(gatk --version | awk '{print \$2}')
    cat <<EOF > versions.yml
    "${task.process}":
      gatk: "\${gatk_version}"
    EOF
    """
}

