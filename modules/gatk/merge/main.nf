process GATK_MERGEVCFS {
    tag { "${meta.id}_${task.process}" }

    label 'process_medium'

    container params.gatk_container
    publishDir params.merged_vcf_outdir, mode: "copy"

    input:
    tuple val(meta), path(vcf_list), path(tbi_list)

    output:
    tuple val(meta.id), 
          path("merged_${meta.id}.vcf.gz"), 
          path("merged_${meta.id}.vcf.gz.tbi"), emit: merged_vcf
    path("versions.yml"), emit: versions

 script:
def avail_mem = task.memory ? task.memory.giga : 3
def vcf_inputs = vcf_list.collect { "-I \"${it}\"" }.join(" \\\n")

"""
echo "Merging VCFs for sample: ${meta.id}"

gatk --java-options "-Xmx${avail_mem}g" MergeVcfs \\
${vcf_inputs} \\
-O "merged_${meta.id}.vcf.gz"

gatk --java-options "-Xmx${avail_mem}g" IndexFeatureFile -I "merged_${meta.id}.vcf.gz"

gatk_version=\$(gatk --version | head -n 1)
cat <<EOF > versions.yml
"${task.process}":
  gatk: \${gatk_version}
EOF
"""



}
