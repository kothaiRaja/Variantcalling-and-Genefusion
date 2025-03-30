process GATK_VCF_TO_TABLE {
	tag "Table creation"
    label 'process_medium'
    container params.gatk_container
	publishDir params.vcf2table_outdir, mode: "copy"

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)

    output:
    tuple val(sample_id), path("${sample_id}.vcf_table.txt"), emit: vcf_table
    path("versions.yml"), emit: versions

    script:
	
	def avail_mem = 3
if (task.memory) {
    avail_mem = task.memory.giga
} else {
    log.info '[GATK VariantsToTable] No memory set — defaulting to 3GB.'
}

	
    """
    echo "Extracting variant table for sample: ${sample_id}"

    gatk --java-options "-Xmx${avail_mem}g" VariantsToTable \\
        -V ${vcf} \\
        -F CHROM -F POS -F ID -F REF -F ALT \\
        -F QUAL -F FILTER \\
        -F AF -F DP -F ANN -F CSQ \\
        -GF GT -GF DP -GF AD -GF GQ \\
        --output ${sample_id}.vcf_table.txt

    # Capture version
    gatk_version=\$(gatk --version | head -n 1)
cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF

    echo "VCF to table conversion complete for ${sample_id}"
    """
}
