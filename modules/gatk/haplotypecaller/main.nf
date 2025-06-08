process GATK_HAPLOTYPE_CALLER {
    tag { "${sample_id}_${task.process}" }

    label 'process_high'

    container params.gatk_container
    publishDir params.haplotype_caller_outdir, mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai), path(interval)
    path(genome)
    path(genome_index)
    path(genome_dict)
    path(known_sites_vcf)
    path(known_sites_vcf_index)

	output:
	tuple val(sample_id), val(strandedness), 
		path("output_${sample_id}_split_${interval.baseName}.vcf.gz"), 
		path("output_${sample_id}_split_${interval.baseName}.vcf.gz.tbi"), emit: vcf_output
	path("versions.yml"), emit: versions

    script:
	def interval_command = interval ? "--intervals ${interval}" : ""
	def avail_mem = 3
if (task.memory) {
    avail_mem = task.memory.giga
} else {
    log.info '[GATK HaplotypeCaller] No memory set — defaulting to 3GB.'
}


    """
    THREADS=${task.cpus}

    echo "Running GATK HaplotypeCaller for sample: ${sample_id}"

    # Run HaplotypeCaller with RNA-seq optimizations
    gatk --java-options "-Xmx${avail_mem}g" HaplotypeCaller \\
		--native-pair-hmm-threads \${THREADS} \\
        --R "${genome}" \\
        --output "output_${sample_id}_split_${interval.baseName}.vcf.gz" \\
        -I "${bam}" \\
        --standard-min-confidence-threshold-for-calling 30.0 \\
        --output-mode EMIT_VARIANTS_ONLY \\
        --dont-use-soft-clipped-bases \\
		--create-output-variant-index true \\
        --dbsnp "${known_sites_vcf}" \\
        --verbosity INFO \\
		${interval_command}

    # Check output
    if [ ! -s "output_${sample_id}_split_${interval.baseName}.vcf.gz" ]; then
        echo "Error: VCF file not generated!" >&2
        exit 1
    fi

    # Capture GATK version
    gatk_version=\$(gatk --version | head -n 1)
	
cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF

    echo "HaplotypeCaller finished for ${sample_id}"
    """
}
