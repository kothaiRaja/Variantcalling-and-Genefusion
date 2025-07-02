process GATK_HAPLOTYPE_CALLER {
    tag { "${meta.id}_${task.process}" }

    label 'process_high'

    container params.gatk_container
    publishDir params.haplotype_caller_outdir, mode: "copy"

    input:
    tuple val(meta), path(bam), path(bai), path(interval)
    path(genome)
    path(genome_index)
    path(genome_dict)
    path(known_sites_vcf)
    path(known_sites_vcf_index)

    output:
	tuple val(meta), path("output_${meta.id}_split_${interval.baseName}.vcf.gz"), path("output_${meta.id}_split_${interval.baseName}.vcf.gz.tbi"), emit: vcf_output
	path("versions.yml"), emit: versions

    script:
    def sample_id = meta.id
    def avail_mem = task.memory ? task.memory.giga : 3
    def interval_command = interval ? "--intervals ${interval}" : ""

    """
    THREADS=${task.cpus}

    echo "Running GATK HaplotypeCaller for sample: ${sample_id}"

    gatk --java-options "-Xmx${avail_mem}g" HaplotypeCaller \\
        --native-pair-hmm-threads \${THREADS} \\
        --R "${genome}" \\
        -I "${bam}" \\
        --output "output_${sample_id}_split_${interval.baseName}.vcf.gz" \\
        --standard-min-confidence-threshold-for-calling 30.0 \\
        --output-mode EMIT_VARIANTS_ONLY \\
        --dont-use-soft-clipped-bases \\
        --create-output-variant-index true \\
        --dbsnp "${known_sites_vcf}" \\
        --verbosity INFO \\
        ${interval_command}

    if [ ! -s "output_${sample_id}_split_${interval.baseName}.vcf.gz" ]; then
        echo "Error: VCF file not generated!" >&2
        exit 1
    fi

    gatk_version=\$(gatk --version | head -n 1)
    cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF

    echo "HaplotypeCaller finished for ${sample_id}"
    """
}
