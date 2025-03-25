process GATK_HAPLOTYPE_CALLER {
    tag { sample_id }
    label 'process_high'

    container params.gatk_container
    publishDir params.haplotype_caller_outdir, mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai)
    path(genome)
    path(genome_index)
    path(genome_dict)
    path(known_sites_vcf)
    path(known_sites_vcf_index)

    output:
    tuple val(sample_id), val(strandedness), 
          path("output_${bam.baseName}.vcf.gz"), 
          path("output_${bam.baseName}.vcf.gz.tbi"), emit: vcf_output
    path("versions.yml"), emit: versions

    script:
    """
    THREADS=${task.cpus}

    echo "Running GATK HaplotypeCaller for sample: ${sample_id}"

    # Run HaplotypeCaller with RNA-seq optimizations
    gatk HaplotypeCaller \\
        --native-pair-hmm-threads \${THREADS} \\
        --reference "${genome}" \\
        --output "output_${bam.baseName}.vcf.gz" \\
        -I "${bam}" \\
        --standard-min-confidence-threshold-for-calling 10.0 \\
        --min-base-quality-score 10 \\
        --output-mode EMIT_VARIANTS_ONLY \\
        --dont-use-soft-clipped-bases true \\
        --disable-read-filter NotDuplicateReadFilter \\
        --dbsnp "${known_sites_vcf}" \\
        --verbosity INFO

    # Check output
    if [ ! -s "output_${bam.baseName}.vcf.gz" ]; then
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
