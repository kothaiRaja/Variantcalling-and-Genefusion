process GATK_HAPLOTYPE_CALLER {
    tag { sample_id }
	
	cpus params.get('gatk_haplotype_caller_cpus', 12)
    memory params.get('gatk_haplotype_caller_memory', '24 GB')
    time params.get('gatk_haplotype_caller_time', '8h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/haplotype_caller", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai) ,val(strandedness)
	path (genome)
	path (genome_index)
	path (genome_dict)
	path (interval_list)

    output:
    tuple val(sample_id), path("output_${sample_id}.vcf.gz"), path("output_${sample_id}.vcf.gz.tbi"), val(strandedness)

    script:
	def intervals_args = interval_list.collect { "--intervals ${it}" }.join(' ')
	
    """
	
	THREADS=${task.cpus}
	
    # Validate Inputs
    if [ ! -s ${bam} ]; then
        echo "Error: BAM file not found or empty." >&2
        exit 1
    fi
    if [ ! -s ${genome} ]; then
        echo "Error: Reference genome not found or empty." >&2
        exit 1
    fi

    # Run HaplotypeCaller
    gatk HaplotypeCaller \
        --native-pair-hmm-threads \$THREADS \
        --reference ${genome} \
        --output output_${sample_id}.vcf.gz \
        -I ${bam} \
        --standard-min-confidence-threshold-for-calling 30.0 \
        --min-base-quality-score 10 \
        --output-mode EMIT_VARIANTS_ONLY \
        ${intervals_args} \
        --verbosity INFO
    """

  
}
