process GATK_HAPLOTYPE_CALLER {
    tag { sample_id }

     container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/haplotype_caller", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai), path(recal_table)
	path (genome)
	path (genome_index)
	path (genome_dict)
	path (interval_list)
	path (known_variants)
	path (known_variants_index)

    output:
    tuple val(sample_id), path("output_${sample_id}.vcf.gz"), path("output_${sample_id}.vcf.gz.tbi") 

    script:
	def intervals_args = interval_list.collect { "--intervals ${it}" }.join(' ')

	
    """
    gatk HaplotypeCaller \
    --native-pair-hmm-threads ${task.cpus} \
    --reference ${genome} \
    --output output_${sample_id}.vcf.gz \
    -I $bam \
    --standard-min-confidence-threshold-for-calling 10.0 \
    --dont-use-soft-clipped-bases true \
    --min-base-quality-score 10 \
    --output-mode EMIT_VARIANTS_ONLY \
    ${intervals_args} \
    --verbosity DEBUG \
	--dbsnp ${known_variants}
	
	# Ensure output VCF file is not empty
    if [ ! -s output_${sample_id}.vcf.gz ]; then
        echo "Error: Output VCF is empty for ${sample_id}" >&2
        exit 1
    fi


    """
}
