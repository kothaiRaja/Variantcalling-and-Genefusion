process GATK_HAPLOTYPE_CALLER {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/haplotype_caller", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai), path(table)
	path (genome)
	path (genome_index)
	path (genome_dict)
	path (interval_list)

    output:
    tuple val(sample_id), path("output_${sample_id}.vcf.gz"), path("output_${sample_id}.vcf.gz.tbi") 

    script:
    """
    gatk HaplotypeCaller \
    --native-pair-hmm-threads ${task.cpus} \
    --reference ${genome} \
    --output output_${sample_id}.vcf.gz \
    -I $bam \
    --standard-min-confidence-threshold-for-calling 5.0 \
    --dont-use-soft-clipped-bases true \
    --min-base-quality-score 10 \
	--emit-ref-confidence GVCF \
	--output-mode EMIT_VARIANTS_ONLY \
    --intervals $interval_list

    """
}