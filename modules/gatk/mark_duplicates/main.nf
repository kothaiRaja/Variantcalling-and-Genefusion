process GATK_MARK_DUPLICATES {
    tag { sample_id }
	
	cpus params.get('gatk_mark_duplicates_cpus', 8)
    memory params.get('gatk_mark_duplicates_memory', '16 GB')
    time params.get('gatk_mark_duplicates_time', '2h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/dedup_bam", mode: "copy"

    input:
    tuple val(sample_id),val(strandedness), path(sorted_bam), path(sorted_bam_index)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_marked_duplicates.bam"), path("${sample_id}_marked_duplicates.bai"), emit: marked_bams_bai
	tuple val(sample_id), val(strandedness), path("${sample_id}_dup_metrics.txt"), emit: marked_bams_bai_metrics

    script:
    """
    THREADS=${task.cpus}

    gatk MarkDuplicates \
        -I ${sorted_bam} \
        -O ${sample_id}_marked_duplicates.bam \
        -M ${sample_id}_dup_metrics.txt \
        --CREATE_INDEX true \
		--REMOVE_DUPLICATES ${params.remove_duplicates ? 'true' : 'false'} \
        --VALIDATION_STRINGENCY ${params.validation_stringency ?: 'LENIENT'}
		
	# Check if output BAM is generated
    if [ ! -s ${sample_id}_marked_duplicates.bam ]; then
        echo "Error: Marked duplicates BAM file not generated for ${sample_id}" >&2
        exit 1
    fi

    # Check if metrics file is generated
    if [ ! -s ${sample_id}_dup_metrics.txt ]; then
        echo "Error: Duplicate metrics file not generated for ${sample_id}" >&2
        exit 1
    fi
    """
}