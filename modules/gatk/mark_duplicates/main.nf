process GATK_MARK_DUPLICATES {
   
    tag { "${meta.id}_${task.process}" }

	label 'process_high'
	
	container params.gatk_container
    publishDir params.markduplicates_outdir, mode: "copy"

    input:
    tuple val(meta), path(sorted_bam), path(sorted_bam_index)

    output:
    tuple val(meta), path("${meta.id}_marked_duplicates.bam"), path("${meta.id}_marked_duplicates.bai"), emit: marked_bams_bai
    tuple val(meta), path("${meta.id}_dup_metrics.txt"), emit: marked_bams_bai_metrics
    path("versions.yml"), emit: versions

    script:
    def sample_id   = meta.id
    def avail_mem   = task.memory ? task.memory.giga : 3
	
    """
	
	
    TTHREADS=${task.cpus}

    gatk --java-options "-Xmx${avail_mem}g" MarkDuplicates \\
        -I ${sorted_bam} \\
        -O ${sample_id}_marked_duplicates.bam \\
        -M ${sample_id}_dup_metrics.txt \\
        --CREATE_INDEX true \\
        --ASSUME_SORTED true \\
        --REMOVE_DUPLICATES ${params.remove_duplicates ? 'true' : 'false'} \\
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

    # Capture version
    gatk_version=\$(gatk --version | head -n 1)

cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF

    """

}