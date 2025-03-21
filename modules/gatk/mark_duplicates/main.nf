process GATK_MARK_DUPLICATES {
   
    tag { sample_id }
	label 'process_high'
	
	container params.gatk_container
    publishDir params.markduplicates_outdir, mode: "copy"

    input:
    tuple val(sample_id),val(strandedness), path(sorted_bam), path(sorted_bam_index)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_marked_duplicates.bam"), path("${sample_id}_marked_duplicates.bai"), emit: marked_bams_bai
	tuple val(sample_id), val(strandedness), path("${sample_id}_dup_metrics.txt"), emit: marked_bams_bai_metrics
	path("versions.yml"), emit: versions

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
	
	#  capture version
	gatk_version=\$(gatk --version | grep -Eo '[0-9.]+' | head -n 1)

	cat <<EOF > versions.yml
	"GATK_MARK_DUPLICATES":
	  gatk: "\${gatk_version}"
	EOF
    """
}