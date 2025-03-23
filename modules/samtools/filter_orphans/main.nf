process SAMTOOLS_FILTER_ORPHANS {
    tag { sample_id }
	label 'process_low'

    container params.samtools_container
    publishDir params.samtools_filter_outdir, mode: "copy"

    input:
    tuple val(sample_id),val(strandedness), path(sorted_bam), path(bai)

    output:
    tuple val(sample_id),val(strandedness), path("${sample_id}_filtered.bam"), path("${sample_id}_filtered.bam.bai"), emit: filtered_sorted_bams
	path("versions.yml"), emit: versions

    script:
    """
    # Filter out orphan reads (retain properly paired reads)
    samtools view -h -F 3844 ${sorted_bam} | samtools view -b - > ${sample_id}_filtered.bam


	# Index the filtered BAM file
    samtools index ${sample_id}_filtered.bam
	
	# Capture samtools version
    samtools_version=\$(samtools --version | head -n 1 | awk '{print \$2}')
    cat <<EOF > versions.yml
    "${task.process}":
      samtools: "\${samtools_version}"
    EOF
    """
}
