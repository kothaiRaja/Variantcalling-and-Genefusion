process SAMTOOLS_FILTER_ORPHANS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/filtered_bam", mode: "copy"

    input:
    tuple val(sample_id),val(strandedness), path(sorted_bam), path(bai)

    output:
    tuple val(sample_id),val(strandedness), path("${sample_id}_filtered.bam"), path("${sample_id}_filtered.bam.bai"), emit: filtered_sorted_bams

    script:
    """
    # Filter out orphan reads (retain properly paired reads)
    samtools view -h -F 3844 ${sorted_bam} | samtools view -b - > ${sample_id}_filtered.bam


	# Index the filtered BAM file
    samtools index ${sample_id}_filtered.bam
    """
}
