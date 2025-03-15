process SAMTOOLS_STATS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_stats_report.*"

    input:
    tuple val(sample_id), val(strandedness), path(sorted_bam), path(bai)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_stats_report.txt")
	
    script:
    """
    # Run samtools stats on the sorted BAM file
    samtools stats ${sorted_bam} > ${sample_id}_stats_report.txt
    """
}