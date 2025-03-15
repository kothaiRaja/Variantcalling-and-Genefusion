process SAMTOOLS_SORT_INDEX {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/sorted_bam", mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam_sorted

    script:
    """
    # Sort the BAM file
    samtools sort -o ${sample_id}_sorted.bam ${bam}

    # Index the sorted BAM file
    samtools index ${sample_id}_sorted.bam
    """
}

