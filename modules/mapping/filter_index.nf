process SAMTOOLS_FILTER_INDEX {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/star/index", mode: "copy"

    input:
    tuple val(replicateId), path(bam)

    output:
    tuple val(replicateId), path("${bam.baseName}.uniq.bam"), path("${bam.baseName}.uniq.bam.bai")
	
    script:
    """
    # Filter for unique alignments and create a new BAM file
    (samtools view -H $bam; samtools view $bam | grep -w 'NH:i:1') \
    | samtools view -Sb - > ${bam.baseName}.uniq.bam

    # Index the BAM file
    samtools index ${bam.baseName}.uniq.bam
    """
}