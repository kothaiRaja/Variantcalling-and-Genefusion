process RNASEQ_MAPPING_STAR {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/star", mode: "copy"

    input:
    path genomeDir
    tuple val(replicateId), path(reads)

    output:
    tuple val(replicateId), path("Aligned.*.sortedByCoord.out.bam")


    script:
    """
    # Align reads to genome
    STAR --genomeDir $genomeDir \
         --readFilesIn ${reads[0]} ${reads[1]} \
         --runThreadN ${task.cpus} \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outSAMtype BAM SortedByCoordinate \
		 --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:GM12878 \
         --outFileNamePrefix Aligned.
    """
}