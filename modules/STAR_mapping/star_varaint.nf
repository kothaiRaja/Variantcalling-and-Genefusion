process STAR_ALIGNMENT {
    tag { sample_id }

   container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/STAR", mode: "copy"

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2)
    path star_index_dir
    

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    """
    # ngs-nf-dev Align reads to genome
  STAR --genomeDir $star_index_dir \
       --readFilesIn ${trimmed_r1} ${trimmed_r2}  \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999

  # Run 2-pass mapping (improve alignmets using table of splice junctions and create a new index)  
  STAR --genomeDir $star_index_dir \
       --readFilesIn ${trimmed_r1} ${trimmed_r2} \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --sjdbFileChrStartEnd SJ.out.tab \
	   --outFileNamePrefix ${sample_id}_ \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$sample_id LB:library PL:illumina PU:machine SM:GM12878
    """
}