process STAR_ALIGN_FUSION {
    tag { "test" }

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
	publishDir "${params.outdir}/star", mode: "copy"
    
	input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2)
    path star_index_dir            
    path gtf                                          

    output:
    tuple val(sample_id), path('*Log.final.out'), emit: log_final
    tuple val(sample_id), path('*Log.out'), emit: log_out
    tuple val(sample_id), path('*Aligned.sortedByCoord.out.bam'), emit: bam_sorted
    tuple val(sample_id), path('*Chimeric.out.sam'), emit: chimeric_sam
    tuple val(sample_id), path('*Log.progress.out'), emit: log_progress
    tuple val(sample_id), path('*SJ.out.tab'), emit: splice_junctions

    
    script:
    """
    # First-pass alignment
    STAR --genomeDir $star_index_dir \
         --readFilesIn ${trimmed_r1} ${trimmed_r2} \
         --runThreadN $task.cpus \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --chimSegmentMin 10 \
         --chimJunctionOverhangMin 15 \
         --chimOutType WithinBAM SeparateSAMold \
         --outSAMtype BAM Unsorted \
		 --outFileNamePrefix ${sample_id}_pass1_

    # Second-pass alignment (uses splice junction table from first pass)
    STAR --genomeDir $star_index_dir \
         --readFilesIn ${trimmed_r1} ${trimmed_r2} \
         --runThreadN $task.cpus \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --sjdbFileChrStartEnd ${sample_id}_pass1_SJ.out.tab \
         --chimSegmentMin 10 \
         --chimJunctionOverhangMin 15 \
         --chimOutType WithinBAM SeparateSAMold \
         --outFileNamePrefix ${sample_id}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattrRGline ID:${sample_id} LB:library PL:illumina PU:machine SM:${sample_id}
    """

}
