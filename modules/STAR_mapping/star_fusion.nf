process STAR_ALIGN_FUSION {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
	publishDir "${params.outdir}/Star_fusion", mode: 'copy'
    
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
    STAR --genomeDir $star_index_dir \
     --readFilesIn ${trimmed_r1} ${trimmed_r2} \
     --runThreadN 8 \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 10 \
     --outFilterMatchNmin 16 \
     --outFilterMatchNminOverLread 0.3 \
     --outFilterScoreMinOverLread 0.3 \
     --chimSegmentMin 12 \
     --chimJunctionOverhangMin 15 \
     --chimOutType WithinBAM SeparateSAMold \
	 --chimMultimapNmax 10 \
     --chimScoreDropMax 50 \
     --chimScoreSeparation 10 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes NH HI AS nM MD NM \
	 --limitBAMsortRAM 32000000000 \
     --outFileNamePrefix ${sample_id}_

    """
}