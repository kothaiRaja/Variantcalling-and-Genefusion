//===========================================STAR Fusion=============================================//
process STAR_ALIGN_FUSION {
    tag { sample_id }
	
	cpus params.get('star_align_fusion_cpus', 12)
    memory params.get('star_align_fusion_memory', '24 GB')
    time params.get('star_align_fusion_time', '5h')

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
	publishDir "${params.outdir}/Star_fusion", mode: 'copy'
    
	input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2), val(strandedness)
    path star_index_dir            
    path gtf                                          

    output:
    output:
    tuple val(sample_id), path('*Log.final.out'), val(strandedness), emit: log_final
    tuple val(sample_id), path('*Log.out'), val(strandedness), emit: log_out
    tuple val(sample_id), path('*Aligned.sortedByCoord.out.bam'), val(strandedness), emit: bam_sorted
    tuple val(sample_id), path('*Chimeric.out.sam'), val(strandedness), emit: chimeric_sam
    tuple val(sample_id), path('*Log.progress.out'), val(strandedness), emit: log_progress
    tuple val(sample_id), path('*SJ.out.tab'), val(strandedness), emit: splice_junctions

    
    script:
    """
	THREADS=${task.cpus}
	
    STAR --genomeDir $star_index_dir \
     --readFilesIn ${trimmed_r1} ${trimmed_r2} \
     --runThreadN \$THREADS \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMatchNmin 16 \
     --outFilterMatchNminOverLread 0.3 \
     --outFilterScoreMinOverLread 0.3 \
     --chimSegmentMin 10 \
     --chimJunctionOverhangMin 10 \
	 --chimScoreJunctionNonGTAG -4 \
	 --chimScoreMin 1 \
     --chimOutType WithinBAM SeparateSAMold \
     --chimScoreDropMax 50 \
     --chimScoreSeparation 10 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes NH HI AS nM MD NM \
	 --limitBAMsortRAM 32000000000 \
     --outFileNamePrefix ${sample_id}_

    """
}
