process STAR_ALIGNMENT {
    tag { sample_id }

    cpus params.get('star_alignment_cpus', 12)
    memory params.get('star_alignment_memory', '32 GB')
    time params.get('star_alignment_time', '6h')

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/star_output", mode: "copy"

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2), val(strandedness)
    path star_index_dir
    path gtf_file

    output:
    tuple val(sample_id), 
          path("${sample_id}_Aligned.sortedByCoord.out.bam"), 
          path("${sample_id}_Log.final.out"), 
          path("${sample_id}_Log.out"), 
          path("${sample_id}_Log.progress.out"),
          path("${sample_id}_Chimeric.out.sam"),
          path("${sample_id}_SJ.out.tab"),
          val(strandedness)

    script:
    """
    echo "Running STAR Two-Pass Alignment for Sample: ${sample_id}"

    THREADS=${task.cpus}
    RAM_LIMIT=32000000000  # 32GB for sorting

    # **First Pass: Discover Splice Junctions**
    STAR --genomeDir ${star_index_dir} \
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
         --outSAMattributes NH HI AS nM MD NM \
         --outFileNamePrefix ${sample_id}_pass1_ \
         --limitBAMsortRAM \$RAM_LIMIT

    # **Second Pass: Final Alignment with Junctions from First Pass**
    STAR --genomeDir ${star_index_dir} \
         --readFilesIn ${trimmed_r1} ${trimmed_r2} \
         --runThreadN \$THREADS \
         --readFilesCommand zcat \
         --sjdbFileChrStartEnd ${sample_id}_pass1_SJ.out.tab \
         --sjdbGTFfile ${gtf_file} \
         --twopassMode None \
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
         --limitBAMsortRAM \$RAM_LIMIT \
         --outSAMattrRGline ID:$sample_id LB:library PL:illumina PU:machine SM:$sample_id \
         --outSAMunmapped Within \
         --quantMode TranscriptomeSAM GeneCounts \
         --outFileNamePrefix ${sample_id}_

    """
}