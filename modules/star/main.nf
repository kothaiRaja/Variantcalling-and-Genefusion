process STAR_ALIGNMENT {
    tag { sample_id }
	label 'process_high'


    container params.star_container
    publishDir params.star_outdir, mode: "copy"

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2), val(strandedness)
    path star_index_dir
    path gtf_file

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sample_id), val(strandedness), path("${sample_id}_Log.final.out"), emit: log_final
    tuple val(sample_id), val(strandedness), path("${sample_id}_Log.out"), emit: log_out
    tuple val(sample_id), val(strandedness), path("${sample_id}_Log.progress.out"), emit: log_progress
    tuple val(sample_id), val(strandedness), path("${sample_id}_Chimeric.out.sam"), optional: true, emit: chimeric_sam
    tuple val(sample_id), val(strandedness), path("${sample_id}_SJ.out.tab"), optional: true, emit: junctions
	path("versions.yml"), emit: versions


    script:
    def extra_args = params.get('star_extra_args', '') 
    def out_sam_type = "BAM SortedByCoordinate"
	def out_sam_attr = "--outSAMattrRGline ID:${sample_id} LB:library PL:${params.get('seq_platform', 'ILLUMINA')} PU:machine SM:${sample_id} CN:${params.get('seq_center', 'Unknown')}"
	def RAM_LIMIT = task.memory.toMega() * 1000000
	// Determine if we should use --outSAMstrandField intronMotif
    def strand_option = ""
    if (strandedness == 'unstranded') {
        strand_option = "--outSAMstrandField intronMotif"
    } 


    """
    echo "Running STAR Alignment for Sample: ${sample_id}"

    THREADS=${task.cpus}
    RAM_LIMIT=${RAM_LIMIT}



    # **First Pass: Discover Splice Junctions**
    STAR --genomeDir ${star_index_dir} \
         --readFilesIn ${trimmed_r1} ${trimmed_r2} \
         --runThreadN \$THREADS \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin ${params.get('star_alignSJoverhangMin', 8)} \
         --alignSJDBoverhangMin ${params.get('star_alignSJDBoverhangMin', 1)} \
         --outFilterMismatchNmax ${params.get('star_outFilterMismatchNmax', 999)} \
         --outFilterMatchNmin ${params.get('star_outFilterMatchNmin', 16)} \
         --outFilterMatchNminOverLread ${params.get('star_outFilterMatchNminOverLread', 0.3)} \
         --outFilterScoreMinOverLread ${params.get('star_outFilterScoreMinOverLread', 0.3)} \
         --outSAMattributes NH HI AS nM MD NM \
         --outFileNamePrefix ${sample_id}_pass1_ \
         --limitBAMsortRAM \$RAM_LIMIT \
		 $strand_option \
         

    # **Second Pass: Final Alignment**
    STAR --genomeDir ${star_index_dir} \
         --readFilesIn ${trimmed_r1} ${trimmed_r2} \
         --runThreadN \$THREADS \
         --readFilesCommand zcat \
         --sjdbFileChrStartEnd ${sample_id}_pass1_SJ.out.tab \
         --sjdbGTFfile ${gtf_file} \
         --twopassMode None \
         --outFilterType BySJout \
         --alignSJoverhangMin ${params.get('star_alignSJoverhangMin', 8)} \
         --alignSJDBoverhangMin ${params.get('star_alignSJDBoverhangMin', 1)} \
         --outFilterMismatchNmax ${params.get('star_outFilterMismatchNmax', 999)} \
         --outFilterMatchNmin ${params.get('star_outFilterMatchNmin', 16)} \
         --outFilterMatchNminOverLread ${params.get('star_outFilterMatchNminOverLread', 0.3)} \
         --outFilterScoreMinOverLread ${params.get('star_outFilterScoreMinOverLread', 0.3)} \
         --chimSegmentMin ${params.get('star_chimSegmentMin', 10)} \
         --chimJunctionOverhangMin ${params.get('star_chimJunctionOverhangMin', 10)} \
         --chimScoreJunctionNonGTAG ${params.get('star_chimScoreJunctionNonGTAG', -4)} \
         --chimScoreMin ${params.get('star_chimScoreMin', 1)} \
         --chimOutType ${params.get('star_chimOutType', 'WithinBAM SeparateSAMold')} \
         --chimScoreDropMax ${params.get('star_chimScoreDropMax', 50)} \
         --chimScoreSeparation ${params.get('star_chimScoreSeparation', 10)} \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM \$RAM_LIMIT \
         --outSAMunmapped Within \
		 --quantMode TranscriptomeSAM GeneCounts \
		 --outFileNamePrefix ${sample_id}_ \
         ${out_sam_attr} \
		 $strand_option \
		$extra_args
		
	#capture the versions
	star_version=\$(STAR --version)
	cat <<EOF > versions.yml
	STAR:
	  version: "${star_version}"
	EOF


    """
}
