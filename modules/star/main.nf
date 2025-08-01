process STAR_ALIGNMENT {
    tag { "${meta.id}_${task.process}" }

    label 'process_high'

    container params.star_container
    publishDir params.star_outdir, mode: "copy"

    input:
    tuple val(meta), path(read1), path(read2)
    path star_index_dir
    path gtf_file

    output:
    tuple val(meta), path("${meta.id}_Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}_Aligned.toTranscriptome.out.bam"), emit: transcriptome_bam
    tuple val(meta), path("${meta.id}_Chimeric.out.sam"), optional: true, emit: chimeric_sam
	tuple val(meta), path("${meta.id}_Chimeric.out.junction"), optional: true, emit: chimeric_junction
    tuple val(meta), path("${meta.id}_SJ.out.tab"), emit: junctions
    tuple val(meta), path("${meta.id}_ReadsPerGene.out.tab"), emit: gene_counts
    tuple val(meta), path("${meta.id}_Log.final.out"), emit: log_final
    tuple val(meta), path("${meta.id}_Log.out"), emit: log_out
    tuple val(meta), path("${meta.id}_Log.progress.out"), emit: log_progress
    path("versions.yml"), emit: versions

    script:
    def sample_id     = meta.id
    def strandedness  = meta.strandedness
    def extra_args    = params.get('star_extra_args', '') 
    def out_sam_attr  = "--outSAMattrRGline ID:${sample_id} LB:library PL:${params.get('seq_platform', 'ILLUMINA')} PU:machine SM:${sample_id} CN:${params.get('seq_center', 'Unknown')}"
    def RAM_LIMIT     = task.memory.toMega() * 1000000
    def strand_option = (strandedness == 'unstranded') ? "--outSAMstrandField intronMotif" : ""
    def sjdb_overhang = params.read_length ? "--sjdbOverhang ${params.read_length - 1}" : ""
    def gtf_flag      = params.star_ignore_sjdbgtf ? "" : "--sjdbGTFfile ${gtf_file}"
    def intron_min    = "--alignIntronMin ${params.get('star_alignIntronMin', 20)}"
    def intron_max    = "--alignIntronMax ${params.get('star_alignIntronMax', 1000000)}"

    """
    echo "Running STAR Alignment for Sample: ${sample_id}"

    THREADS=${task.cpus}
    RAM_LIMIT=${RAM_LIMIT}

    STAR --genomeDir ${star_index_dir} \
         --readFilesIn ${read1} ${read2} \
         --readFilesCommand zcat \
         --runThreadN \$THREADS \
         --twopassMode Basic \
         $sjdb_overhang \
         $gtf_flag \
         $intron_min \
         $intron_max \
         --outFilterType BySJout \
         --alignSJoverhangMin ${params.get('star_alignSJoverhangMin', 8)} \
         --alignSJDBoverhangMin ${params.get('star_alignSJDBoverhangMin', 1)} \
         --outFilterMismatchNmax ${params.get('star_outFilterMismatchNmax', 999)} \
         --outFilterMatchNmin ${params.get('star_outFilterMatchNmin', 16)} \
         --outFilterMatchNminOverLread ${params.get('star_outFilterMatchNminOverLread', 0.3)} \
         --outFilterScoreMinOverLread ${params.get('star_outFilterScoreMinOverLread', 0.3)} \
         --outFilterMismatchNoverReadLmax ${params.get('star_mismatchNoverLmax', 0.04)} \
         --outSAMmapqUnique ${params.get('star_outSAMmapqUnique', 60)} \
         --chimSegmentMin ${params.get('star_chimSegmentMin', 10)} \
         --chimJunctionOverhangMin ${params.get('star_chimJunctionOverhangMin', 10)} \
         --chimScoreJunctionNonGTAG ${params.get('star_chimScoreJunctionNonGTAG', -4)} \
         --chimScoreMin ${params.get('star_chimScoreMin', 1)} \
         --chimOutType WithinBAM HardClip \
         --chimScoreDropMax ${params.get('star_chimScoreDropMax', 50)} \
         --chimScoreSeparation ${params.get('star_chimScoreSeparation', 10)} \
         --limitBAMsortRAM \$RAM_LIMIT \
         --outSAMunmapped Within \
         --quantMode TranscriptomeSAM GeneCounts \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${sample_id}_ \
         ${out_sam_attr} \
         $strand_option \
         $extra_args

    # Capture version
    star_version=\$(STAR --version)

cat <<EOF > versions.yml
"${task.process}":
  STAR: "\${star_version}"
EOF
    """
}
