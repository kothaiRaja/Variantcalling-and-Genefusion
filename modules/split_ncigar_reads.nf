process SPLIT_NCIGAR_READS {
    tag { sample_id }
	
	cpus params.get('split_ncigar_reads_cpus', 8)
    memory params.get('split_ncigar_reads_memory', '16 GB')
    time params.get('split_ncigar_reads_time', '3h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/split_ncigar", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai),val(strandedness)
	path(genome_fasta)
	path (index)
	path (genome_dict)

    output:
    tuple val(sample_id), 
          path("${sample_id}_split.bam"), 
          path("${sample_id}_split.bai"),
		  val(strandedness)

    script:
    """
   gatk SplitNCigarReads \
    -R ${genome_fasta} \
    -I ${bam} \
    -O ${sample_id}_split.bam \
    --create-output-bam-index true \
    --skip-mapping-quality-transform false \
    --max-mismatches-in-overhang 1 \
    --max-bases-in-overhang 50 \
	--process-secondary-alignments true
		
	"""
}