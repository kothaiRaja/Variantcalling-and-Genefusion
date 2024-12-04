process ARRIBA {
    tag { sample_id }
    
    container "https://depot.galaxyproject.org/singularity/arriba%3A2.4.0--hdbdd923_3"
    publishDir "${params.outdir}/arriba", mode: "copy"

    input:
	
	tuple val(sample_id),path(log_final)
	tuple val(sample_id),path(log_out)
	tuple val(sample_id),path(bam)
	tuple val(sample_id),path(chimeric_sam)
	tuple val(sample_id), path(log_progress) 
    tuple val(sample_id), path(splice_junctions) 
    path fasta                                     
    path gtf                                       
    path blacklist                                 
    path known_fusions                             

    output:
    tuple val(sample_id), path("${sample_id}.fusions.tsv")          , emit: fusions
    tuple val(sample_id), path("${sample_id}.fusions.discarded.tsv"), emit: fusions_discarded

    script:
    """
    echo "Running Arriba with the following command:"
    echo arriba \\
         -x $bam \\
         -c ${sample_id}.Chimeric.out.sam \\
         -a $fasta \\
         -g $gtf \\
         -b $blacklist \\
         -k $known_fusions \\
         -o ${sample_id}.fusions.tsv \\
         -O ${sample_id}.fusions.discarded.tsv

    arriba \\
         -x $bam \\
         -c ${sample_id}_Chimeric.out.sam \\
         -a $fasta \\
         -g $gtf \\
         -b $blacklist \\
         -k $known_fusions \\
         -o ${sample_id}.fusions.tsv \\
         -O ${sample_id}.fusions.discarded.tsv
    """
}
