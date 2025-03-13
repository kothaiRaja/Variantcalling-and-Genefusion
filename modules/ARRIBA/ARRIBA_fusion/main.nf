process ARRIBA {
    tag { sample_id }

    cpus params.get('star_align_fusion_cpus', 12)
    memory params.get('star_align_fusion_memory', '24 GB')
    time params.get('star_align_fusion_time', '5h')

    container "https://depot.galaxyproject.org/singularity/arriba%3A2.4.0--hdbdd923_3"
    publishDir "${params.outdir}/ARRIBA", mode: 'copy'

    input:
    tuple val(sample_id), path(log_final), path(log_out), path(bam), path(chimeric_sam), path(log_progress), path(splice_junctions), val(strandedness)
    path fasta                                     
    path gtf                                       
    path blacklist                                 
    path known_fusions

    output:
    tuple val(sample_id), path("${sample_id}.fusions.tsv"), val(strandedness), emit: fusions
    tuple val(sample_id), path("${sample_id}.fusions.discarded.tsv"), val(strandedness), emit: fusions_discarded

    script:
    """
    echo "Running Arriba for Sample: ${sample_id}"

    arriba \
         -x $bam \
         -c $chimeric_sam \
         -a $fasta \
         -g $gtf \
         -b $blacklist \
         -k $known_fusions \
         -o ${sample_id}.fusions.tsv \
         -O ${sample_id}.fusions.discarded.tsv

    echo "Arriba finished for Sample: ${sample_id}"
    """
}