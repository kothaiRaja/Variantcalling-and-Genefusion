process RNASEQ_GATK_SPLITNCIGAR {
    tag "$replicateId"
    label 'mem_large'
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/SNG", mode: 'copy'

    input:
    path genome
    path index
    path genome_dict
    tuple val(replicateId), path(bam)

    output:
    tuple val(replicateId), path('split.bam')

    script:
    """
    gatk SplitNCigarReads \
        -R $genome \
        -I $bam \
        --refactor-cigar-string \
        -O split.bam
    """
}