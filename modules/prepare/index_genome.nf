

//Create a fasta Genome Index

process PREPARE_GENOME_SAMTOOLS {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.outdir}/genome/index", mode: "copy"

    input:
    path genome

    output:
    path "${genome}.fai"

    script:
    """
    samtools faidx ${genome}
    """
}

