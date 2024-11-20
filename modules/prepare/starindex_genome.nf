//Create a genome Index file

process PREPARE_STAR_GENOME_INDEX {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
    publishDir "${params.outdir}/genome/starindex", mode: "copy"

    input:
    path genome

    output:
    path "genome_index" 

    script:
    """
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeDir genome_index \
         --genomeFastaFiles ${genome} \
         --genomeSAindexNbases 11
    """
}

