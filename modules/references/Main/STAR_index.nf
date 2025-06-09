process CREATE_STAR_INDEX {
    tag "Create STAR Genome Index"
    label 'process_medium'

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    input:
    path genome_fasta
    path genome_gtf

    output:
    path "STAR_index", emit: star_index

    script:
    def sjdbOverhang = params.read_length ? params.read_length.toInteger() - 1 : 99

    """
    mkdir -p STAR_index

    echo "Running STAR genomeGenerate..."
    STAR --runMode genomeGenerate \\
         --genomeDir STAR_index \\
         --genomeFastaFiles ${genome_fasta} \\
         --sjdbGTFfile ${genome_gtf} \\
         --sjdbOverhang ${sjdbOverhang} \\
         --runThreadN ${task.cpus}

    STAR --version > STAR_index/STAR.version.txt
    """
}
