process DOWNLOAD_GENOME_INDEX {
    tag "Download Genome Index"
    container null
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa.fai", emit: genome_fai

    script:
    """
    echo "Downloading genome index from provided URL..."
    wget -q -O genome.fa.fai ${params.genome_index_download_url}

    echo "Checking Chromosome Naming Format in genome.fa.fai..."
    FIRST_CHR=\$(awk '{print \$1; exit}' genome.fa.fai)
    echo "  → Detected contig name in index: '\$FIRST_CHR'"
    """
}



process CREATE_GENOME_INDEX {
    tag "Create Genome Index"
	label 'process_medium'
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.14--hb421002_0"
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    input:
    path genome_fa 

    output:
    path("genome.fa.fai"), emit: genome_fai  


    script:
    """
    echo "Creating genome index using samtools..."

    if [[ "${genome_fa}" == *.gz ]]; then
        echo "Input is gzipped. Unzipping..."
        gunzip -c ${genome_fa} > genome.fa
    else
        cp ${genome_fa} genome.fa
    fi

    samtools faidx genome.fa

    echo "Checking Chromosome Naming Format in genome.fa.fai..."
    FIRST_CHR=\$(awk '{print \$1; exit}' genome.fa.fai)
    echo "  → Detected contig name in index: '\$FIRST_CHR'"
    
    """
}
