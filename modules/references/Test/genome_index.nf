process DOWNLOAD_GENOME_INDEX {
    tag "Download Genome Index"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa.fai", emit: genome_fai

    script:
    """
    echo "Downloading genome index from provided URL..."
    wget -q -O genome.fa.fai ${params.genome_index_download_url}
    
    echo "Checking Chromosome Naming Format in genome index..."
    FIRST_CHR=\$(awk '{print \$1; exit}' genome.fa.fai)
    
    if [[ "\$FIRST_CHR" == chr* ]]; then
        echo "Detected 'chr' prefix in genome index. Converting to numeric format..."
        sed -i 's/^chr//' genome.fa.fai
        echo "Genome index chromosome names converted to numeric format."
    else
        echo "Genome index chromosome names are already in numeric format. No changes needed."
    fi
    """
}



process CREATE_GENOME_INDEX {
    tag "Create Genome Index"
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.14--hb421002_0"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    input:
    path genome_fa 

    output:
    path("genome.fa.fai"), emit: genome_fai  


    script:
    """
    echo "Creating genome index using samtools..."
    samtools faidx $genome_fa
    """
}
