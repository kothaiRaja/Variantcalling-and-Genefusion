process DOWNLOAD_GENOME_INDEX {
    tag "Download Genome Index"
    container null
    publishDir "${params.ref_base}/reference", mode: 'copy'

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
    publishDir "${params.ref_base}/reference", mode: 'copy'

    input:
    path genome_fa 

    output:
    path("genome.fa.fai"), emit: genome_fai  


    script:
    """
    echo "Creating genome index using samtools..."


if [[ "genome.fa" == *.gz ]]; then
    echo "Input is gzipped. Unzipping..."
    gunzip -c genome.fa > genome.fa
else
    echo "Input is already uncompressed. Copying..."
    cp genome.fa genome.fa.copy
    mv genome.fa.copy genome.fa
fi

# Index using the exact expected filename
samtools faidx genome.fa



    
    """
}
