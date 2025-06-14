// ========================== Download Reference Genome and Index ========================== //
process CHECK_OR_DOWNLOAD_REF_GENOME {
    tag "Check or Download Reference Genome"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa", emit: genome

    when:
    !file("${params.test_data_dir}/reference/genome.fa").exists()

    script:
    """
    echo "Downloading reference genome from: ${params.genome_download_url}"
    wget -q -O genome.fa.gz ${params.genome_download_url}

    echo "Unzipping reference genome..."
    if file genome.fa.gz | grep -q 'gzip'; then
        gunzip genome.fa.gz
    else
        mv genome.fa.gz genome.fa
    fi

    echo "Checking Chromosome Naming Format in genome.fa..."
    FIRST_CHR=\$(grep '^>' genome.fa | head -1 | sed 's/^>//' | awk '{print \$1}')
    echo "  → Detected contig name: '\$FIRST_CHR'"
    """
}
