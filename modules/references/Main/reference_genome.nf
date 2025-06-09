process CHECK_OR_DOWNLOAD_REF_GENOME {
    tag "Check or Download Reference Genome"
    container null
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa", emit: genome

    when:
    !file("${params.main_data_dir}/reference/genome.fa").exists()

    script:
    """
    echo "Downloading reference genome..."
    wget -q -O genome.fa.gz ${params.genome_download_url}

    if file genome.fa.gz | grep -q 'gzip'; then
        gunzip genome.fa.gz
    else
        mv genome.fa.gz genome.fa
    fi

    echo "Checking chromosome naming format..."
    FIRST_CHR=\$(grep -m1 '^>' genome.fa | sed 's/^>//')

    if [[ "\$FIRST_CHR" =~ ^chr ]]; then
        echo "Detected 'chr' prefix. Converting to numeric format..."
        sed -i 's/^>chr/>/g' genome.fa
    else
        echo "Chromosome names are already in numeric format."
    fi

    md5sum genome.fa > genome.fa.md5
    """
}
