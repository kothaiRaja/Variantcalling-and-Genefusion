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
    wget -q -O genome.fa.gz ${params.genome_download_url}

    if file genome.fa.gz | grep -q 'gzip'; then
        gunzip genome.fa.gz
    else
        mv genome.fa.gz genome.fa
    fi

    echo "Checking Chromosome Naming Format in genome.fa..."

    # Extract first chromosome name
    FIRST_CHR=\$(grep '^>' genome.fa | head -1 | sed 's/>//')

    if [[ "\$FIRST_CHR" == chr* ]]; then
        echo "Detected 'chr' prefix. Converting to numeric format..."
        sed -i 's/>chr/>/g' genome.fa   
        sed -i 's/\\bchr//g' genome.fa  
        echo "Genome chromosome names converted to numeric format."
    else
        echo "Genome chromosome names are already in numeric format. No changes needed."
    fi
    """
}