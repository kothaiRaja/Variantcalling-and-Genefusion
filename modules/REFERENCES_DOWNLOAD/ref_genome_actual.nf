// Process to download the reference genome
process DOWNLOAD_REF_GENOME {
    tag "Download reference genome"
	container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa"
	
	when:
    !file("${params.actual_data_dir}/reference/genome.fa").exists()


    script:
    """
    wget -q -O genome.fa ${params.actual_data_genome}
	
	 # Download the genome file
    wget -q -O genome.fa.gz ${params.actual_data_genome}

    # Rename or delete any existing genome.fa to avoid overwriting
    if [ -e genome.fa ]; then
        rm genome.fa
    fi

    # Check if the file is gzip-compressed and decompress it
    if file genome.fa.gz | grep -q 'gzip'; then
        gunzip genome.fa.gz
    fi
    """
}