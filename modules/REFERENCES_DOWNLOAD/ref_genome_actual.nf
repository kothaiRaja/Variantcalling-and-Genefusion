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
    # Download the genome file
	wget -q -O genome.fa.gz ${params.actual_data_genome}

	# Check if the file is gzip-compressed
	if file genome.fa.gz | grep -q 'gzip'; then
		gunzip genome.fa.gz
	else
		mv genome.fa.gz genome.fa
	fi

    """
}