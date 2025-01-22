// Process to download the test genome
process DOWNLOAD_REF_GENOME {
    tag "Download test genome"
	container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa"
	
	when:
    !file("${params.test_data_dir}/reference/genome.fa").exists()

    script:
    """
   # Download the genome file
	wget -q -O genome.fa.gz ${params.test_data_genome}

	# Check if the file is gzip-compressed
	if file genome.fa.gz | grep -q 'gzip'; then
		gunzip genome.fa.gz
	else
		mv genome.fa.gz genome.fa
	fi
    """
}