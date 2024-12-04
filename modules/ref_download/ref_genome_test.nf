// Process to download the test genome
process DOWNLOAD_TEST_GENOME {
    tag "Download test genome"
	container null
    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "genome.fa"

    script:
    """
    wget -q -O genome.fa ${params.test_data_url}/genome.fa
    """
}