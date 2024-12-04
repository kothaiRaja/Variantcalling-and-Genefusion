// Process to download gtf files
process DOWNLOAD_TEST_GTF {
    tag "Download test GTF"
	container null
    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "annotations.gtf"

    script:
    """
    wget -q -O annotations.gtf ${params.test_data_url}/annotations.gtf
    """
}