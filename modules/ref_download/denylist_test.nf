// Process to download the test denylist BED
process DOWNLOAD_TEST_DENYLIST {
    tag "Download test denylist BED"
	container null	
    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "denylist.bed"

    script:
    """
    wget -q -O denylist.bed ${params.test_data_url}/denylist_chr22_to_22.bed
    """
}