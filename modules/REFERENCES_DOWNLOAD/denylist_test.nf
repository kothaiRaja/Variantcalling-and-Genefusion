// Process to download the test denylist BED
process DOWNLOAD_TEST_DENYLIST {
    tag "Download test denylist BED"
	container null	
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "denylist.bed"

    script:
    """
    wget -q -O denylist.bed ${params.test_data_denylist}
	
	
    """
}