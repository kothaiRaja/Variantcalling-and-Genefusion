// Process to download the denylist BED
process DOWNLOAD_DENYLIST {
    tag "Download denylist BED"
	container null	
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "denylist.bed"

    script:
    """
    wget -q -O denylist.bed ${params.actual_data_denylist}
	
	
    """
}