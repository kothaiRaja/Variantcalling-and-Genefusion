// Process to download the denylist BED
process DOWNLOAD_DENYLIST {
    tag "Download denylist BED"
	container null	
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "denylist.bed"
	
	when:
    !file("${params.actual_data_dir}/reference/denylist.bed").exists()

    script:
    """
    wget -q -O denylist.bed ${params.actual_data_denylist}
	
	# Check if the file is in .gz format and decompress
    if [[ -f denylist.bed.gz ]]; then
        gunzip denylist.bed.gz
    fi
	
	"""
}