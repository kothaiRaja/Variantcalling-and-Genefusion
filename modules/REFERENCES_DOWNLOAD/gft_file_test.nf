// Process to download gtf files
process DOWNLOAD_GTF {
    tag "Download test GTF"
	container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "annotations.gtf"

    script:
    """
    wget -q -O annotations.gtf ${params.test_data_gtf}
    """
}