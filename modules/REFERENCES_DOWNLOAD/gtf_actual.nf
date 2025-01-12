// Process to download gtf files
process DOWNLOAD_GTF {
    tag "Download GTF"
	container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "annotations.gtf"

    script:
    """
    wget -q -O annotations.gtf ${params.actual_data_gtf}
    """
}