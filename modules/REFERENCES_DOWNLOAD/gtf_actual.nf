// Process to download gtf files
process DOWNLOAD_GTF {
    tag "Download GTF"
	container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "annotations.gtf"

	when:
    !file("${params.actual_data_dir}/reference/annotations.gtf").exists()
	
	
    script:
    """
     # Download the GTF file
    wget -q -O annotations.gtf.gz ${params.test_data_gtf}

    # Check if the file is gzipped and uncompress if necessary
    if file annotations.gtf.gz | grep -q 'gzip'; then
        gunzip annotations.gtf.gz
    else
        mv annotations.gtf.gz annotations.gtf
    fi
    """
}