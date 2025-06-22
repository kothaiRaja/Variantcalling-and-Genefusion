process DOWNLOAD_GTF {
    tag "Download GTF"
    publishDir "${params.ref_base}/reference", mode: 'copy'
    container null
    label 'process_high'

    output: 
    path "annotations.gtf.gz", emit: gtf_gz

    when:
    !params.gtf_annotation && !file("${params.ref_base}/reference/annotations.gtf.gz").exists()

    script:
    """
    echo "Downloading GTF annotation from: ${params.gtf_download_url}"
    wget -q -O annotations.gtf.gz ${params.gtf_download_url}
    """
}
