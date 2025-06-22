process CHECK_OR_DOWNLOAD_DENYLIST {
    tag "Denylist BED"
    label 'process_light'
    publishDir "${params.ref_base}/reference", mode: 'copy'
    container null  // or specify if needed

    output:
    path "denylist.bed", emit: denylist

    script:
    """
    echo "Downloading denylist BED file..."
    wget -q -O denylist.bed ${params.denylist_download_url}
    
    echo "First few lines of downloaded denylist:"
    head denylist.bed
    """
}
