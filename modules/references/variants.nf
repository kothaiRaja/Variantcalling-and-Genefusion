process CHECK_OR_DOWNLOAD_ALL_VARIANTS {
    tag "Download 00-All.vcf.gz"
    label 'process_light'
    publishDir "${params.ref_base}/reference", mode: 'copy'
    container null

    output:
    path "00-All.vcf.gz", emit: all_variants_vcf

    script:
    """
    echo "Downloading 00-All.vcf.gz from: ${params.variants_url}"
    wget -q -O 00-All.vcf.gz ${params.variants_url}

    echo "Download complete. Previewing header:"
    zgrep '^##' 00-All.vcf.gz | head
    """
}


