process DOWNLOAD_VEP_CACHE {
    tag "Downloading VEP Cache"
	container null
	publishDir "${params.test_data_dir}/VEP", mode: 'copy'
    output:
    path "vep_cache"

    script:
    """
    mkdir -p vep_cache
    wget -O homo_sapiens_vep_110_GRCh38.tar.gz https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
    
    # Check if the download was successful
    if [ ! -s homo_sapiens_vep_110_GRCh38.tar.gz ]; then
        echo "Error: Failed to download VEP cache file" >&2
        exit 1
    fi

    tar -xvzf homo_sapiens_vep_110_GRCh38.tar.gz -C vep_cache
    rm homo_sapiens_vep_110_GRCh38.tar.gz
    """
}
