process DOWNLOAD_VEP_CACHE {
    tag "VEP Cache: ${params.species} Ensembl ${params.ensembl_release} ${params.genome_assembly}"
    label 'process_high'
    publishDir "${params.main_data_dir}/Tools/VEP", mode: 'copy'
    container null

    output:
    path "vep_cache", emit: vep_cache

    when:
    !file("${params.main_data_dir}/Tools/VEP/vep_cache/${params.species}_vep_${params.ensembl_release}_${params.genome_assembly}")

    script:
    """
    mkdir -p vep_cache
    echo "Downloading VEP cache for ${params.species}, Ensembl ${params.ensembl_release}, Assembly ${params.genome_assembly}..."

    curl -sSLo vep_cache.tar.gz https://ftp.ensembl.org/pub/release-${params.ensembl_release}/variation/indexed_vep_cache/${params.species}_vep_${params.ensembl_release}_${params.genome_assembly}.tar.gz || {
        echo "ERROR: VEP cache download failed"
        exit 1
    }

    tar -xzvf vep_cache.tar.gz -C vep_cache || {
        echo "ERROR: Failed to extract VEP cache"
        exit 1
    }

    rm vep_cache.tar.gz
    """
}

process DOWNLOAD_VEP_PLUGINS {
    tag "Download VEP Plugins"
    label 'process_high'
    publishDir "${params.main_data_dir}/Tools/VEP/plugins", mode: 'copy'
    container null

    output:
    path "plugins", emit: vep_plugins

    when:
    !file("${params.main_data_dir}/Tools/VEP/plugins/Condel.pm") // or any key plugin

    script:
    """
    mkdir -p plugins
    echo "Cloning VEP plugins from Ensembl GitHub..."
    
    git clone https://github.com/Ensembl/VEP_plugins.git plugins

    echo "VEP plugins downloaded to: plugins/"
    """
}
