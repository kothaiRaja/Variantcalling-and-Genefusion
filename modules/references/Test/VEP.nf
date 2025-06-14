process DOWNLOAD_VEP_CACHE {
    tag "Download VEP Cache ${params.species} Ensembl ${params.ensembl_release} ${params.genome_assembly}"
	label 'process_high'
    publishDir "${params.test_data_dir}/Tools/VEP", mode: 'copy'
    container null

    output:
    path "vep_cache", emit: vep_cache  

    script:
    """
    mkdir -p vep_cache  
    echo "Downloading VEP cache for ${params.species}, Ensembl ${params.ensembl_release}, ${params.genome_assembly}..."

    wget -q -O vep_cache.tar.gz https://ftp.ensembl.org/pub/release-${params.ensembl_release}/variation/indexed_vep_cache/${params.species}_vep_${params.ensembl_release}_${params.genome_assembly}.tar.gz


    tar -xzvf vep_cache.tar.gz -C vep_cache --strip-components=1
    rm vep_cache.tar.gz
    """
}

process DOWNLOAD_VEP_PLUGINS {
    tag "Download VEP Plugins"
    label 'process_high'
    publishDir "${params.test_data_dir}/Tools/VEP/plugins", mode: 'copy'
    container null

    output:
    path "plugins", emit: vep_plugins

   

    script:
    """
    mkdir -p plugins
	cd plugins
    echo "Cloning VEP plugins from Ensembl GitHub..."
    
    git clone https://github.com/Ensembl/VEP_plugins.git 

    
    """
}
