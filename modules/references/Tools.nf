process DOWNLOAD_SNPEFF_TOOL {
    tag "Download SnpEff Tool"
    publishDir "${params.ref_base}/Tools/snpEff", mode: 'copy'
    container null
    label 'process_medium'

    output:
    path "snpEff.jar", emit: snpeff_jar
    path "snpEff.config", emit: snpeff_config
 

    script:
    """
    
   
    echo "Downloading SnpEff..."
    wget -q -O snpEff.zip https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download

    echo "Extracting snpEff.jar and snpEff.config..."
    unzip -o -j snpEff.zip 'snpEff/snpEff.jar' 'snpEff/snpEff.config' -d .


    echo "Cleanup..."
    rm snpEff.zip

    echo "SnpEff tool downloaded and extracted to \$(pwd)"
    """
}



process DOWNLOAD_SNPEFF_DB {
    tag "Download SnpEff Database"
    publishDir "${params.ref_base}/Tools/snpEff/data", mode: 'copy'
    container null
    label 'process_medium'

    input:
    val genome
    path snpeff_jar
    path snpeff_config

    output:
    path "${genome}", emit: snpeff_db_dir

    script:
    """
    echo "Preparing SnpEff database download for ${genome}..."

    mkdir -p data
    cp ${snpeff_config} data/snpeff.config
    config_file="data/snpeff.config"

    # Patch the config
    if ! grep -q "^${genome}.genome" \$config_file; then
        echo "${genome}.genome : Homo_sapiens" >> \$config_file
    fi

    if ! grep -q "^database.repository" \$config_file; then
        echo "database.repository : https://snpeff.blob.core.windows.net/databases/" >> \$config_file
    fi

    # Download database
    java -Xmx4g -Xms2g -jar ${snpeff_jar} download ${genome} -dataDir ./ -c \$config_file -v
    """
}



process DOWNLOAD_ARRIBA {
    tag "Download Arriba v${params.arriba_version}"
    container null
    publishDir "${params.ref_base}/Tools/ARRIBA", mode: 'copy'

    output:
    path "arriba_v${params.arriba_version}", emit: arriba_dir
    

    script:
    """
    set -euo pipefail

    URL="https://github.com/suhrig/arriba/releases/download/v${params.arriba_version}/arriba_v${params.arriba_version}.tar.gz"
    TARGET_DIR="arriba_v${params.arriba_version}"

    echo "Downloading Arriba from \$URL..."
    wget -q -O arriba.tar.gz "\$URL"

    echo "Extracting Arriba..."
    mkdir -p "\$TARGET_DIR"
    tar -xzf arriba.tar.gz -C "\$TARGET_DIR" --strip-components=1

    echo "Cleaning up..."
    rm arriba.tar.gz

    echo "Arriba v${params.arriba_version} setup completed."
    """
}

