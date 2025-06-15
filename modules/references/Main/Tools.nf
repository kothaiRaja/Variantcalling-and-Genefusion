process DOWNLOAD_SNPEFF_TOOL {
    tag "Download SnpEff Tool"
    publishDir "${params.main_data_dir}/Tools", mode: 'copy'
    container null  

    output:
    path "${params.snpeff_jar_dir}/snpEff.jar", emit: snpeff_jar
    path "${params.snpeff_jar_dir}/snpEff.config", emit: snpeff_config

    script:
    """
    mkdir -p ${params.snpeff_jar_dir}
    wget -q -O snpEff_latest_core.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip -o -j snpEff_latest_core.zip -d ${params.snpeff_jar_dir}
    rm snpEff_latest_core.zip
    """
}

process DOWNLOAD_SNPEFF_DB {
    tag "Download SnpEff Database"
    publishDir "${params.main_data_dir}/Tools/snpEff", mode: 'copy'
    container null

    input:
    val genome
    path snpeff_jar_path

    output:
    path "${params.snpeff_db_dir}/${genome}"

    script:
"""
mkdir -p ${params.snpeff_db_dir}
data_dir=\$(realpath ${params.snpeff_db_dir})

# Correct config path based on jar location
config_file=\$(dirname ${snpeff_jar_path})/snpEff.config

# Patch config if missing genome entry
if ! grep -q "^${genome}.genome" \$config_file; then
    echo "${genome}.genome : Homo_sapiens" >> \$config_file
fi

# Also ensure repository line exists
if ! grep -q "^database.repository" \$config_file; then
    echo "database.repository : https://snpeff.blob.core.windows.net/databases/" >> \$config_file
fi

# Download the database if not present
if [ ! -d "\$data_dir/${genome}" ]; then
    echo "Downloading SnpEff database for ${genome}..."
    java -Xmx4g -Xms2g -jar ${snpeff_jar_path} download ${genome} -dataDir \$data_dir -c \$config_file -v
else
    echo "SnpEff database for ${genome} already exists. Skipping download."
fi
"""

    
}



process DOWNLOAD_ARRIBA {
    tag "Download Arriba v${params.arriba_version}"
    container null
    publishDir "${params.main_data_dir}/Tools/ARRIBA", mode: 'copy'

    output:
    path "arriba_v${params.arriba_version}", emit: arriba_dir

    when:
    !file("${params.actual_data_dir}/Tools/ARRIBA/arriba_v${params.arriba_version}").exists()

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

