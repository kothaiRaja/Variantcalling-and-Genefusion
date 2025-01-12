process DOWNLOAD_ARRIBA {
    container null
	publishDir "${params.actual_data_dir}/ARRIBA", mode: 'copy' 
    
	output:
    path 'arriba_v2.4.0', emit: 'arriba_dir'

    script:
    """
    # Define the Arriba download URL and target directory
    URL="https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz"
    TARGET_DIR="arriba_v2.4.0"

    # Download Arriba tarball
    wget -O arriba.tar.gz \$URL

    # Create the target directory and extract the tarball
    mkdir -p \$TARGET_DIR
    tar -xzvf arriba.tar.gz -C \$TARGET_DIR --strip-components=1

    # Clean up the downloaded tarball
    rm arriba.tar.gz

    # Output the extracted directory
    echo "Arriba tool extracted to: \$TARGET_DIR"
    """
}