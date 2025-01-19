process DOWNLOAD_SNPEFF_DB {
    tag "Download SnpEff Database"
	publishDir "${params.test_data_dir}/snpEff", mode: 'copy'
	container null
    
	input:
    val genome
    path snpeff_jar_path

    output:
    path "${params.snpeff_db_dir}/${genome}"
	

    script:
    """
    # Ensure the output directory exists first
    mkdir -p ${params.snpeff_db_dir}

    # Use an absolute path for the data directory
    data_dir=\$(realpath ${params.snpeff_db_dir})

    # Download the database
    # Check if the database is already downloaded
    if [ ! -d "\$data_dir/${genome}" ]; then
        echo "Downloading SnpEff database for ${genome}..."
        java -Xmx4g -Xms2g -jar ${snpeff_jar_path} download ${genome} -dataDir \$data_dir -v
    else
        echo "SnpEff database for ${genome} already exists. Skipping download."
    fi
    """
}
