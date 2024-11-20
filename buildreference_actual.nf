nextflow.enable.dsl = 2

// Process to download a single file
process DOWNLOAD_FILE {
    tag "$name" // Tags the process with the file name for better logging
	publishDir "${params.actual_data_dir}", mode: 'copy'
    input:
    tuple val(name), val(url) // Input is a tuple: file name and URL

    output:
    path "$name" publishDir "${params.actual_data_dir}" // Publishes to actual data directory

    script:
    """
    wget -q -O $name $url
    """
}

// Process to extract reads from a tar.gz archive
process EXTRACT_READS {
    tag "Extract reads" // Tags the process for better logging
	publishDir "${params.actual_data_dir}", mode: 'copy'
    input:
    path "reads.tar.gz" // Takes the reads archive as input

    output:
    path "reads" publishDir "${params.actual_data_dir}/reads" // Publishes the extracted reads

    script:
    """
    mkdir -p reads
    tar -xzf reads.tar.gz -C reads
    rm reads.tar.gz
    """
}

// Main workflow
workflow {
    // Step 1: Channel for files to download
    files_to_download = Channel.from(params.actual_urls)

    // Step 2: Download all files
    downloaded_files = files_to_download
        .map { name, url -> tuple(name, url) } // Map each name-URL pair into a tuple
        | DOWNLOAD_FILE

    // Step 3: Extract reads if needed
    extracted_reads = downloaded_files
        .filter { it.name == "reads.tar.gz" } // Filter to find the reads tar.gz file
        | EXTRACT_READS

    // Optional: Emit for debugging or logging
    emit:
    downloaded_files, extracted_reads
}
