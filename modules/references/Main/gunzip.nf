process GUNZIP {
    tag "$archive"
    container "https://depot.galaxyproject.org/singularity/ubuntu:24.04"
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    input:
    path(archive)

    output:
    path(decompressed_file), emit: gunzip

    when:
    archive.toString().endsWith('.gz')

    script:
    def decompressed_file = archive.getBaseName()  // Remove .gz but keep filename

    """
    echo "Decompressing file: $archive"
    gunzip -c $archive > ${decompressed_file}
    
    echo "Decompressed file saved as: ${decompressed_file}"
    """
}