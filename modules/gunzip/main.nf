process GUNZIP {
    tag "$input_file"
    container "https://depot.galaxyproject.org/singularity/ubuntu:24.04"
    publishDir "${params.ref_base}/reference", mode: 'copy'

    input:
    path input_file

    output:
    path("${input_file.baseName}"), emit: gunzip

    script:
    """
    gunzip -c ${input_file} > ${input_file.baseName}
    """
}
