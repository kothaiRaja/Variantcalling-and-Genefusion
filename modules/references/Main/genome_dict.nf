process CREATE_GENOME_DICT {
    tag "Create Genome Dictionary"
    label 'process_medium'

    container "https://depot.galaxyproject.org/singularity/picard%3A2.27.4--hdfd78af_0"
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    input:
    path genome_fa

    output:
    path "*.dict", emit: genome_dict

    script:
    """
    echo "Creating genome dictionary using Picard..."

    BASENAME=\$(basename "${genome_fa}" .fa)
    picard CreateSequenceDictionary R=${genome_fa} O=\${BASENAME}.dict
    """
}
