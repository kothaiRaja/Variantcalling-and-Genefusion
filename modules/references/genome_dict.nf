process CREATE_GENOME_DICT {
    tag "${task.process}"

    publishDir "${params.ref_base}/reference", mode: 'copy'
    container params.gatk_container
    label 'process_medium'

    input:
    path genome_fasta

    output:
    path("*.dict"), emit: genome_dict

    script:
    """
    echo "Creating genome dictionary using GATK (Picard)..."

    # Get the basename of the input FASTA
    fasta_basename=\$(basename "${genome_fasta}")

    # Strip .fa or .fasta to create the .dict name
    dict_output=\$(echo "\$fasta_basename" | sed -E 's/\\.fa(sta)?\$/.dict/')

    # Run GATK CreateSequenceDictionary
    gatk CreateSequenceDictionary -R "${genome_fasta}" -O "\$dict_output"
    """
}
