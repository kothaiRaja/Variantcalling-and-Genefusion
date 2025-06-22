process CREATE_GENOME_DICT {
    tag "${task.process}"

    publishDir "${params.ref_base}/reference", mode: 'copy'
    container params.gatk_container
    label 'process_medium'

    input:
    path genome_fasta

    output:
    path("genome.dict"), emit: genome_dict

    script:
    """
    echo "Creating genome dictionary using GATK (Picard)..."

    # Create dict with GATK and rename to standard name
    gatk CreateSequenceDictionary -R ${genome_fasta} -O temp.dict
    mv temp.dict genome.dict
    """
}
