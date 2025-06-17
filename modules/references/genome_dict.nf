process CREATE_GENOME_DICT {
    tag "${task.process}"

    publishDir "${params.ref_base}/reference", mode: 'copy'
	container params.gatk_container
    label 'process_medium'

    input:
    path genome_fasta

    output:
    path("${genome_fasta.simpleName}.dict"), emit: genome_dict

    script:
    """
    echo "Creating genome dictionary using GATK (Picard)..."
    gatk CreateSequenceDictionary -R ${genome_fasta} -O ${genome_fasta.simpleName}.dict
    """
}
