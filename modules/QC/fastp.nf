process fastp {
    tag "${sample_id}"
    container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--hadf994f_3"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.html"), path("${sample_id}_fastp.json"), emit: fastp_reports

    script:
    """
    fastp \
        -i ${reads[0]} -I ${reads[1]} \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -O ${sample_id}_R2_trimmed.fastq.gz \
        --thread ${task.cpus} \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json
    """
}
