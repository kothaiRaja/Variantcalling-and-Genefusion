nextflow.enable.dsl = 2

process fastqc_raw {
    tag "${sample_id}"
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    publishDir "${params.outdir}/fastqc/raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.zip"), emit: fastqc_zip
    path("*.html"), emit: fastqc_html

    script:
    """
    fastqc --outdir . ${reads.join(' ')}
    """
}

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

process fastqc_trimmed {
    tag "${sample_id}"
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    publishDir "${params.outdir}/fastqc/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.zip"), path("*.html")

    script:
    """
    fastqc --outdir . ${reads.join(' ')}
    """
}

workflow qc {
    reads_channel = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, reads -> tuple(sample_id, reads.sort()) }

    fastqc_raw_results = fastqc_raw(reads_channel)
    fastp_results = fastp(reads_channel)
    trimmed_reads = fastp_results.trimmed_reads
        .map { sample_id, r1, r2 -> tuple(sample_id, [r1, r2]) }
    fastqc_trimmed_results = fastqc_trimmed(trimmed_reads)

    return [
        fastqc_raw_results: fastqc_raw_results,
        fastqc_trimmed_results: fastqc_trimmed_results,
        trimmed_reads: trimmed_reads
    ]
}
