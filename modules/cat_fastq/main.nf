process CONCAT_FASTQ {
    tag { "${meta.id}_${task.process}" }

    container "https://depot.galaxyproject.org/singularity/ubuntu%3A24.04"
    publishDir "${params.outdir}/concatenated", mode: "copy"

    input:
    tuple val(meta), path(fastq1), path(fastq2)

    output:
    tuple val(meta), path("${meta.id}_R1.merged.fastq.gz"), path("${meta.id}_R2.merged.fastq.gz"), emit: merged_reads

    script:
    """
    echo "Processing reads for sample: ${meta.id}"
    echo "FASTQ1: ${fastq1}"
    echo "FASTQ2: ${fastq2}"

    if [[ \$(zcat ${fastq1} ${fastq2} | wc -l) -gt 4 ]]; then
        echo "Merging (even if only one file each — just standardizing)..."
        cat ${fastq1} | gzip > ${meta.id}_R1.merged.fastq.gz
        cat ${fastq2} | gzip > ${meta.id}_R2.merged.fastq.gz
        echo "Merging complete for sample: ${meta.id}"
    else
        echo "Copying files directly for sample: ${meta.id}"
        cp ${fastq1} ${meta.id}_R1.merged.fastq.gz
        cp ${fastq2} ${meta.id}_R2.merged.fastq.gz
    fi
    """
}
