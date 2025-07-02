process TRIM_READS {
    tag { "${meta.id}_${task.process}" }
    label 'process_medium'

    container params.trim_reads_container
    publishDir params.trim_reads_outdir, mode: 'copy' 

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("trimmed_${meta.id}_R1.fastq.gz"), path("trimmed_${meta.id}_R2.fastq.gz"), emit: trimmed_reads
    tuple path("${meta.id}_fastp.html"), path("${meta.id}_fastp.json"), emit: fastp_reports
    path("versions.yml"), emit: versions

    script:
    """
    fastp -i "${r1}" -I "${r2}" \\
        -o "trimmed_${meta.id}_R1.fastq.gz" \\
        -O "trimmed_${meta.id}_R2.fastq.gz" \\
        ${params.fastp_extra} \\
        --html "${meta.id}_fastp.html" \\
        --json "${meta.id}_fastp.json"

# Capture version
fastp_version=\$(fastp --version 2>&1 | head -n1 | awk '{print \$2}')

cat <<EOF > versions.yml
"${task.process}":
 fastp: "\${fastp_version}"
EOF
    """
}
