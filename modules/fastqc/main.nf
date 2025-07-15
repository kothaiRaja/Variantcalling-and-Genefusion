process FASTQC_RAW {
    tag { "${meta.id}_${task.process}" }
    label 'process_low'

    publishDir params.fastqc_outdir, mode: "copy", pattern: "*_fastqc.*"
    container params.fastqc_container

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1_fastqc.zip"), path("${meta.id}_R1_fastqc.html"),
                  path("${meta.id}_R2_fastqc.zip"), path("${meta.id}_R2_fastqc.html"), emit: qc_results
    path("versions.yml"), emit: versions

    script:
    """
    # Run FastQC
    fastqc ${r1} ${r2} --outdir .

    # Extract original FastQC filenames
    r1_base=\$(basename ${r1} | sed 's/\\.gz//; s/\\.fastq//; s/\\.fq//')
    r2_base=\$(basename ${r2} | sed 's/\\.gz//; s/\\.fastq//; s/\\.fq//')

    # Rename outputs to use meta.id
    mv \${r1_base}_fastqc.zip  ${meta.id}_R1_fastqc.zip
    mv \${r1_base}_fastqc.html ${meta.id}_R1_fastqc.html

    mv \${r2_base}_fastqc.zip  ${meta.id}_R2_fastqc.zip
    mv \${r2_base}_fastqc.html ${meta.id}_R2_fastqc.html

    # Capture FastQC version
    fastqc_version=\$(fastqc --version | awk '{print \$2}')
cat <<EOF > versions.yml
"${task.process}":
  fastqc: "\${fastqc_version}"
EOF
    """
}
