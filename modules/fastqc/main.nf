process FASTQC_RAW {
    tag { "${meta.id}_${task.process}" }
    label 'process_low'

    publishDir params.fastqc_outdir, mode: "copy", pattern: "*_fastqc.*"
    container params.fastqc_container

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*_fastqc.zip"), emit: qc_results
    path("*_fastqc.html"), emit: fastqc_reports
    path("versions.yml"), emit: versions

    script:
    """
    fastqc ${r1} ${r2} --outdir .

# Capture FastQC version
fastqc_version=\$(fastqc --version | awk '{print \$2}')
cat <<EOF > versions.yml
"${task.process}":
 fastqc: "\${fastqc_version}"
EOF
    """
}
