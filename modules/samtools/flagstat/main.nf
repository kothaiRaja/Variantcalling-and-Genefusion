process SAMTOOLS_FLAGSTAT {
    tag { "${meta.id}_${task.process}" }
    label 'process_medium'

    container params.samtools_container

    publishDir params.samtools_flagstat_outdir, mode: "copy", pattern: "*_flagstat.*"

    input:
    tuple val(meta), path(sorted_bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_flagstat.txt"), emit: flagstat
    path("versions.yml"), emit: versions

    script:
    def sample_id = meta.id

    """
    # Validate BAM and index file
    if [ ! -s ${sorted_bam} ]; then
        echo "Error: BAM file is empty or missing for ${sample_id}" >&2
        exit 1
    fi

    if [ ! -s ${bai} ]; then
        echo "Error: BAM index (.bai) file is missing for ${sample_id}" >&2
        exit 1
    fi

    # Validate BAM file to ensure it is not corrupted
    samtools quickcheck -v ${sorted_bam} || (echo "BAM file validation failed for ${sample_id}" && exit 1)

    # Run samtools flagstat to generate alignment metrics
    samtools flagstat -@ ${task.cpus} ${sorted_bam} > ${sample_id}_flagstat.txt

    # Generate additional quality metrics with samtools stats
    samtools stats -@ ${task.cpus} ${sorted_bam} > ${sample_id}_stats_report.txt

    # Capture samtools version
    samtools_version=\$(samtools --version | head -n 1 | awk '{print \$2}')
    cat <<EOF > versions.yml
"${task.process}":
  samtools: "\${samtools_version}"
EOF
    """
}
