process SAMTOOLS_SORT_INDEX {
    tag { sample_id }
    label 'process_medium'

    container params.samtools_container
    publishDir params.samtools_outdir, mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam_sorted
    path("versions.yml"), emit: versions


    script:
    """
    echo "Sorting and indexing BAM for sample: ${sample_id}"

    # Sort the BAM file
    samtools sort -o ${sample_id}_sorted.bam "${bam}"

    # Index the sorted BAM file
    samtools index ${sample_id}_sorted.bam

    # Capture the version
    samtools_version=\$(samtools --version | head -n 1 | awk '{print \$2}')
    cat <<EOF > versions.yml
    samtools:
      version: "\${samtools_version}"
    EOF
    """
}
