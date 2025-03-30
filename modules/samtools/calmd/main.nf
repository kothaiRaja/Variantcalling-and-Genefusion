process SAMTOOLS_CALMD {
    tag { "${sample_id}_${task.process}" }
 
    label 'process_medium'

    container params.samtools_container

    publishDir params.merged_calmd_outdir, mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai) 
    path genome_fasta
    path index

    output:
    tuple val(sample_id), val(strandedness), 
          path("${bam.baseName}_calmd.bam"),
          path("${bam.baseName}_calmd.bam.bai"), emit: calmd_bams
    path("versions.yml"), emit: versions

    script:
    """
    echo "Processing BAM file: ${bam} with samtools calmd"

    # Add NM and MD tags using samtools calmd
    samtools calmd -b "${bam}" "${genome_fasta}" > "${bam.baseName}_calmd.bam"

    # Verify BAM file exists before proceeding
    if [ ! -s "${bam.baseName}_calmd.bam" ]; then
        echo "Error: samtools calmd did not produce a valid BAM file!" >&2
        exit 1
    fi

    # Sort before indexing
    samtools sort -o "${bam.baseName}_calmd.sorted.bam" "${bam.baseName}_calmd.bam"
    mv "${bam.baseName}_calmd.sorted.bam" "${bam.baseName}_calmd.bam"

    # Index the BAM file
    samtools index "${bam.baseName}_calmd.bam"

    # Verify BAI file exists
    if [ ! -s "${bam.baseName}_calmd.bam.bai" ]; then
        echo "Error: samtools index did not generate a BAI file!" >&2
        exit 1
    fi

    # Capture samtools version
    samtools_version=\$(samtools --version | head -n 1 | awk '{print \$2}')
cat <<EOF > versions.yml
"${task.process}":
  samtools: "\${samtools_version}"
EOF

    echo "Finished processing ${bam.baseName}_calmd.bam"
    """
}
