process MERGE_BAMS {
    tag "MERGE_BAMS"
    label 'process_medium'

    container params.samtools_container
    publishDir params.merge_bam_outdir, mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam_list), path(bai_list)

    output:
    tuple val(sample_id), val(strandedness), 
          path("${sample_id}_merged.bam"), 
          path("${sample_id}_merged.bam.bai"), emit: merged_bams
    path("versions.yml"), emit: versions

    script:
    """
    echo "Merging BAM files for sample: ${sample_id}"

    # List BAM files into a text file
    ls "${bam_list}" > bam_files.txt

    # Check if any BAM files are listed
    if [ ! -s bam_files.txt ]; then
        echo "ERROR: No BAM files found for sample: ${sample_id}"
        exit 1
    fi

    # Merge BAMs
    samtools merge -@ ${task.cpus} -b bam_files.txt -o "${sample_id}_merged.bam"

    # Index the merged BAM
    samtools index "${sample_id}_merged.bam"

    # Capture Samtools version
    samtools_version=\$(samtools --version | head -n 1 | awk '{print \$2}')
    cat <<EOF > versions.yml
    "${task.process}":
      samtools: "\${samtools_version}"
    EOF

    echo "Merge complete for sample: ${sample_id}"
    """
}
