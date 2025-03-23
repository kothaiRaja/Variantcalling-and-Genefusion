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

# Merge BAMs directly using list expansion
samtools merge -@ ${task.cpus} -o "${sample_id}_merged.bam" ${bam_list.join(' ')}

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
