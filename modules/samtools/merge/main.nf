process MERGE_BAMS {
    tag "MERGE_BAMS"

    cpus params.get('merge_bam_cpus', 4)
    memory params.get('merge_bam_memory', '8 GB')
    time params.get('merge_bam_time', '2h')

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.15.1--h1170115_0"

    input:
    tuple val(sample_id), path(bam_list), path(bai_list), val(strandedness)

    output:
    tuple val(sample_id), 
          path("${sample_id}_merged.bam"), 
          path("${sample_id}_merged.bam.bai"), 
          val(strandedness)

    script:
    """
    echo "Merging BAM files for sample: ${sample_id}"

    # Create a temporary file to store BAM paths
    ls ${bam_list} > bam_files.txt

    # Ensure the file contains BAM files
    if [ ! -s bam_files.txt ]; then
        echo "ERROR: No BAM files found for sample: ${sample_id}"
        exit 1
    fi

    # Merge BAM files
    samtools merge -@ ${task.cpus} -b bam_files.txt -o ${sample_id}_merged.bam

    # Index BAM file
    samtools index ${sample_id}_merged.bam

    echo "Merge complete for sample: ${sample_id}"
    """
}