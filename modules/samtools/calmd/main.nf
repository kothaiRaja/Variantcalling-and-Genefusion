process SAMTOOLS_CALMD {
    tag { "${sample_id}_${bam.baseName}" }  

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"

    publishDir "${params.outdir}/calmd", mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai) 
    path genome_fasta
    path index

    output:
    tuple val(sample_id),val(strandedness), 
          path("${bam.baseName}_calmd.bam"),
          path("${bam.baseName}_calmd.bam.bai")
          

    script:
    """
    echo "Processing BAM file: ${bam} with samtools calmd"

    # Add NM and MD tags using samtools calmd
    samtools calmd -b ${bam} ${genome_fasta} > ${bam.baseName}_calmd.bam

    # Verify BAM file exists before proceeding
    if [ ! -s ${bam.baseName}_calmd.bam ]; then
        echo "Error: samtools calmd did not produce a valid BAM file!" >&2
        exit 1
    fi

    # Ensure BAM is sorted before indexing
    samtools sort -o ${bam.baseName}_calmd.sorted.bam ${bam.baseName}_calmd.bam
    mv ${bam.baseName}_calmd.sorted.bam ${bam.baseName}_calmd.bam

    # Index the BAM file
    samtools index ${bam.baseName}_calmd.bam

    # Verify BAI file exists
    if [ ! -s ${bam.baseName}_calmd.bam.bai ]; then
        echo "Error: samtools index did not generate a BAI file!" >&2
        exit 1
    fi

    echo "Finished processing ${bam.baseName}_calmd.bam"
    """
}
