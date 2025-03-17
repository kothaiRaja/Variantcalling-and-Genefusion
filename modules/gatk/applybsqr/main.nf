process GATK_APPLYBQSR {
    tag { sample_id }

    cpus params.get('gatk_applybqsr_cpus', 8)
    memory params.get('gatk_applybqsr_memory', '32 GB')
    time params.get('gatk_applybqsr_time', '5h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/recalibrated_bams", mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai), path(recal_table), path(interval)
    path(genome_fasta)
    path(index)
    path(dict)

    output:
    tuple val(sample_id),val(strandedness),	
          path("${sample_id}_${interval.baseName}_recalibrated.bam"),
          path("${sample_id}_${interval.baseName}_recalibrated.bai")
          

    script:
    
    def interval_command = interval ? "--intervals ${interval}" : ""

    """
    THREADS=${task.cpus}

    # Step 2: ApplyBQSR
    gatk ApplyBQSR \\
        -R ${genome_fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${sample_id}_${interval.baseName}_recalibrated.bam \\
        $interval_command

    # Check if recalibrated BAM is generated
    if [ ! -s ${sample_id}_${interval.baseName}_recalibrated.bam ]; then
        echo "Error: Recalibrated BAM not generated for ${sample_id}" >&2
        exit 1
    fi

    # Step 3: Index BAM using GATK
    gatk BuildBamIndex \\
        -I ${sample_id}_${interval.baseName}_recalibrated.bam
    """
}