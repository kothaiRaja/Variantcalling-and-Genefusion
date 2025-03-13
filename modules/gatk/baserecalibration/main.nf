process GATK_BASERECALIBRATOR {
    tag { sample_id }
	
    cpus params.get('gatk_recalibration_cpus', 8)
    memory params.get('gatk_recalibration_memory', '32 GB')
    time params.get('gatk_recalibration_time', '5h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/recal_tables", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai), val(strandedness), path(interval)
    path(genome_fasta)
    path(index)
    path(dict)
    path(known_variants)
    path(known_variants_index)

    output:
    tuple val(sample_id), path("${sample_id}_recal_data.table"), val(strandedness)

    script:

    def interval_command = interval ? "--intervals ${interval}" : ""

    """
    THREADS=${task.cpus}

    # Step 1: BaseRecalibrator
    gatk BaseRecalibrator \
        -R ${genome_fasta} \
        -I ${bam} \
        --known-sites ${known_variants} \
        -O ${sample_id}_recal_data.table \
		${interval_command} 
		
       

    # Check if recalibration table is generated
    if [ ! -s ${sample_id}_recal_data.table ]; then
        echo "Error: Recalibration table not generated for ${sample_id}" >&2
        exit 1
    fi
    """
}
