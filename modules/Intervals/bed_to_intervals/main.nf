process BED_TO_INTERVAL_LIST {
    
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/intervals", mode: "copy"

    input:
    tuple val(meta_id), path(bed_file)
	path(genome_fasta)
	path(genome_dict)

    output:
    tuple val(meta_id), path("${bed_file.baseName}.interval_list")

    script:
    """
    gatk BedToIntervalList \
        -I ${bed_file} \
        -O ${bed_file.baseName}.interval_list \
        -SD ${genome_dict}
    """
}