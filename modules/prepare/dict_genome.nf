
//Create genome dictionary with PICARD for GATK4

process PREPARE_GENOME_PICARD {
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/genome/dict", mode: "copy"

    input:
    path genome

    output:
    path "${genome.baseName}.dict"

    script:
    """
    gatk CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
   
	# Correct any UR:file path in the resulting genome.dict file to absolute path
    sed -i 's|UR:file:.*|UR:file:/home/kothai/cq-git-sample/Praktikum/data/genome.fa|' ${genome.baseName}.dict
	"""
}

