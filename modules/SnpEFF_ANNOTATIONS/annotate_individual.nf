process ANNOTATE_INDIVIDUAL_VARIANTS {
    tag "${sample_id}_annotate"

    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_index)
	path (snpEffJar)
	path (snpEffConfig)
	path (snpEffDbDir)
	val (genomedb)

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf"), path("${sample_id}.annotated.summary.html")

    script:
    """
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        ${filtered_vcf} > ${sample_id}.annotated.vcf
		



    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats ${sample_id}.annotated.summary.html \
        ${filtered_vcf} > /dev/null
    """
}
