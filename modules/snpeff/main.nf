process ANNOTATE_VARIANTS {
    tag "Annotate variants"

    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path(snpEffJar)
    path(snpEffConfig)
    path(snpEffDbDir)
    val(genomedb)

    output:
    tuple val(sample_id), path("annotated_${sample_id}.vcf"), emit: annotated_vcf
    path "annotated_${sample_id}.summary.html", emit: summary_html

    script:
    """
    THREADS=${task.cpus}

    # Annotate variants and generate annotated VCF file
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        ${vcf} > annotated_${sample_id}.vcf

    # Generate annotation summary report
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated_${sample_id}.summary.html \
        ${vcf} > /dev/null
    """
}