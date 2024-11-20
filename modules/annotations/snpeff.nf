process ANNOTATE_VARIANTS {
    tag "$sampleId"
    label "mem_large"
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
        tuple val(sampleId), path(vcf)

    output:
        tuple val(sampleId), path("${sampleId}.annotated.vcf"), path("${sampleId}.summary.html")

    script:
    """
    java -Xmx16G -jar ${params.snpeff_jar_dir}/snpEff.jar \\
        -c ${params.snpeff_jar_dir}/snpEff.config \\
        -v ${params.genomedb} \\
        -dataDir ${params.snpeff_db_dir} \\
        ${vcf} > ${sampleId}.annotated.vcf

    java -Xmx16G -jar ${params.snpeff_jar_dir}/snpEff.jar \\
        -c ${params.snpeff_jar_dir}/snpEff.config \\
        -v ${params.genomedb} \\
        -dataDir ${params.snpeff_db_dir} \\
        -stats ${sampleId}.summary.html \\
        ${vcf} > /dev/null
    """
}

