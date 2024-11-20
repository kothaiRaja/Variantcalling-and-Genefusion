process ANNOTATE_VARIANTS {
    tag "$sampleId"
    label "mem_large"
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
        tuple val(sampleId), path(vcf)
        path("/home/kothai/cq-git-sample/Variantcalling-and-Genefusion/data/test/snpEff/snpEff.jar")
        path("/home/kothai/cq-git-sample/Variantcalling-and-Genefusion/data/test/snpEff/snpEff.config")
        path("/home/kothai/cq-git-sample/Variantcalling-and-Genefusion/data/test/snpEff/snpEff/data")

    output:
        tuple val(sampleId), path("${sampleId}.annotated.vcf"), path("${sampleId}.summary.html")

    script:
    """
    # Define absolute paths for SnpEff files
    SNPEFF_JAR="/home/kothai/cq-git-sample/Variantcalling-and-Genefusion/data/test/snpEff/snpEff.jar"
    SNPEFF_CONFIG="/home/kothai/cq-git-sample/Variantcalling-and-Genefusion/data/test/snpEff/snpEff.config"
    SNPEFF_DB_DIR="/home/kothai/cq-git-sample/Variantcalling-and-Genefusion/data/test/snpEff/snpEff/data"
    GENOME="${params.genomedb}"

    # Annotate variants using SnpEff
    java -Xmx16G -jar \$SNPEFF_JAR \\
        -c \$SNPEFF_CONFIG \\
        -v \$GENOME \\
        -dataDir \$SNPEFF_DB_DIR \\
        ${vcf} > ${sampleId}.annotated.vcf

    # Generate a summary HTML file
    java -Xmx16G -jar \$SNPEFF_JAR \\
        -c \$SNPEFF_CONFIG \\
        -v \$GENOME \\
        -dataDir \$SNPEFF_DB_DIR \\
        -stats ${sampleId}.summary.html \\
        ${vcf} > /dev/null
    """
}

