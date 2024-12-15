process ANNOTATE_VARIANTS {
    tag "Annotate variants"
    
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path vcf	// Filtered VCF file
	path index
    path snpEffJar  // Path to the SnpEff JAR file
    path snpEffConfig  // Path to the SnpEff configuration file
    path snpEffDbDir	// Path to the SnpEff database directory
	val genomedb

    output:
    path "annotated.vcf"
	path "annotated.summary.html"

    script:
	
    
    """
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${params.genomedb} \
        -dataDir ${snpEffDbDir} \
        ${vcf} > annotated.vcf

    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${params.genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated.summary.html \
        ${vcf} > /dev/null
    """
}