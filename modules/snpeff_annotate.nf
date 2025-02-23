process ANNOTATE_INDIVIDUAL_VARIANTS {
    tag "${sample_id}_annotate"

    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.resultsdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_index), val(strandedness)
	path (snpEffJar)
	path (snpEffConfig)
	path (snpEffDbDir)
	val (genomedb)

    
    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf"),val(strandedness), emit: annotated_vcf
    path("${sample_id}.annotated.summary.html"), emit: summary_html
    


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

process ANNOTATE_VARIANTS {
    tag "Annotate variants"
    
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path (vcf)	
	path (index)
	path (tsv)
	path (snpEffJar)  
    path (snpEffConfig)  
    path (snpEffDbDir)	
	val (genomedb)
	

    output:
    path "annotated.vcf", emit: annotated_vcf 
	path "annotated.summary.html", emit: summary_html

    script:
	
    
    """
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        ${vcf} > annotated.vcf

    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated.summary.html \
        ${vcf} > /dev/null
    """
}