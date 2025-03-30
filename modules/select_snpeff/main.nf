process ANNOTATE_VARIANTS_SELECTED {
    tag { "${sample_id}_${task.process}" }

    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotated_variants/${sample_id}", mode: "copy"

    cpus params.get('annotate_cpus', 8)
    memory params.get('annotate_memory', '16GB')
    time params.get('annotate_time', '2h')

    input:
    tuple val(sample_id), path(snp_vcf), path(snp_vcf_index), path(indel_vcf), path(indel_vcf_index)
    path(snpEffJar)
    path(snpEffConfig)
    path(snpEffDbDir)
    val(genomedb)

    output:
    tuple val(sample_id), path("${sample_id}_annotated_snps.vcf")
    tuple val(sample_id), path("${sample_id}_annotated_indels.vcf")
    path "annotated_${sample_id}_annotated_snps.summary.html"
    path "annotated_${sample_id}_annotated_indels.summary.html"
    path "annotated_${sample_id}_annotated_snps.csv"
    path "annotated_${sample_id}_annotated_indels.csv"

    script:
    """
    THREADS=${task.cpus}

    # Annotate SNP variants and generate annotated VCF file
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -csvStats annotated_${sample_id}_annotated_snps.csv \
        ${snp_vcf} > ${sample_id}_annotated_snps.vcf

    # Annotate INDEL variants and generate annotated VCF file
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -csvStats annotated_${sample_id}_annotated_indels.csv \
        ${indel_vcf} > ${sample_id}_annotated_indels.vcf

    # Generate annotation summary report for SNPs
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated_${sample_id}_annotated_snps.summary.html \
        ${snp_vcf} > /dev/null
    
    # Generate annotation summary report for INDELs
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated_${sample_id}_annotated_indels.summary.html \
        ${indel_vcf} > /dev/null
    """
}
