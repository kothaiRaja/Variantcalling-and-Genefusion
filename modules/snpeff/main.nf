process ANNOTATE_VARIANTS {
    tag "Annotate variants"
	label 'process_medium'

    container params.annotate_container_snpeff
    publishDir params.annotate_outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path(snpEffJar)
    path(snpEffConfig)
    path(snpEffDbDir)
    val(genomedb)

    output:
    tuple val(sample_id), path("annotated_${sample_id}.vcf"), emit: annotated_vcf
    path "annotated_${sample_id}.summary.html", emit: summary_html
    path "annotated_${sample_id}.csv", emit: annotation_csv 
	path("versions.yml"), emit: versions

    script:
    """
    THREADS=${task.cpus}

    echo "Annotating variants for sample: ${sample_id}"

    # Run SnpEff for annotation (VCF)
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        ${vcf} > annotated_${sample_id}.vcf

    # Generate annotation summary report (HTML)
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated_${sample_id}.summary.html \
        ${vcf} > /dev/null

    # Generate CSV file with annotation statistics
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -csvStats annotated_${sample_id}.csv \
        ${vcf} > /dev/null
		
	
    # Capture SnpEff version
	snpeff_version=\$(java -jar "${snpEffJar}" -version 2>&1 | awk 'NR==1')


    cat <<EOF > versions.yml
    "${task.process}":
      snpEff: "\${snpeff_version}"
    EOF

    echo "Annotation complete for sample: ${sample_id}"

    
    """
}
