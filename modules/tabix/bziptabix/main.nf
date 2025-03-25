process BGZIP_TABIX_ANNOTATIONS {
    tag "Compress & Index annotated variants"

    container "https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0"
    publishDir "${params.outdir}/compressed_annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(annotated_vcf)
    
    output:
    tuple val(sample_id), path("annotated_${sample_id}.vcf.gz"), path("annotated_${sample_id}.vcf.gz.tbi"), emit: compressed_indexed
	path("versions.yml"), emit: versions
    

    script:
    """
    THREADS=${task.cpus}

    # Compress the annotated VCF file
    bgzip --threads \${THREADS} -c ${annotated_vcf} > annotated_${sample_id}.vcf.gz

    # Index the compressed VCF file
    tabix -p vcf annotated_${sample_id}.vcf.gz
	
	# Capture versions
    bgzip_version=\$(bgzip --version 2>&1 | head -n 1)
    tabix_version=\$(tabix --version 2>&1 | head -n 1)

cat <<EOF > versions.yml
"${task.process}":
  bgzip: "\${bgzip_version}"
  tabix: "\${tabix_version}"
EOF

   
    """
}
