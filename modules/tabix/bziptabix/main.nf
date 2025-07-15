process BGZIP_TABIX_ANNOTATIONS {

    tag { "${meta.id}_${task.process}" }
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0"
    publishDir "params.annotated_vcf_compressed", mode: 'copy'

    input:
    tuple val(meta), path(annotated_vcf)

    output:
    tuple val(meta), 
          path("annotated_${meta.id}.vcf.gz"), 
          path("annotated_${meta.id}.vcf.gz.tbi"), 
          emit: compressed_indexed
    path("versions.yml"), emit: versions

    script:
    """
    THREADS=${task.cpus}

    # Compress the annotated VCF file
    bgzip --threads \${THREADS} -c ${annotated_vcf} > annotated_${meta.id}.vcf.gz

    # Index the compressed VCF file
    tabix -p vcf annotated_${meta.id}.vcf.gz

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
