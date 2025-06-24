process ANNOTATE_VARIANTS {
    tag { "${sample_id}_${task.process}" }

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
    tuple val(sample_id), path("snpeff_annotated_${sample_id}.vcf"), emit: annotated_vcf
    path "snpeff_${sample_id}.summary.html", emit: summary_html
    path "snpeff_${sample_id}.snpeff_summary.csv", emit: summary
    path "versions.yml", emit: versions

    script:
    def avail_mem = task.memory?.giga ?: 3

"""
THREADS=${task.cpus}

echo "Annotating variants for sample: ${sample_id}"

config_file=\$(realpath ${snpEffConfig})
data_dir=\$(realpath ${snpEffDbDir})

if ! grep -q "^${genomedb}\\.genome" "\$config_file"; then
    echo "${genomedb}.genome : Homo_sapiens" >> "\$config_file"
fi

if ! grep -q "^database.repository" "\$config_file"; then
    echo "database.repository : https://snpeff.blob.core.windows.net/databases/" >> "\$config_file"
fi

ls -l \$data_dir/${genomedb} || echo "  Warning: SnpEff DB folder not found at expected location."

java -Xmx${avail_mem}G -jar ${snpEffJar} \\
    -v ${genomedb} \\
    -c \$config_file \\
    -dataDir \$data_dir \\
    -stats snpeff_${sample_id}.summary.html \\
    -csvStats snpeff_${sample_id}.snpeff_summary.csv \\
    ${vcf} > snpeff_annotated_${sample_id}.vcf

snpeff_version=\$(java -jar ${snpEffJar} -version | head -n 1)
cat <<EOF > versions.yml
"${task.process}":
  snpeff: "\${snpeff_version}"
EOF

echo "Annotation complete for sample: ${sample_id}"
"""


    
}
