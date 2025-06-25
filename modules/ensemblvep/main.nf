process ANNOTATEVARIANTS_VEP {
    tag { "${sample_id}_${task.process}" }

    label 'process_medium'
    container params.annotate_vep_container
    publishDir params.annotate_vep_outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(input_vcf), path(input_vcf_tbi)
    path vep_cache 
    val genome_assembly
    val cache_version
    val species
	path vep_plugins

    output:
    tuple val(sample_id), path("vep_annotated_${sample_id}.vcf"), emit: annotated_vcf  
    path("vep_annotated_${sample_id}.vep.html"), emit: summary
    path("versions.yml"), emit: versions

    script:
    def plugin_flags = params.vep_plugin_flags ?: ''
	def extra_flags  = params.vep_extra ?: ''

    def avail_mem = task.memory?.giga ?: 3
    def buffer_size = (avail_mem < 8) ? 100 : (avail_mem < 20 ? 200 : 500)

    """
    echo "Running VEP for sample: ${sample_id}"

    if [ ! -d "${vep_cache}" ]; then
        echo "ERROR: VEP cache directory does not exist at: \$(realpath ${vep_cache})" >&2
        exit 1
    fi

    vep \\
        --input_file "${input_vcf}" \\
        --output_file "vep_annotated_${sample_id}.vcf" \\
        --stats_file "vep_annotated_${sample_id}.vep.html" \\
        --cache \\
        --dir_cache "${vep_cache}" \\
        --dir_plugins "${vep_plugins}" \\
        ${plugin_flags} \\
		${extra_flags} \\
        --species "${species}" \\
        --assembly "${genome_assembly}" \\
        --cache_version ${cache_version} \\
        --format vcf \\
        --vcf \\
        --symbol \\
        --protein \\
        --check_existing \\
        --offline \\
        --force_overwrite

vep_version=\$(vep --help 2>&1 | grep 'ensembl-vep' | awk '{print \$NF}')

cat <<-END_VERSIONS > versions.yml
"${task.process}":
  ensemblvep: "\${vep_version}"
END_VERSIONS


    """
}
