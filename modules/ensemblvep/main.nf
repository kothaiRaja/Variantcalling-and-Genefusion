process ANNOTATEVARIANTS_VEP {
    tag "Annotate variants using VEP"

    label 'process_medium'
    container params.annotate_vep_container
    publishDir params.annotate_vep_outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(input_vcf), path(input_vcf_tbi)
    path vep_cache    
    val genome_assembly
    val cache_version
    val species

    output:
    tuple val(sample_id), path("vep_annotated_${sample_id}.vcf"), emit: annotated_vcf  
    path("vep_annotated_${sample_id}.html"), emit: summary_html
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    echo "Running VEP for sample: ${sample_id}"

    if [ ! -d "${vep_cache}" ]; then
        echo "ERROR: VEP cache directory does not exist at: \$(realpath ${vep_cache})" >&2
        exit 1
    fi

    vep \\
        --input_file "${input_vcf}" \\
        --output_file "vep_annotated_${sample_id}.vcf" \\
        --stats_file "vep_annotated_${sample_id}.html" \\
        --cache \\
        --dir_cache "${vep_cache}" \\
        --species "${species}" \\
        --assembly "${genome_assembly}" \\
        --cache_version ${cache_version} \\
        --format vcf \\
        --vcf \\
        --symbol \\
        --protein \\
        --check_existing \\
        --everything \\
        --filter_common \\
        --per_gene \\
        --total_length \\
        --force_overwrite \\
        --offline

    
cat <<-END_VERSIONS > versions.yml
"${task.process}":
  ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
END_VERSIONS
    """
}
