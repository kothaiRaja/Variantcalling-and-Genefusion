process ANNOTATEVARIANTS_VEP {
    tag "Annotate variants using VEP"

    container 'https://depot.galaxyproject.org/singularity/ensembl-vep%3A113.0--pl5321h2a3209d_0'
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(input_vcf), path(input_vcf_tbi)
    path vep_cache    
    val genome_assembly
    val cache_version
    val species

    output:
    tuple val(sample_id), path("vep_annotated_${sample_id}.vcf"), emit: annotated_vcf  
    tuple val(sample_id), path("vep_annotated_${sample_id}.html"), emit: summary_html

    script:
    def args = task.ext.args ?: ''

    """
    
    # Debugging: Log resolved paths
    echo "Input VCF: \$(realpath ${input_vcf})"
    echo "VEP Cache Directory: \$(realpath ${vep_cache})"

    if [ ! -d "${vep_cache}" ]; then
        echo "ERROR: VEP cache directory does not exist at: \$(realpath ${vep_cache})" >&2
        exit 1
    fi

    ls -lh \$(realpath ${vep_cache})

    vep \
    --input_file "${input_vcf}" \
    --output_file "vep_annotated_${sample_id}.vcf" \
    --stats_file "vep_annotated_${sample_id}.html" \
    --cache \
    --dir_cache "${vep_cache}" \
    --format vcf \
    --vcf \
    --symbol \
    --protein \
    --force_overwrite \
    --check_existing \
    --everything \
    --filter_common \
    --per_gene \
    --total_length \
    --offline
    """ 

   
}
