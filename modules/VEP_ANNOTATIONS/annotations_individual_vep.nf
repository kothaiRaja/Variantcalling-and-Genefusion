process ANNOTATE_INDIVIDUAL_VARIANTS_VEP {
    tag "${sample_id}_vep_annotate"

    container "https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_index)  // Input VCF and its index
    path vep_cache                                                  // VEP cache directory
    path clinvar_vcf                                                // ClinVar VCF file
    path clinvar_index                                              // ClinVar Tabix index file

    output:
    tuple val(sample_id), path("${sample_id}.vep.annotated.vcf"), path("${sample_id}.vep.summary.html")

    script:
    """
    # Annotate using Ensembl VEP with ClinVar
    vep --input_file ${filtered_vcf} \
        --output_file ${sample_id}.vep.annotated.vcf \
        --stats_file ${sample_id}.vep.summary.html \
        --cache \
        --dir_cache ${vep_cache} \
        --assembly GRCh38 \
        --format vcf \
        --vcf \
        --symbol \
        --protein \
        --force_overwrite \
        --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN

    # Validate that the annotated VCF is not empty
    if [ ! -s ${sample_id}.vep.annotated.vcf ]; then
        echo "Error: VEP annotation output is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}