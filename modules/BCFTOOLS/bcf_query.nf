process BCFTOOLS_QUERY {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.19--h8b25389_1"
    publishDir "${params.outdir}/variant_summary", mode: "copy"

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_vcf_index)

    output:
    path("filtered_variants_summary_${sample_id}.txt")

    script:
    """
    # Generate a summary of filtered variants
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${filtered_vcf} > filtered_variants_summary_${sample_id}.txt

    # Validate the output
    if [ ! -s filtered_variants_summary_${sample_id}.txt ]; then
        echo "Error: Summary of filtered variants is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}