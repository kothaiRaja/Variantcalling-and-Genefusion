process ANNOTATEVARIANTS_VEP {
    container 'https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0'
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path input_vcf          // Input VCF file
    path input_vcf_tbi      // Tabix index file for the VCF
	path tsv
    path vep_cache          // VEP cache directory
    path clinvar_vcf        // ClinVar VCF file
    path clinvar_vcf_tbi    // Tabix index file for ClinVar
	

    output:
    path "annotated_variants.vcf"          // Annotated VCF output
    path "annotated_variants.html"         // HTML summary report

    script:
    """
    # Debugging: Log resolved paths and existence
    echo "Input VCF: \$(realpath ${input_vcf})"
    echo "Input VCF Index: \$(realpath ${input_vcf_tbi})"
    echo "VEP Cache Directory: \$(realpath ${vep_cache})"
    echo "ClinVar VCF: \$(realpath ${clinvar_vcf})"
    echo "ClinVar Index: \$(realpath ${clinvar_vcf_tbi})"

    if [ ! -d "${vep_cache}" ]; then
        echo "Error: VEP cache directory does not exist at: \$(realpath ${vep_cache})" >&2
        exit 1
    fi

    if [ ! -f "${clinvar_vcf}" ]; then
        echo "Error: ClinVar VCF file not found at: \$(realpath ${clinvar_vcf})" >&2
        exit 1
    fi

    if [ ! -f "${clinvar_vcf_tbi}" ]; then
        echo "Error: ClinVar index file not found at: \$(realpath ${clinvar_vcf_tbi})" >&2
        exit 1
    fi

    # Run VEP annotation
    vep \
        --input_file ${input_vcf} \
        --output_file annotated_variants.vcf \
        --stats_file annotated_variants.html \
        --cache \
        --dir_cache ${vep_cache} \
        --assembly GRCh38 \
        --format vcf \
        --vcf \
        --symbol \
        --protein \
        --force_overwrite \
        --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN
    """
}
