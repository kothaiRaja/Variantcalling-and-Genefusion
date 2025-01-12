process ANNOTATEVARIANTS_VEP {
    container 'https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0'

    input:
    path "input.vcf.gz"          // Input VCF file
    path "input.vcf.gz.tbi"      // Tabix index file for the VCF
	path tsv
    path "vep_cache"             // VEP cache directory
    path "clinvar.vcf.gz"        // ClinVar VCF file
    path "clinvar.vcf.gz.tbi"    // Tabix index file for ClinVar

    output:
    path "annotated_variants.vcf"          // Annotated VCF output
    path "annotated_variants.html"         // HTML summary report

    script:
    """
    vep --input_file input.vcf.gz \
        --output_file annotated_variants.vcf \
        --stats_file annotated_variants.html \
        --cache \
        --dir_cache vep_cache \
        --assembly GRCh38 \
        --format vcf \
        --vcf \
        --symbol \
        --protein \
        --force_overwrite \
        --custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN
    """
}