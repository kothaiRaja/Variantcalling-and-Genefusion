process ANNOTATEVARIANTS_VEP {
    tag "Annotate variants using VEP"

    container 'https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0'
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(input_vcf), path(input_vcf_tbi)
    path vep_cache
    path clinvar_vcf        
    path clinvar_vcf_tbi    
    val genome_assembly
    val cache_version
    val species

    output:
    tuple val(sample_id), path("vep_annotated_${sample_id}.vcf"), emit: annotated_vcf  
    tuple val(sample_id), path("vep_annotated_${sample_id}.html"), emit: summary_html

    script:
    def args = task.ext.args ?: ''

    """
    # Debugging: Log resolved paths and check file existence
    echo "Input VCF: \$(realpath "${input_vcf}")"
    echo "Input VCF Index: \$(realpath "${input_vcf_tbi}")"
    echo "VEP Cache Directory: \$(realpath "${vep_cache}")"
    echo "ClinVar VCF: \$(realpath "${clinvar_vcf}")"
    echo "ClinVar Index: \$(realpath "${clinvar_vcf_tbi}")"
    echo "Genome Assembly: ${genome_assembly}"
    echo "Cache Version: ${cache_version}"
    echo "Species: ${species}"

    # Check if required input files and directories exist
    if [ ! -d "${vep_cache}" ]; then
        echo "Error: VEP cache directory does not exist at: \$(realpath "${vep_cache}")" >&2
        exit 1
    fi

    if [ ! -f "${clinvar_vcf}" ]; then
        echo "Error: ClinVar VCF file not found at: \$(realpath "${clinvar_vcf}")" >&2
        exit 1
    fi

    if [ ! -f "${clinvar_vcf_tbi}" ]; then
        echo "Error: ClinVar index file not found at: \$(realpath "${clinvar_vcf_tbi}")" >&2
        exit 1
    fi

    # Run VEP annotation
    vep \\
        --input_file "${input_vcf}" \\
        --output_file "vep_annotated_${sample_id}.vcf" \\
		$args \\
        --stats_file "vep_annotated_${sample_id}.html" \\
        --cache \\
        --dir_cache "${vep_cache}" \\
        --assembly "${genome_assembly}" \\
        --cache_version "${cache_version}" \\
        --species "${species}" \\
        --format vcf \\
        --vcf \\
        --symbol \\
        --protein \\
        --force_overwrite \\
        --check_existing \\
		--fields "Allele,Consequence,IMPACT,CLIN_SIG,CLNDN" \\
		--custom "${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN" 

    """
}
