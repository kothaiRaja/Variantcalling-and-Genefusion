process VCF_TO_TABLE {
    tag "Convert VCF to Table"

    container 'https://depot.galaxyproject.org/singularity/cyvcf2%3A0.31.1--py312h68a07e8_1'
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path vcf_file
	path html
	path script_file

    output:
    path "filtered_variants.csv"
    

    script:
    """
	python ${script_file} ${vcf_file}
	
    """
}
