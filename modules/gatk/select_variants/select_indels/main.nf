process SELECT_INDELs {
    tag "${sample_id}_select_indels"
	
	container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/selected_variants/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path genome
	path index
	path dict

    output:
    tuple val(sample_id), path("${sample_id}_indels.vcf.gz"), path("${sample_id}_indels.vcf.gz.tbi")

    script:
    """
    gatk SelectVariants \
        -R ${genome} \
        -V ${vcf_file} \
        --select-type-to-include INDEL \
        -O ${sample_id}_indels.vcf.gz

    gatk IndexFeatureFile -I ${sample_id}_indels.vcf.gz
    """
}
