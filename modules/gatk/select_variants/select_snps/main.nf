process SELECT_SNPs {
    tag "${sample_id}_select_snps"
    label 'process_medium'
    
    container params.gatk_container
    publishDir params.snp_select_outdir, mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path genome
    path index
    path dict

    output:
    tuple val(sample_id), path("${sample_id}_snps.vcf.gz"), path("${sample_id}_snps.vcf.gz.tbi"), emit: selected_snps
    path("versions.yml"), emit: versions

    script:
    """
    echo "Selecting SNPs for sample: ${sample_id}"

    gatk SelectVariants \\
        -R "${genome}" \\
        -V "${vcf_file}" \\
        --select-type-to-include SNP \\
        -O "${sample_id}_snps.vcf.gz"

    gatk IndexFeatureFile -I "${sample_id}_snps.vcf.gz"

    # Capture version
    gatk_version=\$(gatk --version | head -n 1)
	
cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
