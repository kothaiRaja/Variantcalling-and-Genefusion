process GATK_VARIANT_FILTER {
    tag "variant_filter"

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/variant_filter", mode: "copy"

    input:
    path(vcf_file)
	path(vcf_index)  // Input VCF and index
	path(tsv_file)
    path genome 	// Reference genome
	path genome_index
	path genome_dict

    output:
          path("final.vcf.gz") 
          path("final.vcf.gz.tbi")  // Filtered VCF and index

    script:
    """
    gatk VariantFiltration \
    -R ${genome} \
    -V ${vcf_file} \
    --cluster-window-size 35 --cluster-size 3 \
    --filter-name "LowQual" --filter-expression "QUAL < 30.0" \
    --filter-name "LowQD" --filter-expression "QD < 2.0" \
    --filter-name "HighFS" --filter-expression "FS > 60.0" \
    -O final.vcf.gz

    """
}