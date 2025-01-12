process GATK_VARIANT_FILTER {
    tag "${sample_id}_variant_filter"

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/variant_filter", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path genome
    path genome_index
    path genome_dict

    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi")

    script:
    """
    gatk VariantFiltration \
        -R ${genome} \
        -V ${vcf_file} \
        --cluster-window-size 35 \
        --cluster-size 3 \
        --filter-name "LowQual" --filter-expression "QUAL < 20.0" \
        --filter-name "LowQD" --filter-expression "QD < 1.5" \
        --filter-name "HighFS" --filter-expression "FS > 60.0" \
        --filter-name "LowMQ" --filter-expression "MQ < 30.0" \
        --filter-name "HighSOR" --filter-expression "SOR > 4.0" \
        --filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum < -5.0" \
        --filter-name "LowBaseQRankSum" --filter-expression "BaseQRankSum < -3.0" \
        -O ${sample_id}_filtered.vcf.gz

    # Validate the filtered VCF
    if [ ! -s ${sample_id}_filtered.vcf.gz ]; then
        echo "Error: Filtered VCF is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}
