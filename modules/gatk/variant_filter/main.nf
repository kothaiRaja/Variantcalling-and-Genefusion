process GATK_VARIANT_FILTER {
    tag { "${meta.id}_${task.process}" }
    label 'process_high'

    container params.gatk_container
    publishDir params.variant_filter_outdir, mode: "copy"

    input:
    tuple val(meta), path(vcf_file), path(vcf_index)
    path genome
    path genome_index
    path genome_dict

    output:
    tuple val(meta), path("filtered_${meta.id}.vcf.gz"), path("filtered_${meta.id}.vcf.gz.tbi"), emit: filtered_vcf
    tuple val(meta), path("hardfiltered_${meta.id}.vcf.gz"), path("hardfiltered_${meta.id}.vcf.gz.tbi"), emit: hardfiltered_vcf
    path("versions.yml"), emit: versions

    script:
    def mem   = task.memory?.giga ?: 6
    def win   = params.gatk_vf_window_size ?: 35
    def clu   = params.gatk_vf_cluster_size ?: 3
    def minDP = params.min_dp ?: 6
    def minGQ = params.min_gq ?: 20

    """
    set -euo pipefail
    SAMPLE='${meta.id}'
    INVCF="${vcf_file}"

    # Quick count before
    gatk CountVariants -V "\$INVCF" > pre_${meta.id}.count.txt || true

    # 1) Site-level hard filters (flags only)
    gatk --java-options "-Xmx${mem}g" VariantFiltration \
      -R "${genome}" \
      -V "\$INVCF" \
      --cluster-window-size ${win} \
      --cluster-size ${clu} \
      --filter-name "QD2_SNP"         --filter-expression "vc.isSNP()  && QD  < 2.0" \
      --filter-name "QD2_INDEL"       --filter-expression "vc.isIndel()&& QD  < 2.0" \
      --filter-name "FS30_SNP"        --filter-expression "vc.isSNP()  && FS  > 30.0" \
      --filter-name "FS200_INDEL"     --filter-expression "vc.isIndel()&& FS  > 200.0" \
      --filter-name "SOR3_SNP"        --filter-expression "vc.isSNP()  && SOR > 3.0" \
      --filter-name "SOR10_INDEL"     --filter-expression "vc.isIndel()&& SOR > 10.0" \
      --filter-name "MQ40_SNP"        --filter-expression "vc.isSNP()  && MQ  < 40.0" \
      --filter-name "MQRankSum12.5"   --filter-expression "vc.isSNP()  && MQRankSum < -12.5" \
      --filter-name "ReadPosRankSum8" --filter-expression "ReadPosRankSum < -8.0" \
      --filter-name "BaseQRankSum8"   --filter-expression "BaseQRankSum < -8.0" \
      -O "hardfiltered_${meta.id}.vcf.gz"

    gatk IndexFeatureFile -I "hardfiltered_${meta.id}.vcf.gz"
    gatk CountVariants -V "hardfiltered_${meta.id}.vcf.gz" > hard_${meta.id}.count.txt || true

    # 2) Minimal genotype filters (DP & GQ)
    gatk --java-options "-Xmx${mem}g" VariantFiltration \
      -R "${genome}" \
      -V "hardfiltered_${meta.id}.vcf.gz" \
      --genotype-filter-name "LowDP" --genotype-filter-expression "DP < ${minDP}" \
      --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ < ${minGQ}" \
      --set-filtered-genotype-to-no-call \
      -O "tmp_${meta.id}.gtfiltered.vcf.gz"

    gatk IndexFeatureFile -I "tmp_${meta.id}.gtfiltered.vcf.gz"
    gatk CountVariants -V "tmp_${meta.id}.gtfiltered.vcf.gz" > gt_${meta.id}.count.txt || true

    # 3) Keep only PASS sites; drop non-variants & unused ALTs
    gatk --java-options "-Xmx${mem}g" SelectVariants \
      -R "${genome}" \
      -V "tmp_${meta.id}.gtfiltered.vcf.gz" \
      --exclude-filtered true \
      --remove-unused-alternates true \
      --exclude-non-variants true \
      -O "filtered_${meta.id}.vcf.gz"

    gatk IndexFeatureFile -I "filtered_${meta.id}.vcf.gz"
    gatk CountVariants -V "filtered_${meta.id}.vcf.gz" > post_${meta.id}.count.txt || true

    # Versions
    gatk_version=\$(gatk --version | head -n 1)
    cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
