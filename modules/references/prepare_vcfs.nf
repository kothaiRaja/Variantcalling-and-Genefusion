process SPLIT_AND_INDEX_VCF {
    tag "Split SNPs and INDELs"
    label 'process_high'

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.ref_base}/reference", mode: 'copy'

    input:
    path all_variants_vcf

    output:
    path "known_snps.vcf.gz", emit: snps_vcf
    path "known_snps.vcf.gz.tbi", emit: snps_index
    path "known_indels.vcf.gz", emit: indels_vcf
    path "known_indels.vcf.gz.tbi", emit: indels_index

    script:
    """
    THREADS=${task.cpus}

    echo "Extracting SNPs..."
    bcftools view -v snps ${all_variants_vcf} -Oz -o known_snps.vcf.gz --threads \$THREADS
    tabix -p vcf known_snps.vcf.gz

    echo "Extracting INDELs..."
    bcftools view -v indels ${all_variants_vcf} -Oz -o known_indels.vcf.gz --threads \$THREADS
    tabix -p vcf known_indels.vcf.gz
    """
}

process INDEX_VCF {
    tag "Index VCF if missing"
    label 'process_light'

    container "https://depot.galaxyproject.org/singularity/htslib%3A1.1--h60f3df9_6"
    publishDir "${params.ref_base}/reference", mode: 'copy'

    input:
    path vcf_file
    val type

    output:
    path "known_${type}.vcf.gz.tbi", emit: index

    when:
    !file("known_${type}.vcf.gz.tbi").exists()

    script:
    """
    echo "Indexing known_${type}.vcf.gz ..."
    tabix -p vcf ${vcf_file}
    """
}