process DOWNLOAD_VEP_CACHE {
    tag "Download VEP Cache"
    publishDir "${params.test_data_dir}/Tools/VEP", mode: 'copy'
    container null

    output:
    path "vep_cache", emit: vep_cache  // Ensure the output is correctly captured

    script:
    """
    mkdir -p vep_cache  # Extract inside the work directory first
    echo "Downloading VEP cache for GRCh38..."

    wget -q -O vep_cache.tar.gz https://ftp.ensembl.org/pub/release-109/variation/indexed_vep_cache/homo_sapiens_vep_109_GRCh38.tar.gz
    tar -xzvf vep_cache.tar.gz -C vep_cache --strip-components=1
    rm vep_cache.tar.gz
    """
}


process DOWNLOAD_CLINVAR {
    tag "Downloading ClinVar VCF"
    container null
    publishDir "${params.test_data_dir}/Tools/VEP", mode: 'copy'

    when:
    !file("${params.clinvar}").exists() || !file("${params.clinvartbi}").exists()

    output:
    path "clinvar.vcf.gz", emit: clinvar_vcf
    path "clinvar.vcf.gz.tbi", emit: clinvar_tbi

    script:
    """
    wget -O clinvar.vcf.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

    if [ ! -s clinvar.vcf.gz ]; then
        echo "Error: Failed to download ClinVar VCF file" >&2
        exit 1
    fi

    wget -O clinvar.vcf.gz.tbi ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

    if [ ! -s clinvar.vcf.gz.tbi ]; then
        echo "Error: Failed to download ClinVar VCF index file" >&2
        exit 1
    fi
    """
}
