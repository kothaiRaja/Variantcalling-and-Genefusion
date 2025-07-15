process VCF2MAF {

    tag { "${meta.id}_${task.process}" }
    label 'process_low'

    container params.maftools_container
    publishDir params.maftools_outdir, mode: "copy"

    input:
    tuple val(meta), path(vcf)
    path fasta
    path vep_cache

    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}_${vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, '')}"
    def args = "--species homo_sapiens --ncbi-build GRCh38"
    def vep_path = "--vep-path \$(dirname \$(which vep))"
    def vep_data = "--vep-data $vep_cache"

    """
    set -euo pipefail

    echo "VCF FILE: $vcf"
    echo "FASTA: $fasta"
    echo "VEP CACHE: $vep_cache"
    echo "Sample ID: ${meta.id}"

    vcf2maf.pl \\
      --input-vcf $vcf \\
      --output-maf ${prefix}.maf \\
      --tumor-id ${meta.id} \\
      --vcf-tumor-id ${meta.id} \\
      --ref-fasta $fasta \\
      $vep_path \\
      $vep_data \\
      $args

    VCF2MAF_VERSION=\$(vcf2maf.pl --help 2>&1 | grep -i 'vcf2maf' | grep -oP 'v[0-9.]+' || echo "not_detected")
    VEP_VERSION=\$(vep --help 2>&1 | grep -i 'ensembl-vep' | grep -oP '[0-9.]+' || echo "not_detected")

cat <<-END_VERSIONS > versions.yml
"${task.process}":
 vcf2maf: "\${VCF2MAF_VERSION}"
 ensembl-vep: "\${VEP_VERSION}"
END_VERSIONS
    """
}
