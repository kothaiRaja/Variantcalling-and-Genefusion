process RESET_READGROUPS {
    tag { "${meta.id}_${task.process}" }
    label 'process_medium'

    container params.gatk_container
    publishDir params.rg_reset_outdir, mode: "copy"

    input:
    tuple val(meta), path(merged_bam), path(merged_bai)

    output:
    tuple val(meta),
          path("${meta.id}_merged_fixedRG.bam"),
          path("${meta.id}_merged_fixedRG.bai"), emit: fixed_bams
    path("versions.yml"), emit: versions

    script:
    def sample_id = meta.id
    def rgid = sample_id
    def rglb = "library"
    def rgpl = params.get('seq_platform', 'ILLUMINA')
    def rgpu = "machine"
    def rgsm = sample_id
    def rgcn = params.get('seq_center', 'Unknown')

    return """
    echo "Resetting read groups for ${sample_id}"

    gatk AddOrReplaceReadGroups \\
        -I ${merged_bam} \\
        -O ${sample_id}_merged_fixedRG.bam \\
        -RGID ${rgid} \\
        -RGLB ${rglb} \\
        -RGPL ${rgpl} \\
        -RGPU ${rgpu} \\
        -RGSM ${rgsm} \\
        -RGCN ${rgcn}

    gatk BuildBamIndex \\
        -I ${sample_id}_merged_fixedRG.bam

    gatk_version=\$(gatk --version | head -n 1)

cat <<EOF > versions.yml
"${task.process}":
  gatk: "\${gatk_version}"
EOF
    """
}
