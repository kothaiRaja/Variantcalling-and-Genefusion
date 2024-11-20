process RNASEQ_GATK_RECALIBRATE {
    tag "$replicateId"
    label "mem_large"
    container "/home/kothai/cq-git-sample/Praktikum/gatk4_4.2.6.0--hdfd78af_0.sif"
    publishDir "${params.outdir}/recalibrate", mode: 'copy'

    input:
    path genome
    path index
    path dict
    tuple val(replicateId), path(bam), path(bai)
    tuple path(variants_file), path(variants_file_index)

    output:
    tuple val(replicateId), path("${replicateId}.final.uniq.bam")

    script:
    """
    gatk BaseRecalibrator \
        -R $genome \
        -I $bam \
        --known-sites $variants_file \
        -O final.rnaseq.grp \
        --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true"

    gatk ApplyBQSR \
        -R $genome \
        -I $bam \
        --bqsr-recal-file final.rnaseq.grp \
        -O ${replicateId}.final.uniq.bam \
        --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
    """
}