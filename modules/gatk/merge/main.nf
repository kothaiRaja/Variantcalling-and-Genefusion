process GATK_MERGEVCFS {
    tag { sample_id }

    cpus params.get('gatk_merge_vcfs_cpus', 4)
    memory params.get('gatk_merge_vcfs_memory', '8GB')
    time params.get('gatk_merge_vcfs_time', '4h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/merged_vcf", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_list), path(tbi_list)

    output:
    tuple val(sample_id), path("merged_${sample_id}.vcf.gz"), path("merged_${sample_id}.vcf.gz.tbi")

    script:
    """
    gatk MergeVcfs \\
        ${vcf_list.collect { "-I ${it}" }.join(" \\\n")} \\
        -O merged_${sample_id}.vcf.gz

    gatk IndexFeatureFile -I merged_${sample_id}.vcf.gz
    """
}