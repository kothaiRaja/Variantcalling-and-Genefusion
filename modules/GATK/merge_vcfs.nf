process GATK_MERGE_VCFS {
    tag "Combine GVCFs"
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/combined_gvcf", mode: 'copy'

    input:
    path genome
	path genome_index
	path genome_dict
    tuple val(sample_id), path (vcfs), path(vcf_index)  // List of GVCFs

    output:
    path "merged_output.vcf.gz", emit: vcf          // Merged VCF file
    path "merged_output.vcf.gz.tbi", emit: vcf_idx  // Index for the merged VCF file

    script:
    """
    gatk MergeVcfs \
        ${vcfs.collect { "-I $it" }.join(' ')} \
        --OUTPUT merged_output.vcf.gz
    """
}