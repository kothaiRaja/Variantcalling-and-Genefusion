process RNASEQ_CALL_VARIANTS {
    tag "$sampleId"
    label "mem_xlarge"
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
        path genome         // Path to the genome reference (e.g., genome.fa)
        path genome_fai     // Path to the .fai index file for the genome
        path genome_dict    // Path to the genome dictionary file (e.g., genome.dict)
        tuple val(sampleId), path(bam), path(bai)

    output:
        tuple val(sampleId), path("final_${sampleId}.vcf")

    script:
    """
	
	# Print paths for debugging
    echo "Genome: ${genome}"
    echo "BAM file: ${bam}"
    echo "Output: output_${sampleId}.vcf.gz"
	
	
    # Adjust genome dictionary path to absolute if necessary
    sed -i 's|UR:file:.*|UR:file:/home/kothai/cq-git-sample/Praktikum/data/genome.fa|' ${genome_dict}

    # Variant calling with HaplotypeCaller
    gatk HaplotypeCaller \
        --native-pair-hmm-threads ${task.cpus} \
        --reference ${genome} \
        --output output_${sampleId}.vcf.gz \
        -I ${bam} \
        --standard-min-confidence-threshold-for-calling 20.0 \
        --dont-use-soft-clipped-bases

    # Variant filtering with VariantFiltration
    gatk VariantFiltration \
        -R ${genome} -V output_${sampleId}.vcf.gz \
        --cluster-window-size 35 --cluster-size 3 \
        --filter-name FS --filter-expression "FS > 30.0" \
        --filter-name QD --filter-expression "QD < 2.0" \
        -O final_${sampleId}.vcf
    """
}