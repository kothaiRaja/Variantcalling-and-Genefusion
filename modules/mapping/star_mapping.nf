nextflow.enable.dsl = 2

// Define the STAR alignment process
process RNASEQ_MAPPING_STAR {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/star", mode: "copy"

    input:
    path genomeDir
    tuple val(replicateId), path(reads)

    output:
    tuple val(replicateId), path("Aligned.*.sortedByCoord.out.bam")

    script:
    """
    STAR --genomeDir $genomeDir \
         --readFilesIn ${reads[0]} ${reads[1]} \
         --runThreadN ${task.cpus} \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:GM12878 \
         --outFileNamePrefix Aligned.
    """
}

process SAMTOOLS_FLAGSTAT {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/flagstat", mode: "copy"

    input:
    tuple val(replicateId), path(bam)

    output:
    path "${replicateId}_flagstat.txt"

    script:
    """
    samtools flagstat ${bam} > ${replicateId}_flagstat.txt
    """
}

process SAMTOOLS_FILTER_INDEX {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/star/index", mode: "copy"

    input:
    tuple val(replicateId), path(bam)

    output:
    tuple val(replicateId), path("${bam.baseName}.uniq.bam"), path("${bam.baseName}.uniq.bam.bai")

    script:
    """
    (samtools view -H $bam; samtools view $bam | grep -w 'NH:i:1') \
    | samtools view -Sb - > ${bam.baseName}.uniq.bam
    samtools index ${bam.baseName}.uniq.bam
    """
}

process ADD_READ_GROUP {
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/star/RG", mode: "copy"

    input:
    tuple val(replicateId), path(bamFile), path(baiFile)

    output:
    tuple val(replicateId), path("rg_added.bam")

    script:
    """
    gatk AddOrReplaceReadGroups \
        -I $bamFile \
        -O rg_added.bam \
        -RGID '$replicateId' \
        -RGLB 'library' \
        -RGPL 'illumina' \
        -RGPU 'machine' \
        -RGSM 'GM12878'
    """
}

workflow star_mapping {
    // Take input channels
    take:
    trimmed_reads
    genome_index

    // Channel wiring for processes
    star_mapped_bam = RNASEQ_MAPPING_STAR(genome_index, trimmed_reads)
    flagstat_results = SAMTOOLS_FLAGSTAT(star_mapped_bam)
    filtered_bam = SAMTOOLS_FILTER_INDEX(star_mapped_bam)
    final_bam = ADD_READ_GROUP(filtered_bam)

    // Return results
    return: [
        final_bam: final_bam,
        flagstat_results: flagstat_results
    ]
}
