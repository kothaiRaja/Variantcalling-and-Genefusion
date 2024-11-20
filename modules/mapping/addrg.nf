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