include { DOWNLOAD_GTF }         from '../../modules/gtf_annotation/main.nf'
include { GUNZIP as GUNZIP_GTF } from '../../modules/gunzip/main.nf'
include { GENERATEEXONS_BED }    from '../../modules/references/exons_bed.nf'

workflow DOWNLOAD_GTF_ANNOTATION {

    main:

    gtf_ch        = Channel.empty()
    exons_bed_ch  = Channel.empty()

    // ======== Step 1: GTF file handling ========
    if (params.gtf_annotation) {
        println " Using GTF from config: ${params.gtf_annotation}"

        gtf_path_ch = Channel.fromPath(params.gtf_annotation, checkIfExists: true)

        if (params.gtf_annotation.endsWith('.gz')) {
            GUNZIP_GTF(gtf_path_ch)
            gtf_ch = GUNZIP_GTF.out.gunzip
        } else {
            gtf_ch = gtf_path_ch
        }

    } else if (file("${params.ref_base}/reference/annotations.gtf").exists()) {
        println " Reusing annotation.gtf from ref_base"
        gtf_ch = Channel.fromPath("${params.ref_base}/reference/annotations.gtf", checkIfExists: true)

    } else {
        println " Downloading annotation.gtf.gz and unzipping..."
        DOWNLOAD_GTF()
        GUNZIP_GTF(DOWNLOAD_GTF.out.gtf)
        gtf_ch = GUNZIP_GTF.out.gunzip
    }

    // ======== Step 2: Exons BED file handling ========
    if (params.exons_bed) {
        println " Using exons BED from config: ${params.exons_bed}"
        exons_bed_ch = Channel.fromPath(params.exons_bed, checkIfExists: true)

    } else if (file("${params.ref_base}/reference/exons.bed").exists()) {
        println " Reusing exons.bed from ref_base"
        exons_bed_ch = Channel.fromPath("${params.ref_base}/reference/exons.bed", checkIfExists: true)

    } else {
        println "️ Generating exons.bed from GTF..."
        GENERATEEXONS_BED(gtf_ch)
        exons_bed_ch = GENERATEEXONS_BED.out.exons_bed
    }

    emit:
    gtf        = gtf_ch
    exons_bed  = exons_bed_ch
}
