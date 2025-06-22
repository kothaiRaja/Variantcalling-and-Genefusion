include { DOWNLOAD_GTF } from '../modules/gtf_annotation/main.nf'
include { GUNZIP as GUNZIP_GTF} from '../modules/gunzip/main.nf'
include { GENERATEEXONS_BED } from '../modules/references/exons_bed.nf'

workflow DOWNLOAD_GTF_ANNOTATION {

    main:

    gtf_ch       = Channel.empty()
    exons_bed_ch = Channel.empty()

    // ======== Step 1: GTF file handling ========
    if (params.gtf_annotation) {
        println " Using GTF from config: ${params.gtf_annotation}"
        gtf_ch = Channel
            .fromPath(params.gtf_annotation, checkIfExists: true)
            .map { file -> [ file ] }
            .collect()
    }
    else if (file("${params.ref_base}/reference/annotation.gtf").exists()) {
        println "️  Reusing annotation.gtf from ref_base"
        gtf_ch = Channel
            .fromPath("${params.ref_base}/reference/annotation.gtf", checkIfExists: true)
            .map { file -> [ file ] }
            .collect()
    }
    else {
        println " Downloading annotation.gtf.gz and gunzipping..."
        dl = DOWNLOAD_GTF()
        gtf_ch = GUNZIP_GTF(dl)
        
    }

    // ======== Step 2: Exons BED file handling ========
    if (params.exons_bed) {
        println " Using exons BED from config: ${params.exons_bed}"
        exons_bed_ch = Channel
            .fromPath(params.exons_bed, checkIfExists: true)
            .map { file -> [ file ] }
            .collect()
    }
    else if (file("${params.ref_base}/reference/exons.bed").exists()) {
        println "️  Reusing exons.bed from ref_base"
        exons_bed_ch = Channel
            .fromPath("${params.ref_base}/reference/exons.bed", checkIfExists: true)
            .map { file -> [ file ] }
            .collect()
    }
    else {
        println "️  Generating exons.bed from GTF..."
        exons_bed_ch = GENERATEEXONS_BED(gtf_ch).collect()
         
    }

    emit:
    gtf        = gtf_ch
    exons_bed  = exons_bed_ch
}
