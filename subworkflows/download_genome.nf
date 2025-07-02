include { DOWNLOAD_REF_GENOME } from '../modules/reference_genome/main.nf'
include { GUNZIP } from '../modules/gunzip/main.nf'
include { CREATE_GENOME_INDEX } from '../modules/references/genome_index.nf'
include { CREATE_GENOME_DICT  } from '../modules/references/genome_dict.nf'

workflow DOWNLOAD_REFERENCE_GENOME {

    main:

    genome_fa_ch    = Channel.empty()
    genome_index_ch = Channel.empty()
    genome_dict_ch  = Channel.empty()

    // ================= Reference Genome FASTA =================
    if (params.reference_genome) {
        println " Reference FASTA provided: ${params.reference_genome}"

        if (params.reference_genome.endsWith('.gz')) {
            println " Detected compressed reference genome. Will decompress..."
            def gz_fa_ch = Channel.fromPath(params.reference_genome, checkIfExists: true)
            genome_fa_ch = GUNZIP(gz_fa_ch)
        } else {
            genome_fa_ch = Channel.fromPath(params.reference_genome, checkIfExists: true)
        }

    } else if (file("${params.ref_base}/reference/genome.fa").exists()) {
        println " Reference genome already exists at ${params.ref_base}/reference/genome.fa"
        genome_fa_ch = Channel.fromPath("${params.ref_base}/reference/genome.fa", checkIfExists: true)

    } else {
        println " Downloading reference genome..."
        def download_ref = DOWNLOAD_REF_GENOME()
        genome_fa_ch = GUNZIP(download_ref)
    }

    // ================= Genome Index (.fai) =================
    if (params.reference_genome_index) {
        println " Reference Index (.fai) provided: ${params.reference_genome_index}"
        genome_index_ch = Channel.fromPath(params.reference_genome_index, checkIfExists: true)

    } else if (file("${params.reference_genome}.fai").exists()) {
        println " Reference Index found: ${params.reference_genome}.fai"
        genome_index_ch = Channel.fromPath("${params.reference_genome}.fai", checkIfExists: true)

    } else {
        println " Creating genome index (.fai)..."
        genome_index_ch = CREATE_GENOME_INDEX(genome_fa_ch)
    }

    // ================= Genome Dict (.dict) =================
    if (params.reference_genome_dict) {
        println " Reference Dict (.dict) provided: ${params.reference_genome_dict}"
        genome_dict_ch = Channel.fromPath(params.reference_genome_dict, checkIfExists: true)

    } else if (file("${params.reference_genome}".replaceAll(/\.fa(sta)?$/, ".dict")).exists()) {
        def dict_path = "${params.reference_genome}".replaceAll(/\.fa(sta)?$/, ".dict")
        println " Reference Dict found: ${dict_path}"
        genome_dict_ch = Channel.fromPath(dict_path, checkIfExists: true)

    } else {
        println " Creating genome dict (.dict)..."
        genome_dict_ch = CREATE_GENOME_DICT(genome_fa_ch)
    }

    emit:
    genome       = genome_fa_ch
    genome_index = genome_index_ch
    genome_dict  = genome_dict_ch
}
