include { DOWNLOAD_REF_GENOME } from '../../modules/reference_genome/main.nf'
include { GUNZIP }             from '../../modules/gunzip/main.nf'
include { CREATE_GENOME_INDEX } from '../../modules/references/genome_index.nf'
include { CREATE_GENOME_DICT  } from '../../modules/references/genome_dict.nf'

workflow DOWNLOAD_REFERENCE_GENOME {

    main:

    genome_fa_ch    = Channel.empty()
    genome_index_ch = Channel.empty()
    genome_dict_ch  = Channel.empty()

    // ================= Reference Genome FASTA =================
    if (params.reference_genome) {
        println " Reference FASTA provided: ${params.reference_genome}"

        genome_input_ch = Channel.fromPath(params.reference_genome, checkIfExists: true)

        if (params.reference_genome.endsWith('.gz')) {
            GUNZIP(genome_input_ch)
            genome_fa_ch = GUNZIP.out.gunzip
        } else {
            genome_fa_ch = genome_input_ch.map { file -> [ file ] }.collect()
        }

    } else if (file("${params.ref_base}/reference/genome.fa").exists()) {
        println " Reference genome already exists at ${params.ref_base}/reference/genome.fa"

        genome_fa_ch = Channel
            .fromPath("${params.ref_base}/reference/genome.fa", checkIfExists: true)
            .map { file -> [ file ] }
            .collect()

    } else {
        println " Downloading reference genome..."
        DOWNLOAD_REF_GENOME()
        GUNZIP(DOWNLOAD_REF_GENOME.out.genome)
        genome_fa_ch = GUNZIP.out.gunzip
    }

    // ================= Genome Index (.fai) =================
    if (params.reference_genome_index) {
        println " Reference Index (.fai) provided: ${params.reference_genome_index}"
        genome_index_ch = Channel
            .fromPath(params.reference_genome_index, checkIfExists: true)
            .map { file -> [ file ] }
            .collect()

    } else if (file("${params.ref_base}/reference/genome.fa.fai").exists()) {
        println " Reference Index found locally: ${params.ref_base}/reference/genome.fa.fai"
        genome_index_ch = Channel
            .fromPath("${params.ref_base}/reference/genome.fa.fai", checkIfExists: true)
            .map { file -> [ file ] }
            .collect()

    } else {
        println "️ Creating genome index (.fai)..."
        CREATE_GENOME_INDEX(genome_fa_ch)
        genome_index_ch = CREATE_GENOME_INDEX.out.genome_index
    }

    // ================= Genome Dict (.dict) =================
    if (params.reference_genome_dict) {
        println " Reference Dict (.dict) provided: ${params.reference_genome_dict}"
        genome_dict_ch = Channel
            .fromPath(params.reference_genome_dict, checkIfExists: true)
            .map { file -> [ file ] }
            .collect()

    } else if (file("${params.ref_base}/reference/genome.dict").exists()) {
        println " Reference Dict found locally: ${params.ref_base}/reference/genome.dict"
        genome_dict_ch = Channel
            .fromPath("${params.ref_base}/reference/genome.dict", checkIfExists: true)
            .map { file -> [ file ] }
            .collect()

    } else {
        println " Creating genome dict (.dict)..."
        CREATE_GENOME_DICT(genome_fa_ch)
        genome_dict_ch = CREATE_GENOME_DICT.out.genome_dict
    }

    emit:
    genome       = genome_fa_ch
    genome_index = genome_index_ch
    genome_dict  = genome_dict_ch
}
