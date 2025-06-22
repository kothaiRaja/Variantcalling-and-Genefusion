include { CREATE_STAR_INDEX 		} from '../modules/references/STAR_index.nf'

workflow BUILD_STAR_INDEX {

    take:
    reference_genome
    gtf_annotation

    main:

    star_index_ch = Channel.empty()

    if (params.star_genome_index && file(params.star_genome_index).exists()) {
        println "Using provided STAR genome index: ${params.star_genome_index}"
        star_index_ch = Channel.value(file(params.star_genome_index))
    } else {
        println "Generating STAR genome index..."
        star_index_ch     = CREATE_STAR_INDEX(reference_genome, gtf_annotation)
         
    }

    emit:
    star_index = star_index_ch
}
