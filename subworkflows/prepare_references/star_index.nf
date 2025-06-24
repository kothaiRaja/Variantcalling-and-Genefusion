include { CREATE_STAR_INDEX 		} from '../../modules/references/STAR_index.nf'

workflow BUILD_STAR_INDEX {

    take:
    reference_genome
    gtf_annotation

    main:

    star_index_ch = Channel.empty()

    if (params.star_genome_index && file(params.star_genome_index).exists()) {
        println " Using STAR index from params: ${params.star_genome_index}"
        star_index_ch = Channel
            .fromPath(params.star_genome_index, checkIfExists: true)
            .map { file -> [ file ] }
            .collect()
    }
    else if (file("${params.ref_base}/reference/STAR_index").exists()) {
        println " Reusing STAR index from: ${params.ref_base}/reference/STARindex"
        star_index_ch = Channel
            .fromPath("${params.ref_base}/reference/STAR_index", checkIfExists: true)
            .map { file -> [ file ] }
            .collect()
    }
    else {
        println "  STAR index not found — generating now..."
        star_index_ch = CREATE_STAR_INDEX(reference_genome, gtf_annotation)
    }


    emit:
    star_index = star_index_ch
}
