include { DOWNLOAD_ARRIBA } from '../../modules/references/Tools.nf'

workflow ARRIBA_SETUP {

    main:

    arriba_dir_ch = Channel.empty()

    def local_arriba_dir = "${params.ref_base}/Tools/ARRIBA/arriba_v${params.arriba_version}"

    // ========== CASE 1: User provides arriba_dir ==========
    if (params.arriba_dir && file(params.arriba_dir).exists()) {
        println " Using user-provided Arriba directory: ${params.arriba_dir}"
        arriba_dir_ch = Channel.value(file(params.arriba_dir))
    }

    // ========== CASE 2: Reuse previously downloaded Arriba ==========
    else if (file(local_arriba_dir).exists()) {
        println " Reusing Arriba from local cache: ${local_arriba_dir}"
        arriba_dir_ch = Channel.value(file(local_arriba_dir))
    }

    // ========== CASE 3: Download Arriba ==========
    else {
        println "️  Downloading Arriba v${params.arriba_version}..."
        arriba_dir_ch = DOWNLOAD_ARRIBA()
    }

    emit:
    arriba_dir = arriba_dir_ch
}
