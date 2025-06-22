include { DOWNLOAD_ARRIBA } from '../modules/references/Tools.nf'

workflow ARRIBA_SETUP {

    main:

    // ========== CASE 1: User provides arriba_dir ==========
    if (params.arriba_dir && file(params.arriba_dir).exists()) {
        println " Using user-provided Arriba tool directory: ${params.arriba_dir}"
        arriba_dir_ch = Channel.value(file(params.arriba_dir))
    } 

    // ========== CASE 2: Need to download Arriba and extract files ========== //
    else {
        println "  Downloading Arriba and extracting reference files"

        arriba_dir_ch = DOWNLOAD_ARRIBA()

               }

    emit:
    arriba_dir       = arriba_dir_ch
    
}
