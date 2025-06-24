include { DOWNLOAD_VEP_CACHE } from '../../modules/references/VEP.nf'
include { DOWNLOAD_VEP_PLUGINS } from '../../modules/references/VEP.nf'

workflow VEP_SETUP {

    main:

    vep_cache_ch   = Channel.empty()
    vep_plugins_ch = Channel.empty()

    def local_cache_path   = "${params.ref_base}/Tools/VEP/cache"
    def local_plugins_path = "${params.ref_base}/Tools/VEP/plugins"

    

    // ========== VEP CACHE ==========
    if (params.vep_cache && file(params.vep_cache).exists()) {
        println " Using provided VEP cache: ${params.vep_cache}"
        vep_cache_ch = Channel
            .fromPath(params.vep_cache, checkIfExists: true)
            .map { [it] }
            .collect()

    } else if (file(local_cache_path).exists()) {
        println " Reusing VEP cache from: ${local_cache_path}"
        vep_cache_ch = Channel
            .fromPath(local_cache_path, checkIfExists: true)
            .map { [it] }
            .collect()

    } else {
        println "️ Downloading VEP cache..."
        vep_cache_ch = DOWNLOAD_VEP_CACHE().vep_cache
    }

    // ========== VEP PLUGINS ==========
    if (params.vep_plugins && file(params.vep_plugins).exists()) {
        println " Using provided VEP plugins: ${params.vep_plugins}"
        vep_plugins_ch = Channel
            .fromPath(params.vep_plugins, checkIfExists: true)
            .map { [it] }
            .collect()

    } else if (file(local_plugins_path).exists()) {
        println " Reusing VEP plugins from: ${local_plugins_path}"
        vep_plugins_ch = Channel
            .fromPath(local_plugins_path, checkIfExists: true)
            .map { [it] }
            .collect()

    } else {
        println "️ Downloading VEP plugins..."
        vep_plugins_ch = DOWNLOAD_VEP_PLUGINS().vep_plugins
    }

    emit:
    vep_cache   = vep_cache_ch
    vep_plugins = vep_plugins_ch
}
