include { DOWNLOAD_VEP_CACHE } from '../modules/references/VEP.nf'
include { DOWNLOAD_VEP_PLUGINS } from '../modules/references/VEP.nf'

workflow VEP_SETUP {

    
	

    main:

vep_cache_ch   = Channel.empty()
vep_plugins_ch = Channel.empty()

if (!params.vep_enable) {
    println "Skipping VEP setup as params.vep_enable is set to false"

    emit:
    vep_cache   = vep_cache_ch
    vep_plugins = vep_plugins_ch
    return
}

// VEP cache logic
if (params.vep_cache) {
    try {
        vep_cache_ch = Channel
            .fromPath(params.vep_cache, checkIfExists: true)
            .map { [it] }
            .collect()
        println "Using provided VEP cache: ${params.vep_cache}"
    } catch (Exception e) {
        println "Provided VEP cache path invalid → downloading cache"
        vep_cache_ch = DOWNLOAD_VEP_CACHE()
    }
} else {
    println "No VEP cache path provided → downloading"
    vep_cache_ch = DOWNLOAD_VEP_CACHE()
}

// VEP plugin logic
if (params.vep_plugins) {
    try {
        vep_plugins_ch = Channel
            .fromPath(params.vep_plugins, checkIfExists: true)
            .map { [it] }
            .collect()
        println "Using provided VEP plugins: ${params.vep_plugins}"
    } catch (Exception e) {
        println "Provided VEP plugin path invalid → downloading plugins"
        vep_plugins_ch = DOWNLOAD_VEP_PLUGINS().vep_plugins
    }
} else {
    println "No VEP plugin path provided → downloading"
    vep_plugins_ch = DOWNLOAD_VEP_PLUGINS().vep_plugins
}

emit:
vep_cache   = vep_cache_ch
vep_plugins = vep_plugins_ch

   
}
