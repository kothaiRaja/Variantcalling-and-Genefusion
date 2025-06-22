include { DOWNLOAD_VEP_CACHE } from '../modules/references/VEP.nf'
include { DOWNLOAD_VEP_PLUGINS } from '../modules/references/VEP.nf'

workflow VEP_SETUP {

    
	

    main:

    vep_cache_ch   = Channel.empty()
    vep_plugins_ch = Channel.empty()
	
	if ( !params.vep_enable ) {
        println " Skipping VEP setup as params.vep_enable is set to false"
		
		emit:
        vep_cache   = vep_cache_ch
        vep_plugins = vep_plugins_ch
        return
        
    }

    if (
        params.vep_cache && file(params.vep_cache).exists() &&
        params.vep_plugins && file(params.vep_plugins).exists()
    ) {
        println " Using provided VEP cache and plugins"
		println "  VEP cache:   ${params.vep_cache}"
		println "  VEP plugins: ${params.vep_plugins}"


        vep_cache_ch = Channel
            .fromPath(params.vep_cache, checkIfExists: true)
            .map { [it] }
            .collect()

        vep_plugins_ch = Channel
            .fromPath(params.vep_plugins, checkIfExists: true)
            .map { [it] }
            .collect()

    } else {
        println " Downloading VEP cache and plugins..."

      

        vep_cache_ch = DOWNLOAD_VEP_CACHE()
         

        plugins = DOWNLOAD_VEP_PLUGINS()
        vep_plugins_ch = plugins.vep_plugins
    }

    emit:
    vep_cache   = vep_cache_ch
    vep_plugins = vep_plugins_ch
}
