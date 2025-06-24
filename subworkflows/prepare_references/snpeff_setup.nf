include { DOWNLOAD_SNPEFF_TOOL } from '../../modules/references/Tools.nf'
include { DOWNLOAD_SNPEFF_DB   } from '../../modules/references/Tools.nf'

workflow SNPEFF_SETUP {

    take:
    genome_id

    main:

    snpeff_jar_ch    = Channel.empty()
    snpeff_config_ch = Channel.empty()
    snpeff_db_ch     = Channel.empty()

    def local_base   = "${params.ref_base}/Tools/snpEff"
    def local_jar    = "${local_base}/snpEff.jar"
    def local_config = "${local_base}/snpEff.config"
    def local_db_dir = "${local_base}/snpEff/data"
    

    // ========== Step 1: Check for snpEff jar + config ==========
    def jar_and_config_available = false

    if (params.snpeff_jar && params.snpeff_config &&
        file(params.snpeff_jar).exists() && file(params.snpeff_config).exists()) {

        println " Using snpEff jar and config from config:"
        println "   JAR:    ${params.snpeff_jar}"
        println "   CONFIG: ${params.snpeff_config}"

        snpeff_jar_ch    = Channel.fromPath(params.snpeff_jar, checkIfExists: true).map { [it] }.collect()
        snpeff_config_ch = Channel.fromPath(params.snpeff_config, checkIfExists: true).map { [it] }.collect()
        jar_and_config_available = true

    } else if (file(local_jar).exists() && file(local_config).exists()) {

        println " Reusing snpEff jar and config from local cache:"
        println "   JAR:    ${local_jar}"
        println "   CONFIG: ${local_config}"

        snpeff_jar_ch    = Channel.fromPath(local_jar, checkIfExists: true).map { [it] }.collect()
        snpeff_config_ch = Channel.fromPath(local_config, checkIfExists: true).map { [it] }.collect()
        jar_and_config_available = true
    }

    if (!jar_and_config_available) {
        println "️ snpEff jar or config not found — downloading tool..."
        def tool_outputs = DOWNLOAD_SNPEFF_TOOL()
        snpeff_jar_ch    = tool_outputs.snpeff_jar
        snpeff_config_ch = tool_outputs.snpeff_config
    }

    // ========== Step 2: Check for snpEff database ==========
    if (params.snpeff_db_dir) {
        println " Using snpEff DB from config: ${params.snpeff_db_dir}"
        snpeff_db_ch = Channel.fromPath(params.snpeff_db_dir, checkIfExists: true).map { [it] }.collect()

    } else if (file(local_db).exists()) {
        println " Reusing local snpEff DB from ref_base: ${local_db}"
        snpeff_db_ch = Channel.fromPath(local_db_dir, checkIfExists: true).map { [it] }.collect()

    } else {
        println "️ snpEff database not found — downloading database for ${genome_id}..."
        snpeff_db_ch = DOWNLOAD_SNPEFF_DB(genome_id,snpeff_jar_ch,snpeff_config_ch)
    }

    emit:
    snpeff_jar    = snpeff_jar_ch
    snpeff_config = snpeff_config_ch
    snpeff_db_dir = snpeff_db_ch
}
