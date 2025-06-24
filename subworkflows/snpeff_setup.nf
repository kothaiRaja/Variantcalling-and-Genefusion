include { DOWNLOAD_SNPEFF_TOOL } from '../modules/references/Tools.nf'
include { DOWNLOAD_SNPEFF_DB   } from '../modules/references/Tools.nf'

workflow SNPEFF_SETUP {

    take:
    genome_id

    main:

    snpeff_jar_ch = Channel.empty()
    snpeff_config_ch = Channel.empty()
    snpeff_db_ch = Channel.empty()

    def local_base     = "${params.ref_base}/Tools/snpEff"
    def local_db_path  = "${local_base}/data/${genome_id}/snpEffectPredictor.bin"
    def local_jar      = "${local_base}/snpEff.jar"
    def local_config   = "${local_base}/snpEff.config"

    if (
        params.snpeff_jar && file(params.snpeff_jar).exists() &&
        params.snpeff_config && file(params.snpeff_config).exists() &&
        params.snpeff_db_dir && file("${params.snpeff_db_dir}/${genome_id}/snpEffectPredictor.bin").exists()
    ) {
        println " Using snpEff from provided parameters"
        snpeff_jar_ch    = Channel.value(file(params.snpeff_jar))
        snpeff_config_ch = Channel.value(file(params.snpeff_config))
        snpeff_db_ch     = Channel.value(file(params.snpeff_db_dir))

    } else if (
        file(local_jar).exists() &&
        file(local_config).exists() &&
        file(local_db_path).exists()
    ) {
        println " Using snpEff from local cache in ref_base"
        snpeff_jar_ch    = Channel.value(file(local_jar))
        snpeff_config_ch = Channel.value(file(local_config))
        snpeff_db_ch     = Channel.value(file("${local_base}/data"))

    } else {
        println "️  Downloading snpEff tool and database..."
        tool_outputs     = DOWNLOAD_SNPEFF_TOOL()
        snpeff_jar_ch    = tool_outputs.snpeff_jar
        snpeff_config_ch = tool_outputs.snpeff_config
        snpeff_db_ch     = DOWNLOAD_SNPEFF_DB(genome_id, tool_outputs)
    }

    emit:
    snpeff_jar    = snpeff_jar_ch
    snpeff_config = snpeff_config_ch
    snpeff_db_dir = snpeff_db_ch
}
