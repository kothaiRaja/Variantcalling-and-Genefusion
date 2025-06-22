include { DOWNLOAD_SNPEFF_TOOL } from '../modules/references/Tools.nf'
include { DOWNLOAD_SNPEFF_DB   } from '../modules/references/Tools.nf'


workflow SNPEFF_SETUP {

    take:
    genome_id

    main:

    if (
        params.snpeff_jar && file(params.snpeff_jar).exists() &&
        params.snpeff_config && file(params.snpeff_config).exists() &&
        params.snpeff_db_dir && file("${params.snpeff_db_dir}/${genome_id}/snpEffectPredictor.bin").exists()) 
		{
        println "Using user-provided snpEff resources"
		println "  snpEff JAR:     ${params.snpeff_jar}"
		println "  snpEff config:  ${params.snpeff_config}"
		println "  snpEff DB dir:  ${params.snpeff_db_dir}/${genome_id}"

        snpeff_jar_ch    = Channel.value(file(params.snpeff_jar))
        snpeff_config_ch = Channel.value(file(params.snpeff_config))
        snpeff_db_ch = Channel.value(file(params.snpeff_db_dir))


    } else {
        println "Downloading snpEff tool and database"

        tool_outputs = DOWNLOAD_SNPEFF_TOOL()

        snpeff_jar_ch    = tool_outputs.snpeff_jar
        snpeff_config_ch = tool_outputs.snpeff_config

        //  Combine the 2 path channels into a tuple
        

        snpeff_db_ch = DOWNLOAD_SNPEFF_DB(genome_id, tool_outputs)

     
    }
	
	emit:
    snpeff_jar    = snpeff_jar_ch
    snpeff_config = snpeff_config_ch
    snpeff_db_dir = snpeff_db_ch

}