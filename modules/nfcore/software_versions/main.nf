process CUSTOM_DUMPSOFTWAREVERSIONS {
    tag 'container_Version'
	label 'process_low'

    container params.multiqc_container
	publishDir "${params.software_versions_outdir}", mode: 'copy'


    input:
    path versions
    path dump_script

    output:
    path "software_versions.yml", emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml", emit: versions

    script:
    """
    python ${dump_script}
    cp software_versions.yml versions.yml
    """
}