process FASTQC_RAW {
    tag { "${sample_id}_${task.process}" }
	label 'process_low'
	
	publishDir params.fastqc_outdir, mode: "copy", pattern: "*_fastqc.*"
    container params.fastqc_container

    input:
    tuple val(sample_id), path(r1), path(r2), val(strandedness)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), val(strandedness), emit: qc_results
	path("*_fastqc.html"), emit: fastqc_reports
	path("versions.yml"), emit: versions
 

    script:
    """
    fastqc ${r1} ${r2} --outdir .
	
	#Capture Versions
	fastqc_version=\$(fastqc --version | awk '{print \$2}')
	cat <<EOF > versions.yml
	"${task.process}":
	  fastqc: "\${fastqc_version}"
	EOF


    """
}