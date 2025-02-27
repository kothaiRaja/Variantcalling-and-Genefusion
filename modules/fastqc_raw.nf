process FASTQC_RAW {
    tag { sample_id }
	
	cpus params.get('fastqc_raw_cpus', 2)
    memory params.get('fastqc_raw_memory', '4 GB')
    time params.get('fastqc_raw_time', '1h')
	
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_fastqc.*"
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"

    input:
    tuple val(sample_id), path(r1), path(r2), val(strandedness)

    output:
    tuple val(sample_id), path("*.zip"), path("*.html"), val(strandedness), emit: qc_results

    script:
    """
    fastqc ${r1} ${r2} --outdir .
    """
}