process trim_reads {
    tag { sample_id }
	container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--hadf994f_3"
    publishDir "${params.outdir}/fastp", mode: "copy"
	storeDir "${params.test_data_dir}/fastp"
    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("trimmed_${sample_id}_R1.fastq.gz"), path("trimmed_${sample_id}_R2.fastq.gz"), emit: trimmed_reads
	tuple val(sample_id), path("${sample_id}_fastp.html"), path("${sample_id}_fastp.json"), emit: fastp_reports

    

    script:
    """
    fastp -i ${r1} -I ${r2} \
          -o trimmed_${sample_id}_R1.fastq.gz \
          -O trimmed_${sample_id}_R2.fastq.gz \
          --detect_adapter_for_pe \
          --html ${sample_id}_fastp.html \
          --json ${sample_id}_fastp.json
    """
}