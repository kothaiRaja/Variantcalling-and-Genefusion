process TRIM_READS {
    tag { "${sample_id}_${task.process}" }
	label 'process_medium'
	
    container params.trim_reads_container 
    publishDir params.trim_reads_outdir


    input:
    tuple val(sample_id), path(r1), path(r2), val(strandedness)

    output:
    tuple val(sample_id), path("trimmed_${sample_id}_R1.fastq.gz"), path("trimmed_${sample_id}_R2.fastq.gz"), val(strandedness), emit: trimmed_reads
    tuple path ("${sample_id}_fastp.html"), path("${sample_id}_fastp.json"), emit: fastp_reports
	path("versions.yml"), emit: versions

    script:
    """
    fastp -i ${r1} -I ${r2} \
		-o trimmed_${sample_id}_R1.fastq.gz \
		-O trimmed_${sample_id}_R2.fastq.gz \
		--detect_adapter_for_pe ${params.trim_reads_detect_adapters} \
		--adapter_sequence ${params.trim_reads_adapter_sequence} \
		--adapter_sequence_r2 ${params.trim_reads_adapter_sequence_r2} \
		--length_required ${params.trim_reads_length_required} \
		--cut_front ${params.trim_reads_cut_front} \
		--cut_tail ${params.trim_reads_cut_tail} \
		--cut_window_size ${params.trim_reads_cut_window_size} \
		--cut_mean_quality ${params.trim_reads_cut_mean_quality} \
		--html ${sample_id}_fastp.html \
		--json ${sample_id}_fastp.json
	  
	#Capture the version
	fastp_version=\$(fastp --version 2>&1 | head -n1 | awk '{print \$2}')

	cat <<EOF > versions.yml
	"${task.process}":
	  fastp: "\${fastp_version}"
	EOF


    """
}