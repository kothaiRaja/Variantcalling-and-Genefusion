process ARRIBA {
    tag { sample_id }
    label 'process_high'

    container params.arriba_container
    publishDir params.arriba_outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path gtf
    path blacklist
    path known_fusions

    output:
  
    tuple val(sample_id), path("*.fusions.tsv"), emit: fusions
    tuple val(sample_id), path("*.fusions.discarded.tsv"), emit: fusions_discarded
	path("versions.yml"), emit: versions

    script:

"""
echo "Running Arriba for Sample: ${sample_id}"

arriba \\
     -x "${bam}" \\
     -a "${fasta}" \\
     -g "${gtf}" \\
     -b "${blacklist}" \\
     -k "${known_fusions}" \\
     -o "${sample_id}.fusions.tsv" \\
     -O "${sample_id}.fusions.discarded.tsv"

echo "Arriba finished for Sample: ${sample_id}"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
  arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\\s//')
END_VERSIONS
"""
}
