process ARRIBA {

    tag { "${meta.id}_${task.process}" }
    label 'process_high'

    container params.arriba_container
    publishDir params.arriba_outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path gtf
    path blacklist
    path known_fusions

    output:
    tuple val(meta), path("*.fusions.tsv"), emit: fusions
    tuple val(meta), path("*.fusions.discarded.tsv"), emit: fusions_discarded
    path("versions.yml"), emit: versions

    script:
    """
    echo "Running Arriba for Sample: ${meta.id}"

    arriba \\
         -x "${bam}" \\
         -a "${fasta}" \\
         -g "${gtf}" \\
         -b "${blacklist}" \\
         -k "${known_fusions}" \\
         -o "${meta.id}.fusions.tsv" \\
         -O "${meta.id}.fusions.discarded.tsv"

    echo "Arriba finished for Sample: ${meta.id}"

    arriba_version=\$(arriba -h 2>&1 | grep 'Version:' | sed 's/Version: //')

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
  arriba: "\${arriba_version}"
END_VERSIONS
    """
}
