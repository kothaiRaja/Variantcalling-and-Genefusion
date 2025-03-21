process ARRIBA {
    tag { sample_id }
    label 'process_high'

    container params.arriba_container
    publishDir params.arriba_outdir, mode: 'copy'

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai)
    tuple val(sample_id), val(strandedness), path(chimeric_sam)
    path fasta
    path gtf
    path blacklist
    path known_fusions

    output:
    tuple val(sample_id), path("${sample_id}.fusions.tsv"), val(strandedness), emit: fusions
    tuple val(sample_id), path("${sample_id}.fusions.discarded.tsv"), val(strandedness), emit: fusions_discarded
    path("versions.yml"), emit: versions

    script:
    """
    echo "Running Arriba for Sample: ${sample_id}"

    arriba \\
         -x "${bam}" \\
         -c "${chimeric_sam}" \\
         -a "${fasta}" \\
         -g "${gtf}" \\
         -b "${blacklist}" \\
         -k "${known_fusions}" \\
         -o "${sample_id}.fusions.tsv" \\
         -O "${sample_id}.fusions.discarded.tsv"

    echo "Arriba finished for Sample: ${sample_id}"

    # Capture Arriba version
    arriba_version=\$(arriba --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' || echo "unknown")

    cat <<EOF > versions.yml
    "${task.process}":
      arriba: "\${arriba_version}"
    EOF
    """
}
