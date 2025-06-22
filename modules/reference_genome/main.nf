process DOWNLOAD_REF_GENOME {
    tag "Download reference genome"
    label 'process_high'
    container null

    // Always copy to ref_base for reusability
    publishDir "${params.ref_base}/reference", mode: 'copy'

    output:
    path("genome.fa.gz"), emit: genome_gz

    // Only run if no path is given AND genome.fa.gz not already downloaded
    when:
    !params.reference_genome &&
    !file("${params.ref_base}/reference/genome.fa.gz").exists()

    script:
    """
    echo "Downloading reference genome from: ${params.genome_download_url}"
    wget -q -O genome.fa.gz ${params.genome_download_url}
    """
}
