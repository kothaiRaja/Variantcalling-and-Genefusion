process ARRIBA_VISUALIZATION {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/r-base%3A4.4.1"  
    publishDir "${params.outdir}/ARRIBA_VISUALIZATION", mode: 'copy'

    input:
    tuple val(sample_id), path(fusions_tsv)
	tuple val (sample_id), path(discarded_tsv)
    path r_script  
    path fasta  
    path gtf    

    output:
    path "*.fusion_plot_discarded.pdf"

    script:
    """
    # Extract prefix (filename without extension)
    PREFIX=\$(basename ${discarded_tsv} .fusions.discarded.tsv)

    # Debug: Log resolved paths
    echo "Resolved discarded file path: \$(realpath ${discarded_tsv})"
    echo "Resolved annotation file path: \$(realpath ${gtf})"
    

    # Verify the resolved paths exist
    if [ ! -f "\$(realpath ${discarded_tsv})" ]; then
        echo "Error: Discarded file not found at: \$(realpath ${discarded_tsv})"
        exit 1
    fi
    if [ ! -f "\$(realpath ${gtf})" ]; then
        echo "Error: Annotation file not found at: \$(realpath ${gtf})"
        exit 1
    fi

    # Run the R script with escaped paths
    Rscript draw_fusions.R \
        --fusions="\$(realpath ${discarded_tsv})" \
        --annotation="\$(realpath ${gtf})" \
        --output="\${PREFIX}.fusion_plot_discarded.pdf"
    
    """
}

