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
    path "*.fusion_plot.pdf"

    script:
    """
    # Extract prefix (filename without extension)
    PREFIX=\$(basename ${fusions_tsv} .fusions.tsv)

    # Debug: Log resolved paths
    echo "Resolved fusions file path: \$(realpath ${fusions_tsv})"
    echo "Resolved annotation file path: \$(realpath ${gtf})"
    

    # Verify the resolved paths exist
    if [ ! -f "\$(realpath ${fusions_tsv})" ]; then
        echo "Error: Fusions file not found at: \$(realpath ${fusions_tsv})"
        exit 1
    fi
    if [ ! -f "\$(realpath ${gtf})" ]; then
        echo "Error: Annotation file not found at: \$(realpath ${gtf})"
        exit 1
    fi
	
	# Verify the fusions file is not empty
    if [ ! -s ${fusions_tsv} ]; then
        echo "No fusions found for ${sample_id}. Skipping visualization." >&2
        touch "\${PREFIX}.fusion_plot.pdf"
        exit 0
    fi

    # Run the R script
    Rscript draw_fusions.R \
        --fusions "\$(realpath ${fusions_tsv})" \
        --annotation "\$(realpath ${gtf})" \
        --output "\${PREFIX}.fusion_plot.pdf"
    
    """
}

