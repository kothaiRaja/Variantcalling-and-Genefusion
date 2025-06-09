process MAF_VISUALIZATION {

    tag { "${sample_id}_${task.process}" }
    label 'process_low'

    container params.maftools_visualisation_container
    publishDir params.maftools_visual_outdir, mode: "copy"

    input:
    tuple val(sample_id), path(maf_file)
    path r_script

    output:
    path "plots_${sample_id}", emit: maf_plots
    path "versions.yml", emit: versions

    script:
    """
    # Run the R plotting script
    Rscript ${r_script} ${maf_file} ${sample_id}

    # Get maftools version only
    MAFTOOLS_VERSION=\$(Rscript -e "cat(as.character(packageVersion('maftools')))")

    # Save version info
cat <<-EOF > versions.yml
"${task.process}":
  maftools: "\$MAFTOOLS_VERSION"
EOF
    """
}
