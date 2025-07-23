process MAF_VISUALIZATION {
    
    label 'process_low'

    container params.maftools_visualisation_container
    publishDir "${params.maftools_visual_outdir}", mode: "copy"

    input:
    tuple val(meta), path(maf_file), path(r_script)

    output:
    path "plots/*", emit: maf_plots
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p plots

    Rscript ${r_script} ${maf_file} ${meta}

    MAFTOOLS_VERSION=\$(Rscript -e "cat(as.character(packageVersion('maftools')))")

    cat <<-EOF > versions.yml
"${task.process}":
  maftools: "\$MAFTOOLS_VERSION"
EOF
    """
}
