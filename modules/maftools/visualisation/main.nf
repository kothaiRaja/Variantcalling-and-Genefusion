process MAF_VISUALIZATION {
    tag { "${sample_id}" }
    label 'process_low'
    container params.maftools_visualisation_container
    publishDir "${params.maftools_visual_outdir}/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(maf_file),path(r_script)
   

    output:
    path "plots/*", emit: maf_plots
    path "versions.yml", emit: versions

    script:
    """
    # Create the exact directory structure Nextflow expects
    mkdir -p plots
    
    # Run with explicit output directory
    Rscript ${r_script} ${maf_file} ${sample_id}
    
    MAFTOOLS_VERSION=\$(Rscript -e "cat(as.character(packageVersion('maftools')))")
    cat <<-EOF > versions.yml
    "${task.process}":
      maftools: "\$MAFTOOLS_VERSION"
    EOF
    """
}