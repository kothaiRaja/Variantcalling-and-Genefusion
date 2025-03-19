nextflow.enable.dsl = 2

//  Include required processes
include { ARRIBA } from '../modules/ARRIBA/ARRIBA_fusion/main.nf'
include { ARRIBA_VISUALIZATION } from '../modules/ARRIBA/ARRIBA_visualisation/main.nf'

workflow GENE_FUSION {
    
    take:
		STAR_bam_output
        STAR_chimeric_output
		reference_genome
        gtf_annotation		
        arriba_blacklist      
        arriba_known_fusions  
        arriba_scripts_dir    

    main:
        
        log.info "Starting Gene Fusion Detection Workflow..."

        

        // **Step 2: ARRIBA Fusion Detection**
        arriba_results = ARRIBA(
			STAR_bam_output,
            STAR_chimeric_output,
            reference_genome,
            gtf_annotation,
            arriba_blacklist,
            arriba_known_fusions
        )

        // **Step 3: ARRIBA Visualization**
        fusion_visuals = ARRIBA_VISUALIZATION(
            arriba_results[0],
            arriba_scripts_dir,
            reference_genome,
            gtf_annotation
        )

        log.info " Gene Fusion Detection Workflow Completed."

    emit:
        fusion_results = arriba_results.fusions
		discarded_results = arriba_results.fusions_discarded
        fusion_visualizations = fusion_visuals
}
