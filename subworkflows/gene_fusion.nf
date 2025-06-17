nextflow.enable.dsl = 2

// Include required processes
include { ARRIBA } from '../modules/ARRIBA/ARRIBA_fusion/main.nf'
include { ARRIBA_VISUALIZATION } from '../modules/ARRIBA/ARRIBA_visualisation/main.nf'

workflow GENE_FUSION {

    take:
        STAR_bam_output
        reference_genome
        gtf_annotation
        arriba_blacklist
        arriba_known_fusions
		cytobands
		protein_domains

    main:

        log.info " Starting Gene Fusion Detection Workflow..."

        ch_versions = Channel.empty()

        // Step 1: Use sample_id only (no strandedness)
        arriba_input_bam_ch = STAR_bam_output
            .map { sample_id, _, bam, bai ->
                tuple(sample_id, bam, bai)
            }
            .view { " arriba_input_bam_ch: ${it}" }

        // Step 2: Run ARRIBA
        ARRIBA(
            arriba_input_bam_ch,
            reference_genome,
            gtf_annotation,
            arriba_blacklist,
            arriba_known_fusions,
			cytobands,
			protein_domains
        )

        ch_versions = ch_versions.mix(ARRIBA.out.versions.first())

        // Step 3: Join for ARRIBA_VISUALIZATION
        fusion_viz_input_ch = arriba_input_bam_ch
            .join(ARRIBA.out.fusions, by: 0)
            

        // Step 4: Visualization
        ARRIBA_VISUALIZATION(
            fusion_viz_input_ch,
            gtf_annotation
        )

        fusion_visual_ch = ARRIBA_VISUALIZATION.out.fusion_plot
        ch_versions = ch_versions.mix(ARRIBA_VISUALIZATION.out.versions.first())

        log.info " Gene Fusion Detection Workflow Completed."

    emit:
        fusion_results        = ARRIBA.out.fusions
        discarded_results     = ARRIBA.out.fusions_discarded
        fusion_visualizations = fusion_visual_ch
        versions              = ch_versions
}