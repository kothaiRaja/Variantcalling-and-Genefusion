nextflow.enable.dsl = 2

// Import processes related to intervals
include { BED_TO_INTERVAL_LIST } from '../modules/Intervals/bed_to_intervals/main.nf'
include { SCATTER_INTERVAL_LIST } from '../modules/Intervals/scattered_intervals/main.nf'

workflow INTERVAL_PROCESSING {
    take:
    bed_file_ch
	reference_genome
	reference_genome_dict
	
	 main:
	 
	ch_versions = Channel.empty()
    log.info "Starting Interval Processing..."


	// **Step 1: Convert BED to Interval List**
		interval_list = BED_TO_INTERVAL_LIST(bed_file_ch, reference_genome, reference_genome_dict)
		
		interval_list_ch = BED_TO_INTERVAL_LIST.out.interval_list
		ch_versions = ch_versions.mix(BED_TO_INTERVAL_LIST.out.versions.first())
		
		

	// Step 2: Scatter if enabled, else use the full list
    scattered_intervals_ch = Channel.empty()

    if (params.scatterintervals) {
        log.info " Scattering intervals for parallel execution..." 
		scattered_intervals = SCATTER_INTERVAL_LIST(interval_list_ch, params.reference_genome_dict)
		scattered_intervals_ch = SCATTER_INTERVAL_LIST.out.scattered_intervals
								.map { meta, file -> file }  
								.flatten()
        
        scattered_intervals_ch.view { file -> " Scattered interval: $file" }
        ch_versions = ch_versions.mix(SCATTER_INTERVAL_LIST.out.versions.first())
    } else {
        log.info " Using full interval list without scattering..."
        scattered_intervals_ch = interval_list_ch
    }

	
	emit:
    intervals = scattered_intervals_ch 
	versions  = ch_versions
}