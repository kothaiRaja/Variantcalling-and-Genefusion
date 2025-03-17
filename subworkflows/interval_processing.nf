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
    log.info "Starting Interval Processing..."


	// **Step 1: Convert BED to Interval List**
		interval_list_ch = BED_TO_INTERVAL_LIST(bed_file_ch, reference_genome, reference_genome_dict)
		
		

	// **Step 2: Check if Scatter Intervals is Enabled**	
		scattered_intervals_ch = Channel.empty() 

		if (params.scatterintervals) {
    // Scatter the interval list
    scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, reference_genome_dict)
        .map { meta, file -> file }  
        .flatten()
    
    println " Using **scattered intervals** for parallel processing"
    
    // **VIEW statement to debug output**
    scattered_intervals_ch.view { file -> "Generated Scattered Interval: ${file}" }

	} else {
    // Use the full interval list directly
    scattered_intervals_ch = interval_list_ch.map { meta, file -> file }  
    println " Using **full exons interval list** without scattering"
	}
	
	emit:
    intervals = scattered_intervals_ch 
}