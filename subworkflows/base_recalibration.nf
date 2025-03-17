nextflow.enable.dsl = 2

// Import required processes
include { GATK_BASERECALIBRATOR } from '../modules/gatk/baserecalibration/main.nf'
include { GATK_APPLYBQSR } from '../modules/gatk/applybsqr/main.nf'

workflow BASE_RECALIBRATION {
    take:
    bam_input_ch    // BAMs after merging or CALMD step
    intervals_ch    // Scattered intervals from INTERVAL_PROCESSING
	reference_genome
	reference_genome_index
	reference_genome_dict
	merged_vcf
	merged_vcf_index

    main:
    log.info "Starting Base Recalibration Workflow..."
	
	// Pair split BAMs with Scattered Intervals using combine()

	bam_input_ch.map { tuple(it[0], it[1], it[2], it[3]) }
			.combine(intervals_ch)
			.map { sample_id, strandedness, bam, bai, interval -> 
			tuple(sample_id, strandedness, bam, bai, interval)
       
    }
	.set { ch_recalibrated_bam_bai_interval}
	
		//Run GATK Recalibration 
	recalibrated_bams_table = GATK_BASERECALIBRATOR(ch_recalibrated_bam_bai_interval, reference_genome, reference_genome_index, reference_genome_dict, merged_vcf, merged_vcf_index )
       
	recalibrated_bams_table.view{"Recalibrated Bams Output: $it" }
	
	//Applying BSQR 
	
	// Join Calmd BAMs with Recalibrated Tables (Removing Redundant Strandedness)
ch_applybqsr = bam_input_ch.map { sample_id, strandedness, bam, bai ->  
    tuple(sample_id, strandedness, bam, bai) 
}.join(recalibrated_bams_table.map { sample_id, _ , recal_table ->  // Ignore strandedness
    tuple(sample_id, recal_table)
}, by: 0)  // Join by sample_id

// Debug View to check joined data
ch_applybqsr.view { "Joined output: $it "}

// Combine with Scattered Intervals
ch_applybqsr.combine(intervals_ch)
    .map { sample_id, strandedness, bam, bai, recal_table, interval -> 
        tuple(sample_id, strandedness, bam, bai, recal_table, interval)
    }
    .set { ch_applybqsr_bam_bai_interval }

// Debug View to check final channel
ch_applybqsr_bam_bai_interval.view { "Joined combined channel: $it"}


//Run BSQR
	
	// Step 1: Apply BQSR on scattered BAMs
bams_base_recalibrated = GATK_APPLYBQSR(
    ch_applybqsr_bam_bai_interval, 
    reference_genome, 
    reference_genome_index, 
    reference_genome_dict 
)

bams_base_recalibrated.view { "recalibrated_bsqr_bams:$it " }

emit:
    recalibrated_bams = bams_base_recalibrated
	
}
