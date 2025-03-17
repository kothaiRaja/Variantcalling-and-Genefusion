nextflow.enable.dsl = 2

// Import processes
include { SPLIT_NCIGAR_READS } from '../modules/gatk/splitncigar/main.nf'
include { MERGE_BAMS } from '../modules/samtools/merge/main.nf'
include { SAMTOOLS_CALMD } from '../modules/samtools/calmd/main.nf'


workflow SPLIT_MERGE_BAMS {
    take:
    bam_input_ch   // Input BAMs channel (already processed BAMs)
    intervals_ch   // Scattered interval list
	reference_genome
	reference_genome_index
	reference_genome_dict

    main:
    log.info "Starting Split & Merge BAMs Workflow..."
split_bams_ch = bam_input_ch
    .map { tuple(it[0], it[1], it[2], it[3]) }
    .combine(intervals_ch)
    .map { sample_id, strandedness, bam, bai, interval -> 
        tuple(sample_id, strandedness, bam, bai, interval)
    }
    .set { ch_splitncigar_bam_bai_interval }

// Run SPLIT_NCIGAR_READS process using the paired inputs
split_bams = SPLIT_NCIGAR_READS(
    ch_splitncigar_bam_bai_interval,
    reference_genome, 
    reference_genome_index, 
    reference_genome_dict
)

split_bams.view { "SplitNcigar Output : $it" }

split_bams
    .groupTuple() // Groups by sample ID
    .map { sample_id, strandedness_list, bams, bais -> 
        tuple(sample_id, strandedness_list.unique()[0], bams.flatten(), bais.flatten()) 
    }
    .set { ch_merged_bams }


merged_bams = MERGE_BAMS(ch_merged_bams)

merged_bams.view { "Merged BAMs Output : $it" }

// **Step 3: Apply samtools calmd**
    calmd_bams_ch = SAMTOOLS_CALMD(merged_bams, reference_genome, reference_genome_index)

	calmd_bams_ch.view { "Calmd BAM Output: $it" }


// **Emit Outputs**
    emit:
    split_bams  = split_bams
    merged_bams = merged_bams
	merged_calmd_bams 	= calmd_bams_ch

}