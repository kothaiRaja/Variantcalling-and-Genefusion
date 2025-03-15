nextflow.enable.dsl = 2

// Import processes
include { STAR_ALIGNMENT } from '../modules/star/main.nf'
include { SAMTOOLS_SORT_INDEX } from '../modules/samtools/sort/main.nf'
include { SAMTOOLS_STATS } from '../modules/samtools/stats/main.nf'
include { SAMTOOLS_FILTER_ORPHANS } from '../modules/samtools/filter_orphans/main.nf'
include { SAMTOOLS_FLAGSTAT } from '../modules/samtools/flagstat/main.nf'

workflow STAR_ALIGN {
    take:
    trimmed_reads_ch   
    star_index      
    gtf_file        

    main:
    
    // Step 1: STAR Alignment
    star_aligned_ch = STAR_ALIGNMENT(trimmed_reads_ch, star_index, gtf_file)


    // Step 2: Sort and Index BAM Files
    indexed_bams_ch = SAMTOOLS_SORT_INDEX(star_aligned_ch[0])
	indexed_bams_ch.view { "sorted_indexed BAMs Output: $it" }
	
	//Step 3: samtools stats
	alignment_stats_ch = SAMTOOLS_STATS(indexed_bams_ch)
		
		
	// Step 4: Remove orphan reads
    filtered_bams_ch = SAMTOOLS_FILTER_ORPHANS(indexed_bams_ch)
	filtered_bams_ch.view { "Filtered BAMs Output: $it" }

    // Step 5: Generate alignment statistics
    flagstat_ch = SAMTOOLS_FLAGSTAT(filtered_bams_ch)
	
	

    emit:
    bam_sorted   = indexed_bams_ch
    chimeric_reads = star_aligned_ch.chimeric_sam // Emit Chimeric Reads for Fusion Detection
    flagstats    = flagstat_ch 	 // Emit flagstat output
	align_stats  = alignment_stats_ch
    star_logs       = star_aligned_ch.log_final
	star_log_out = star_aligned_ch.log_out
}
