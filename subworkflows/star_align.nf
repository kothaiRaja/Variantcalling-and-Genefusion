// =========================
// SUBWORKFLOW: STAR_ALIGN
// =========================

nextflow.enable.dsl=2

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

    // Define output channels
    star_bam_ch = Channel.empty()
    chimeric_reads_ch = Channel.empty()
    flagstats_ch = Channel.empty()
    align_stats_ch = Channel.empty()
    star_logs_ch = Channel.empty()
    filtered_bams_ch = Channel.empty()
	ch_versions = Channel.empty()
	chimeric_junction_ch = Channel.empty()

    // CASE 1: Perform STAR alignment if no pre-aligned BAMs are provided

        log.info "Running STAR alignment..."
		star_aligned_ch = STAR_ALIGNMENT(trimmed_reads_ch, star_index, gtf_file)
        star_bam_ch      = star_bam_ch.mix(star_aligned_ch.bam)
        chimeric_reads_ch = chimeric_reads_ch.mix(star_aligned_ch.chimeric_sam)
        star_logs_ch      = star_logs_ch.mix(star_aligned_ch.log_final)
		ch_versions = ch_versions.mix(STAR_ALIGNMENT.out.versions.first())
		
		
	

    
	

    // STEP: Sort and index BAMs
    sorted_bams = SAMTOOLS_SORT_INDEX(star_bam_ch)
	sorted_bams_ch = SAMTOOLS_SORT_INDEX.out.bam_sorted 
	ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX.out.versions.first())
	

    // STEP: Collect stats
    align_stats = SAMTOOLS_STATS(sorted_bams_ch)
	align_stats_ch = align_stats_ch.mix(SAMTOOLS_STATS.out.stats)
	ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
	
	

    // STEP: Filter orphans
    filtered_bams = SAMTOOLS_FILTER_ORPHANS(sorted_bams_ch)
	filtered_bams_ch = filtered_bams_ch.mix(SAMTOOLS_FILTER_ORPHANS.out.filtered_sorted_bams)
	ch_versions = ch_versions.mix(SAMTOOLS_FILTER_ORPHANS.out.versions.first())
	

    // STEP: Flagstat
    flagstats_filtered_bam = SAMTOOLS_FLAGSTAT(filtered_bams_ch)
	flagstats_ch = flagstats_ch.mix(SAMTOOLS_FLAGSTAT.out.flagstat)
	ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())
	

    emit:
    bam_sorted     = sorted_bams_ch
    chimeric_reads = chimeric_reads_ch
	chimeric_junction = chimeric_junction_ch
    flagstats      = flagstats_ch
    align_stats    = align_stats_ch
    star_logs      = star_logs_ch
    filtered_bams  = filtered_bams_ch
	versions       = ch_versions
}