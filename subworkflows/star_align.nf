nextflow.enable.dsl=2

include { STAR_ALIGNMENT }          from '../modules/star/main.nf'
include { SAMTOOLS_SORT_INDEX }     from '../modules/samtools/sort/main.nf'
include { SAMTOOLS_STATS }          from '../modules/samtools/stats/main.nf'
include { SAMTOOLS_FILTER_ORPHANS } from '../modules/samtools/filter_orphans/main.nf'
include { SAMTOOLS_FLAGSTAT }       from '../modules/samtools/flagstat/main.nf'

workflow STAR_ALIGN {
    take:
    trimmed_reads
    star_index
    gtf_file

    main:
    def ch_versions = Channel.empty()
	



    // STAR alignment
    star_output = STAR_ALIGNMENT(trimmed_reads, star_index, gtf_file)

    // Outputs from STAR
    star_bam_ch          = star_output.bam
    chimeric_reads_ch    = star_output.chimeric_sam
    chimeric_junction_ch = star_output.chimeric_junction
    log_final_ch   = star_output.log_final
    ch_versions          = ch_versions.mix(star_output.versions)

    // Sort BAM
    SAMTOOLS_SORT_INDEX(star_bam_ch)
    sorted_bams_ch = SAMTOOLS_SORT_INDEX.out.bam_sorted
    ch_versions    = ch_versions.mix(SAMTOOLS_SORT_INDEX.out.versions)

    // Alignment Stats
    SAMTOOLS_STATS(sorted_bams_ch)
    stats_ch       = SAMTOOLS_STATS.out.stats
    ch_versions    = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    // Filter Orphans
    SAMTOOLS_FILTER_ORPHANS(sorted_bams_ch)
    filtered_bams_ch = SAMTOOLS_FILTER_ORPHANS.out.filtered_sorted_bams
    ch_versions      = ch_versions.mix(SAMTOOLS_FILTER_ORPHANS.out.versions)

    // Flagstat
    SAMTOOLS_FLAGSTAT(filtered_bams_ch)
    flagstat_ch    = SAMTOOLS_FLAGSTAT.out.flagstat
    ch_versions  = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)
	
	


    emit:
    bam_sorted         = sorted_bams_ch
    chimeric_reads     = chimeric_reads_ch
    chimeric_junction  = chimeric_junction_ch
     flagstats      = flagstat_ch
    align_stats    = stats_ch
    star_logs      = log_final_ch
    filtered_bams  = filtered_bams_ch
	versions       = ch_versions

    
}
