nextflow.enable.dsl = 2

// Import necessary processes
include { SAMTOOLS_FILTER_ORPHANS } from '../modules/samtools/filter_orphans/main.nf'
include { GATK_MARK_DUPLICATES } from '../modules/gatk/mark_duplicates/main.nf'
include { SAMTOOLS_SORT_INDEX } from '../modules/samtools/sort/main.nf'
include { SAMTOOLS_STATS } from '../modules/samtools/stats/main.nf'
include { SAMTOOLS_FLAGSTAT } from '../modules/samtools/flagstat/main.nf'

def get_sample_metadata(bam_path) {
    def sample_id = bam_path.baseName.replace('.bam', '') // Extract sample name from BAM file

    // Check for strandedness from params.strandedness_map
    def strandedness = params.strandedness_map.containsKey(sample_id) ? 
                        params.strandedness_map[sample_id] : 
                        null

    // **If not found in map, check in sample sheet**
    if (!strandedness && params.samplesheet) {
        def sample_data = file(params.samplesheet)
            .readLines()
            .drop(1) // Skip header
            .collectEntries { line ->
                def fields = line.split(',') // Adjust if CSV uses different delimiter (e.g., `split('\t')` for TSV)
                [(fields[0].trim()): fields[3].trim()] // Extract `sample_id` (col 1) -> `strandedness` (col 4)
            }

        strandedness = sample_data.get(sample_id, 'unstranded') // Default to "unstranded"
    }

    return [sample_id, bam_path, strandedness]
}

workflow MARK_DUPLICATES {
    take:
    bam_input_ch   

    main:
    log.info "Starting MarkDuplicates Workflow..."

    // **Step 1: Determine Input BAMs**
    filtered_bams_ch = Channel.empty()

    if (params.input_bam) {
        log.info "Using user-provided BAM files for MarkDuplicates..."
        bam_channel = Channel.fromPath(params.input_bam)
					.map { bam_file -> get_sample_metadata(bam_file) }
		indexed_bams_ch = SAMTOOLS_SORT_INDEX(bam_channel)
		alignment_stats_ch = SAMTOOLS_STATS(indexed_bams_ch)
        filtered_bams_ch = SAMTOOLS_FILTER_ORPHANS(bam_channel)
		flagstat_ch = SAMTOOLS_FLAGSTAT(filtered_bams_ch)
    } else {
        log.info "Using filtered BAMs from STAR alignment..."
        filtered_bams_ch = bam_input_ch 
    }

    // **Step 2: Mark Duplicates**
    dedup_bams_ch = GATK_MARK_DUPLICATES(filtered_bams_ch)

    // **Emit Outputs**
    emit:
    marked_bams_bai = dedup_bams_ch.marked_bams_bai
	marked_bams_bai_metrics = dedup_bams_ch.marked_bams_bai_metrics
}
