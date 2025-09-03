nextflow.enable.dsl = 2

// Import processes
include { SPLIT_NCIGAR_READS } from '../modules/gatk/splitncigar/main.nf'
include { MERGE_BAMS }         from '../modules/samtools/merge/main.nf'
include { RESET_READGROUPS }   from '../modules/gatk/AddRG_Group/main.nf'
include { SAMTOOLS_CALMD }     from '../modules/samtools/calmd/main.nf'

workflow SPLIT_MERGE_BAMS {
    take:
    bam_input_ch        
    intervals_ch       
    reference_genome
    reference_genome_index
    reference_genome_dict

    main:

    ch_versions = Channel.empty()
//    log.info " Starting Split & Merge BAMs Workflow..."

    // STEP 1: Combine BAMs with intervals (scatter)
   
bam_input_ch
    .combine(intervals_ch)             
    .map { meta, bam, bai, interval ->
        def m = meta.clone()
        m.sample = meta.sample ?: meta.id     
        m.id     = m.sample                   
        m.shard  = interval.baseName          
        tuple(m, bam, bai, interval)
    }
    .set { ch_splitncigar_bam_bai_interval }


//    ch_splitncigar_bam_bai_interval.view { " SPLIT INPUT: $it" }

    // STEP 2: Run SplitNCigarReads per interval
    split_bams = SPLIT_NCIGAR_READS(
        ch_splitncigar_bam_bai_interval,
        reference_genome,
        reference_genome_index,
        reference_genome_dict
    )

    split_bams_ch = SPLIT_NCIGAR_READS.out.split_interval_bams
    ch_versions = ch_versions.mix(SPLIT_NCIGAR_READS.out.versions)

//    split_bams_ch.view { " Split BAM: $it" }

    split_bams_ch
    .map { m, bam, bai -> 
        tuple(m.sample, [m, bam, bai]) 
    }
    .groupTuple()
  .map { sample_id, entries ->
      def m = entries[0][0].clone()
      m.id = m.sample        
      def bams = entries.collect { it[1] }
      def bais = entries.collect { it[2] }
      tuple(m, bams, bais)
  }
  .set { ch_merged_bams }

 ch_merged_bams.view { " Merged Input: $it" }
  

    // STEP 4: Merge BAMs
    merged_bams = MERGE_BAMS(ch_merged_bams)
    merged_bams_ch = MERGE_BAMS.out.merged_bams
    ch_versions = ch_versions.mix(MERGE_BAMS.out.versions)

//    merged_bams_ch.view { " Merged BAM: $it" }

    // STEP 5: Reset Read Groups
    reset_bams = RESET_READGROUPS(merged_bams_ch)
    reset_rg_bams_ch = RESET_READGROUPS.out.fixed_bams
    ch_versions = ch_versions.mix(RESET_READGROUPS.out.versions)

//    reset_rg_bams_ch.view { " Fixed RG BAM: $it" }

    // STEP 6: Apply samtools calmd
    calmd_bams = SAMTOOLS_CALMD(reset_rg_bams_ch, reference_genome, reference_genome_index)
    calmd_bams_ch = SAMTOOLS_CALMD.out.calmd_bams
    ch_versions = ch_versions.mix(SAMTOOLS_CALMD.out.versions)

    // OUTPUT
    emit:
    split_bams          = split_bams_ch
    merged_bams         = merged_bams_ch
    fixed_bams          = reset_rg_bams_ch
    merged_calmd_bams   = calmd_bams_ch
    versions            = ch_versions
}
