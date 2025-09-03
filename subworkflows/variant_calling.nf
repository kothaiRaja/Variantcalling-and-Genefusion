nextflow.enable.dsl = 2

// Import variant calling & filtering modules
include { GATK_HAPLOTYPE_CALLER } from '../modules/gatk/haplotypecaller/main.nf'
include { GATK_MERGEVCFS       } from '../modules/gatk/merge/main.nf'
include { BCFTOOLS_STATS       } from '../modules/bcftools/stats/main.nf'
include { GATK_VARIANT_SELECT_FILTER } from '../modules/gatk/variant_select_filter/main.nf'
include { GATK_VARIANT_FILTER        } from '../modules/gatk/variant_filter/main.nf'
include { BGZIP_TABIX_ANNOTATIONS as BGZIP_TABIX_VCF } from '../modules/tabix/bziptabix/main.nf'
include { BCFTOOLS_QUERY } from '../modules/bcftools/query/main.nf'

workflow VARIANT_CALLING {
    take:
    recalibrated_bams                   
    intervals                           
    reference_genome
    reference_genome_index
    reference_genome_dict
    known_snps_vcf
    known_snps_vcf_index

    main:
    ch_versions = Channel.empty()
	
	// Pairing recalibrated BAMs with intervals 
	ch_haplotypecaller_input = recalibrated_bams
	.combine(intervals)
	.map { t ->
      def m   = t[0]
      def bam = t[1]
      def bai = t[2]
      def iv  = t[-1]    

      def mm = m.clone()
      mm.sample = mm.sample ?: mm.id
      mm.id     = mm.sample          
      mm.shard  = iv.baseName        

      tuple(mm, bam, bai, iv)
  }


    // Collect known sites 
    ch_known_sites_vcf       = known_snps_vcf.collect()
    ch_known_sites_vcf_index = known_snps_vcf_index.collect()

    // HaplotypeCaller per (sample, interval)
    GATK_variant_caller = GATK_HAPLOTYPE_CALLER(
        ch_haplotypecaller_input,
        reference_genome,
        reference_genome_index,
        reference_genome_dict,
        ch_known_sites_vcf,
        ch_known_sites_vcf_index
    )
    ch_versions = ch_versions.mix(GATK_HAPLOTYPE_CALLER.out.versions.first())

    // Group shards back per sample for merging
    ch_merged_vcf_input = GATK_HAPLOTYPE_CALLER.out.vcf_output
        .map { meta, vcf, tbi ->
            tuple([id: meta.sample, sample: meta.sample, strandedness: meta.strandedness], vcf, tbi)
        }
        .groupTuple()

    // Merge per-sample VCFs
    merged_vcf = GATK_MERGEVCFS(ch_merged_vcf_input)
    ch_versions = ch_versions.mix(GATK_MERGEVCFS.out.versions.first())

    // Normalize to (meta, vcf, tbi) for downstream
    ch_merged_vcfs = GATK_MERGEVCFS.out.merged_vcf.map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }

    // Stats on merged VCF
    BCFTOOLS_STATS(ch_merged_vcfs)
    bcftools_stats_ch = BCFTOOLS_STATS.out.stats_output
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    // Filter strategy
    if (params.variant_filter_mode == 'select') {
        GATK_VARIANT_SELECT_FILTER(ch_merged_vcfs, reference_genome, reference_genome_index, reference_genome_dict)
        ch_versions = ch_versions.mix(GATK_VARIANT_SELECT_FILTER.out.versions.first())

        // compress & index if module emits uncompressed VCF
        BGZIP_TABIX_VCF(
            GATK_VARIANT_SELECT_FILTER.out.filtered_vcf.map { meta, vcf -> tuple(meta, vcf) }
        )
        ch_versions = ch_versions.mix(BGZIP_TABIX_VCF.out.versions.first())
        ch_compressed_vcf = BGZIP_TABIX_VCF.out.compressed_indexed

    } else if (params.variant_filter_mode == 'global') {
        GATK_VARIANT_FILTER(ch_merged_vcfs, reference_genome, reference_genome_index, reference_genome_dict)
        ch_versions = ch_versions.mix(GATK_VARIANT_FILTER.out.versions.first())
        ch_compressed_vcf = GATK_VARIANT_FILTER.out.filtered_vcf
    } else {
        ch_compressed_vcf = ch_merged_vcfs
    }

    // Bcftools query summary on final VCFs
    BCFTOOLS_QUERY(ch_compressed_vcf)
    bcftools_query_ch = BCFTOOLS_QUERY.out.query_output
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

    emit:
    final_variants  = ch_compressed_vcf
    bcftools_stats  = bcftools_stats_ch
    bcftools_query  = bcftools_query_ch
    versions        = ch_versions
}
