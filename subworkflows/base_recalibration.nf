nextflow.enable.dsl = 2

include { GATK_BASERECALIBRATOR } from '../modules/gatk/baserecalibration/main.nf'
include { GATK_APPLYBQSR       } from '../modules/gatk/applybsqr/main.nf'

workflow BASE_RECALIBRATION {
  take:
    bam_input
    recalib_ch
    intervalsch
    reference_genome
    reference_genome_index
    reference_genome_dict
    known_snps_vcf
    known_snps_index
    known_indels_vcf
    known_indels_index

  main:
    ch_versions = Channel.empty()

    // 1) Known sites as singletons
    ch_known_sites_vcf   = known_snps_vcf.concat(known_indels_vcf).collect()
    ch_known_sites_index = known_snps_index.concat(known_indels_index).collect()

    // 2) BaseRecalibrator once per sample 
    baserecalibrator_results = GATK_BASERECALIBRATOR(
      recalib_ch,
      reference_genome, reference_genome_index, reference_genome_dict,
      ch_known_sites_vcf, ch_known_sites_index
    )
    ch_versions = ch_versions.mix(baserecalibrator_results.versions)

    // 3) JOIN BY SAMPLE KEY 
    bams_by_sample = bam_input.map { m, bam, bai ->
      def sample = (m.sample ?: m.id)
      tuple(sample, m, bam, bai)
    }
    tables_by_sample = baserecalibrator_results.recal_table.map { m, table ->
      def sample = (m.sample ?: m.id)
      tuple(sample, table)
    }
    applybqsr_joined_ch = bams_by_sample.join(tables_by_sample)
      .map { sample, m, bam, bai, table ->
        def mm = m.clone()
        mm.sample = sample
        mm.id     = sample
        tuple(mm, bam, bai, table)
      }

    // 4) Pair with intervals 
	interval_only = intervalsch.map { iv_meta, iv_path -> iv_path }

	// pair BAM+table with the single interval path
	ch_bqsr_apply_input = applybqsr_joined_ch
	.combine(interval_only)
	.map { s, bam, bai, table, iv ->
      tuple([id: s.sample, sample: s.sample], bam, bai, table, iv)
  }


    // 5) ApplyBQSR
    GATK_APPLYBQSR(
      ch_bqsr_apply_input,
      reference_genome, reference_genome_index, reference_genome_dict
    )
    ch_versions = ch_versions.mix(GATK_APPLYBQSR.out.versions)

    emit:
      recalibrated_bams = GATK_APPLYBQSR.out.recalibrated_bam
      versions          = ch_versions
}
