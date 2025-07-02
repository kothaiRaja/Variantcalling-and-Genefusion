
nextflow.enable.dsl = 2

// Import variant calling & filtering modules
include { GATK_HAPLOTYPE_CALLER } from '../modules/gatk/haplotypecaller/main.nf'
include { GATK_MERGEVCFS } from '../modules/gatk/merge/main.nf'
include { BCFTOOLS_STATS } from '../modules/bcftools/stats/main.nf'
include { GATK_VARIANT_SELECT_FILTER } from '../modules/gatk/variant_select_filter/main.nf'
include { GATK_VARIANT_FILTER } from '../modules/gatk/variant filter/main.nf'
include { BGZIP_TABIX_ANNOTATIONS as BGZIP_TABIX_VCF  } from '../modules/tabix/bziptabix/main.nf'
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

//    log.info " Starting Variant Calling Workflow..."

ch_haplotypecaller_input = recalibrated_bams
    .combine(intervals) 
    .map { meta, bam, bai, interval_meta, interval_file ->
        def new_meta = meta.clone()
        new_meta.id     = meta.id + "_" + interval_file.baseName
        new_meta.sample = meta.id  
        tuple(new_meta, bam, bai, interval_file)
    }


	ch_known_sites_vcf = known_snps_vcf.collect()
	ch_known_sites_vcf_index = known_snps_vcf_index.collect()
	






	



	
	GATK_variant_caller = GATK_HAPLOTYPE_CALLER(ch_haplotypecaller_input, reference_genome, reference_genome_index, reference_genome_dict, ch_known_sites_vcf, ch_known_sites_vcf_index)
	
	ch_GATK_variant_caller = GATK_HAPLOTYPE_CALLER.out.vcf_output 
	ch_versions = ch_versions.mix(GATK_HAPLOTYPE_CALLER.out.versions.first())


	ch_GATK_variant_caller.view { "HaplotypeCaller_output: $it " }



	ch_merged_vcf_input = GATK_HAPLOTYPE_CALLER.out.vcf_output
    .map { meta, vcf, tbi ->
        def new_meta = [id: meta.sample, sample: meta.sample, strandedness: meta.strandedness]
        tuple(new_meta, vcf, tbi)
    }
    .groupTuple() 








	// Debug Output
	ch_merged_vcf_input.view { "Grouped VCFs for Merging: $it" }


	merged_vcf = GATK_MERGEVCFS(
		ch_merged_vcf_input
	)
	
	merged_vcf_ch = GATK_MERGEVCFS.out.merged_vcf
	ch_versions = ch_versions.mix(GATK_MERGEVCFS.out.versions.first())
	
	// Correct assignment of the new channel
	ch_merged_vcfs = merged_vcf_ch.map { meta, vcf, tbi -> 
		tuple(meta, vcf, tbi)
	}

	// View the output
	ch_merged_vcfs.view { it -> "Merged_vcf_output: $it" }

	

	bcftools_stats = BCFTOOLS_STATS(ch_merged_vcfs)
	
	bcftools_stats_ch = BCFTOOLS_STATS.out.stats_output
	ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())
	


	

	// Choose filtering mode dynamically
if (params.variant_filter_mode == "select") {
    variant_filter = GATK_VARIANT_SELECT_FILTER(ch_merged_vcfs, reference_genome, reference_genome_index, reference_genome_dict)
    ch_variant_filter_uncompressed = GATK_VARIANT_SELECT_FILTER.out.filtered_vcf
    ch_versions = ch_versions.mix(GATK_VARIANT_SELECT_FILTER.out.versions.first())

    // Compress and index (bgzip + tabix)
    compressed_vcf = BGZIP_TABIX_VCF(
        ch_variant_filter_uncompressed.map { sample_id, vcf -> tuple(sample_id, vcf) }
    )
    ch_compressed_vcf = BGZIP_TABIX_VCF.out.compressed_indexed
    ch_versions = ch_versions.mix(BGZIP_TABIX_VCF.out.versions.first())

} else if (params.variant_filter_mode == "global") {
    variant_filter = GATK_VARIANT_FILTER(ch_merged_vcfs, reference_genome, reference_genome_index, reference_genome_dict)
    ch_compressed_vcf = GATK_VARIANT_FILTER.out.filtered_vcf
    ch_versions = ch_versions.mix(GATK_VARIANT_FILTER.out.versions.first())
}


	
        
	ch_compressed_vcf.view { "Filtered VCF (compressed & indexed): $it" }

	

	
	bcftools_query = BCFTOOLS_QUERY(ch_compressed_vcf)
	
	bcftools_query_ch = BCFTOOLS_QUERY.out.query_output
	ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())
	


    emit:
    final_variants = ch_compressed_vcf
	bcftools_stats = bcftools_stats_ch
	bcftools_query = bcftools_query_ch
	versions  = ch_versions
		
	
	
	}
	
	
	
	
