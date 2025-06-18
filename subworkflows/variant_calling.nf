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
    recalibrated_bams_ch  
	intervals_ch 
    reference_genome
	reference_genome_index
	reference_genome_dict
	known_variants
	known_variants_index
	
	main:
	
	ch_versions = Channel.empty()

    log.info " Starting Variant Calling Workflow..."
	
	recalibrated_bams_ch
	.map { tuple(it[0], it[1], it[2], it[3]) }
    .combine(intervals_ch)
    .map { sample_id, strandedness, bam, bai, interval -> 
        tuple(sample_id, strandedness, bam, bai, interval)
    }
    .set { ch_haplotypecaller_input }

	
	GATK_variant_caller = GATK_HAPLOTYPE_CALLER(ch_haplotypecaller_input,reference_genome, reference_genome_index, reference_genome_dict, known_variants, known_variants_index)
	
	ch_GATK_variant_caller = GATK_HAPLOTYPE_CALLER.out.vcf_output 
	ch_versions = ch_versions.mix(GATK_HAPLOTYPE_CALLER.out.versions.first())


	ch_GATK_variant_caller.view { "HaplotypeCaller_output: $it " }



	ch_GATK_variant_caller
		.groupTuple() // Groups by sample ID
		.map { sample_id, strandedness_list, vcf, tbi -> 
			tuple(sample_id, strandedness_list.unique()[0], vcf.flatten(), tbi.flatten()) 
    }
		.set { ch_haplotypecaller_grouped }

	// Debug Output
	ch_haplotypecaller_grouped.view { "Grouped VCFs for Merging: $it" }


	merged_vcf = GATK_MERGEVCFS(
		ch_haplotypecaller_grouped.map { sample_id, strandedness, vcf_list, tbi_list -> 
			tuple(sample_id, vcf_list, tbi_list)
		}
	)
	
	merged_vcf_ch = GATK_MERGEVCFS.out.merged_vcf
	ch_versions = ch_versions.mix(GATK_MERGEVCFS.out.versions.first())
	

	bcftools_stats = BCFTOOLS_STATS(merged_vcf_ch)
	
	bcftools_stats_ch = BCFTOOLS_STATS.out.stats_output
	ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())
	


	// Correct assignment of the new channel
	ch_merged_vcfs = merged_vcf_ch.map { sample_id, vcf, tbi -> 
		tuple(sample_id, vcf, tbi)
	}

	// View the output
	ch_merged_vcfs.view { it -> "Merged_vcf_output: $it" }


	// Choose filtering mode dynamically
if (params.variant_filter_mode == "select") {
    variant_filter = GATK_VARIANT_SELECT_FILTER(ch_merged_vcfs, reference_genome, reference_genome_index, reference_genome_dict)
    ch_variant_filter_uncompressed = GATK_VARIANT_SELECT_FILTER.out.filtered_vcf
    ch_versions = ch_versions.mix(GATK_VARIANT_SELECT_FILTER.out.versions.first())

    // Compress and index (bgzip + tabix)
    compressed_vcf = BGZIP_TABIX_VCF(
        ch_variant_filter_uncompressed.map { sample_id, vcf -> tuple(sample_id, vcf) }
    )
    ch_compressed_vcf = compressed_vcf.out.compressed_indexed
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
	
	
	
	
