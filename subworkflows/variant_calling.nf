nextflow.enable.dsl = 2

// Import variant calling & filtering modules
include { GATK_HAPLOTYPE_CALLER } from '../modules/gatk/haplotypecaller/main.nf'
include { GATK_MERGEVCFS } from '../modules/gatk/merge/main.nf'
include { BCFTOOLS_STATS } from '../modules/bcftools/stats/main.nf'
include { GATK_VARIANT_FILTER } from '../modules/gatk/variant filter/main.nf'
include { BCFTOOLS_QUERY } from '../modules/bcftools/query/main.nf'
include { SELECT_SNPs } from '../modules/gatk/select_variants/select_snps/main.nf'
include { SELECT_INDELs } from '../modules/gatk/select_variants/select_indels/main.nf'

workflow VARIANT_CALLING {
    take:
    recalibrated_bams_ch  // Recalibrated BAM files from ApplyBQSR
    reference_genome
	reference_genome_index
	reference_genome_dict
	known_variants
	known_variants_index
	
	main:
    log.info "🔹 Starting Variant Calling Workflow..."
	
	ch_GATK_variant_caller = GATK_HAPLOTYPE_CALLER(recalibrated_bams_ch,reference_genome, reference_genome_index, reference_genome_dict, known_variants, known_variants_index)

	ch_GATK_variant_caller.view { "HaplotypeCaller_output: $it " }



	ch_GATK_variant_caller
		.groupTuple() // Groups by sample ID
		.map { sample_id, strandedness_list, vcf, tbi -> 
			tuple(sample_id, strandedness_list.unique()[0], vcf.flatten(), tbi.flatten()) 
    }
		.set { ch_haplotypecaller_grouped }

	// Debug Output
	ch_haplotypecaller_grouped.view { "Grouped VCFs for Merging: $it" }


	merged_vcf_ch = GATK_MERGEVCFS(
		ch_haplotypecaller_grouped.map { sample_id, strandedness, vcf_list, tbi_list -> 
			tuple(sample_id, vcf_list, tbi_list)
		}
	)

	bcftools_stats_ch = BCFTOOLS_STATS(merged_vcf_ch)


	// Correct assignment of the new channel
	ch_merged_vcfs = merged_vcf_ch.map { sample_id, vcf, tbi -> 
		tuple(sample_id, vcf, tbi)
	}

	// View the output
	ch_merged_vcfs.view { it -> "Merged_vcf_output: $it" }


	ch_variant_filter = GATK_VARIANT_FILTER(ch_merged_vcfs, reference_genome, reference_genome_index, reference_genome_dict)
        
	ch_variant_filter.view {"Filtered_vcf: $it"}

	bcftools_query_ch = BCFTOOLS_QUERY(ch_variant_filter)


	if (params.select_variants) {
		selected_snps_ch = SELECT_SNPs(ch_variant_filter, reference_genome, reference_genome_index, reference_genome_dict)
		selected_snps_ch.view {"SNPs output: $it"}
		selected_indels_ch = SELECT_INDELs(ch_variant_filter, reference_genome, reference_genome_index, reference_genome_dict)
		selected_indels_ch.view {"INDELS output: $it"}
		
	 // Combine SNPs & INDELs before annotation
        ch_selected_variants = selected_snps_ch.join(selected_indels_ch, by: 0)
        
    }

    emit:
    final_variants = ch_variant_filter
    selected_snps = selected_snps_ch ?: Channel.empty()
    selected_indels = selected_indels_ch ?: Channel.empty()
    selected_variants = ch_selected_variants
		
	
	
	}
	
	
	
	
