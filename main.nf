nextflow.enable.dsl = 2

// Include the QC Subworkflow
include { STAR_ALIGNMENT } from './modules/STAR_mapping/star_varaint.nf'
include { SAMTOOLS_SORT_INDEX } from './modules/SAMTOOLS_sort_index/samtools_sort_index.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/SAMTOOLS_flagstat/samtools_stat.nf'
include { GATK_MARK_DUPLICATES } from './modules/GATK_mark_duplicates/mark_duplicates.nf'
include { SPLIT_NCIGAR_READS } from './modules/Split_NCIGAR_reads/splitNCIGAR.nf'
include { GATK_RECALIBRATION } from './modules/GATK_recalibration/recalibration.nf'
include { BED_TO_INTERVAL_LIST } from './modules/BED_to_INTERVAL_list/Interval_list.nf'
include { SCATTER_INTERVAL_LIST  } from './modules/Scatter_interval_list/scatter_intervals.nf'
include { GATK_HAPLOTYPE_CALLER  } from './modules/GATK_Haplotype_caller/variant_call.nf'
include {GATK_MERGE_VCFS  } from './modules/GATK_merge_VCFs/merge_vcfs.nf'
include { GATK_VARIANT_FILTER  } from './modules/GATK_varaint_filter/variant_filter.nf'
include { ANNOTATE_VARIANTS  } from './modules/SnpEFF_Annotations/annotations.nf'




// Include statements for ARRIBA fusion
//include { STAR_ALIGN_FUSION } from './modules/STAR_mapping/star_fusion.nf'
//include { ARRIBA } from './modules/ARRIBA/arriba_fusion.nf'







workflow {
    
//PartA: Variant calling
    // Load trimmed reads
    trimmed_reads_ch = Channel.fromFilePairs(params.reads, flat: true)
        
    // Run STAR Alignment
    aligned_bams = STAR_ALIGNMENT(trimmed_reads_ch, params.star_index_dir)
	
	// Sort and index BAM files
    sorted_bams = SAMTOOLS_SORT_INDEX(aligned_bams)
	
	// Generate alignment statistics
    alignment_stats = SAMTOOLS_FLAGSTAT(sorted_bams)
	
	// Mark duplicates
    marked_bams = GATK_MARK_DUPLICATES(sorted_bams)
	
	// Split N CIGAR reads
    split_bams = SPLIT_NCIGAR_READS(marked_bams.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict) })
	
	// Recalibrate and Apply BQSR in one step
    recalibrated_bams = GATK_RECALIBRATION(
        split_bams.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict, params.filtered_vcf, params.filtered_vcf_index) })
		
	// Convert BED to interval list
    interval_list_ch = BED_TO_INTERVAL_LIST(params.denylist, params.genome, params.genome_dict)
	
	// Scatter the Interval List
    scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, params.genome_dict)
	
	//GATK HaplotypeCaller
	gvcf_output = GATK_HAPLOTYPE_CALLER(recalibrated_bams, params.genome,params.fasta_index,params.genome_dict,  scattered_intervals_ch)
	

	// Combine GVCFs
	merged_vcf = GATK_MERGE_VCFS(params.genome,params.fasta_index, params.genome_dict,gvcf_output)

	//// Variant Filtering
    filtered_vcf = GATK_VARIANT_FILTER(merged_vcf,params.genome,params.fasta_index, params.genome_dict )
	
	// Pass all required inputs to ANNOTATE_VARIANTS
    annotated_vcf = ANNOTATE_VARIANTS(filtered_vcf, file('./data/test/snpEff/snpEff.jar'),
        file('./data/test/snpEff/snpEff.config'),
        file('./data/test/snpEff/data') , params.genomedb) 

//PART B:

	// PART 2: STAR RNA-Seq Mapping
      //star_align_ch = STAR_ALIGN_FUSION( trimmed_reads_ch,params.star_index_dir, params.gtf_file)
	  
	  
	  
	  //Part 3: ARRIBA
	  //ARRIBA_ch = ARRIBA(star_align_ch, params.genome, params.gtf_file, params.denylist, params.known_fusions)

}