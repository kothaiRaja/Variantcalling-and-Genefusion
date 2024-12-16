nextflow.enable.dsl = 2

// Include the QC Subworkflow
include { STAR_ALIGNMENT } from './modules/STAR_mapping/star_varaint.nf'
include { SAMTOOLS_SORT_INDEX } from './modules/SAMTOOLS/samtools_sort_index.nf'
include { SAMTOOLS_FILTER_ORPHANS } from './modules/SAMTOOLS/filter_orphans.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/SAMTOOLS/samtools_stat.nf'
include { GATK_MARK_DUPLICATES } from './modules/GATK/mark_duplicates.nf'
include { SPLIT_NCIGAR_READS } from './modules/GATK/splitNCIGAR.nf'
include { SAMTOOLS_CALMD } from './modules/SAMTOOLS/call_md.nf'
include { GATK_RECALIBRATION } from './modules/GATK/recalibration.nf'
include { BED_TO_INTERVAL_LIST } from './modules/PREPARE_INTERVALS/Interval_list.nf'
include { SCATTER_INTERVAL_LIST  } from './modules/PREPARE_INTERVALS/scatter_intervals.nf'
include { GATK_HAPLOTYPE_CALLER  } from './modules/GATK/variant_call.nf'
include {BCFTOOLS_MERGE  } from './modules/BCFTOOLS/merge.nf'
include {BCFTOOLS_STATS  } from './modules/BCFTOOLS/bcf_stats.nf'
include { GATK_VARIANT_FILTER  } from './modules/GATK/variant_filter.nf'
include { ANNOTATE_VARIANTS  } from './modules/SnpEFF_ANNOTATIONS/annotations.nf'




// Include statements for ARRIBA fusion
include { STAR_ALIGN_FUSION } from './modules/STAR_mapping/star_fusion.nf'
include { ARRIBA } from './modules/ARRIBA/arriba_fusion.nf'
include { ARRIBA_VISUALIZATION  } from './modules/ARRIBA/arriba_visualisation.nf'

//Include statement for multiqc
include { MULTIQC_REPORT } from './modules/QUALITY_CONTROL/multiqc.nf'


workflow {
//==============PART A: Variant_calling===============//    

//Step 1: Load trimmed reads
    trimmed_reads_ch = Channel.fromFilePairs(params.reads, flat: true)
        
//Step 2: Run STAR Alignment
    aligned_bams = STAR_ALIGNMENT(trimmed_reads_ch, params.star_index_dir)
	
//Step 3: Sort and index BAM files
    sorted_bams = SAMTOOLS_SORT_INDEX(aligned_bams)
	
//Step 5: Filter the Orphan reads
	filtered_bams = SAMTOOLS_FILTER_ORPHANS(sorted_bams)
	
//Step 6: Generate alignment statistics
    alignment_stats = SAMTOOLS_FLAGSTAT(filtered_bams)
	
//Step 7: Mark duplicates
    marked_bams = GATK_MARK_DUPLICATES(filtered_bams)
	
//Step 8: Split N CIGAR reads
    split_bams = SPLIT_NCIGAR_READS(marked_bams.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict) })
	
//Step 9: SAMTOOLS CALMD process
	calmd_ch = SAMTOOLS_CALMD(split_bams , params.genome, params.fasta_index )
	
//Step 10: Recalibrate and Apply BQSR in one step
    recalibrated_bams = GATK_RECALIBRATION(
        calmd_ch.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict, params.filtered_vcf, params.filtered_vcf_index) })
		
//Step 11: Convert BED to interval list
    interval_list_ch = BED_TO_INTERVAL_LIST(params.denylist, params.genome, params.genome_dict)
	
//Step 12: Scatter the Interval List
    scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, params.genome_dict)
	
//Step 13: GATK HaplotypeCaller
	gvcf_output = GATK_HAPLOTYPE_CALLER(recalibrated_bams, params.genome, params.fasta_index, params.genome_dict, scattered_intervals_ch)
				  .map { it[1] } 
				  .collect()
	
	gvcf_output.view { "Raw GVCF output: $it" }


//Step 14: Merge VCF files
        merged_vcf_ch = BCFTOOLS_MERGE(gvcf_output)

	
//Step 15: provide stats
	bcftools_stats_ch = BCFTOOLS_STATS(merged_vcf_ch)

//Step 16: Variant Filtering
    filtered_vcf = GATK_VARIANT_FILTER(merged_vcf_ch,params.genome,params.fasta_index, params.genome_dict )
	
//Step 17: Dynamically set paths based on the mode (test or actual)
    def snpEffJar = file(params.mode == 'test' ? './data/test/snpEff/snpEff.jar' : './data/actual/snpEff/snpEff.jar')
    def snpEffConfig = file(params.mode == 'test' ? './data/test/snpEff/snpEff.config' : './data/actual/snpEff/snpEff.config')
    def snpEffDbDir = file(params.mode == 'test' ? './data/test/snpEff/snpEff/data' : './data/actual/snpEff/data')

//Step 18:  ANNOTATE_VARIANTS
    annotated_vcf = ANNOTATE_VARIANTS(filtered_vcf, snpEffJar, snpEffConfig, snpEffDbDir, params.genomedb)
	
//===============Workflow for ARRIBA===================//

//Step 19: STAR alignment for fusion detection
	star_align_fusion_ch = STAR_ALIGN_FUSION(trimmed_reads_ch, params.star_index_dir, params.gtf_file )

//Step 20: Fusion detection using ARRIBA
	ARRIBA_ch = ARRIBA(star_align_fusion_ch, params.genome, params.gtf_file, params.test_blacklist_fusion, params.test_knownfusion)
	
//Step 21: Visualization step
    fusion_visuals = ARRIBA_VISUALIZATION(ARRIBA_ch, params.scripts_dir, params.genome, params.gtf_file)
	
	
collect_results_ch = Channel.fromPath("${params.outdir}", checkIfExists: true)


	
    // Run MultiQC
    MULTIQC_REPORT(collect_results_ch)
}
