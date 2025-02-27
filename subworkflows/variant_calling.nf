nextflow.enable.dsl = 2

//  Include required modules
include { STAR_ALIGNMENT } from '../modules/star_alignment.nf'
include { SAMTOOLS_SORT_INDEX } from '../modules/samtools_sort_index.nf'
include { SAMTOOLS_FILTER_ORPHANS } from '../modules/samtools_filter_orphans.nf'
include { SAMTOOLS_FLAGSTAT } from '../modules/samtools_flagstat.nf'
include { GATK_MARK_DUPLICATES } from '../modules/gatk_mark_duplicates.nf'
include { SPLIT_NCIGAR_READS } from '../modules/split_ncigar_reads.nf'
include { SAMTOOLS_CALMD } from '../modules/samtools_calmd.nf'
include { GATK_RECALIBRATION } from '../modules/gatk_recalibration.nf'
include { BED_TO_INTERVAL_LIST } from '../modules/bed_to_interval_list.nf'
include { SCATTER_INTERVAL_LIST } from '../modules/scatter_interval_list.nf'
include { GATK_HAPLOTYPE_CALLER } from '../modules/gatk_haplotype_caller.nf'
include { GATK_VARIANT_FILTER } from '../modules/gatk_variant_filter.nf'
include { BCFTOOLS_STATS } from '../modules/bcftools_stats.nf'
include { BCFTOOLS_MERGE } from '../modules/bcftools_merge.nf'
include { BCFTOOLS_QUERY } from '../modules/bcftools_query.nf'

//  Include required annotation processes
include { ANNOTATE_INDIVIDUAL_VARIANTS } from '../modules/snpeff_annotate.nf'
include { ANNOTATE_INDIVIDUAL_VARIANTS_VEP } from '../modules/ensemblvep_annotate.nf'
include { ANNOTATE_VARIANTS } from '../modules/snpeff_annotate.nf'
include { ANNOTATEVARIANTS_VEP } from '../modules/ensemblvep_annotate.nf'
include { EXTRACT_VCF } from '../modules/extract_vcf.nf'
include { EXTRACT_individual_VCF } from '../modules/extract_individual_vcf.nf'

workflow VARIANT_CALLING {
    take:
        trimmed_reads_ch
        star_index
        reference_genome
        genome_index
        genome_dict
        merged_vcf
        merged_vcf_index
        denylist_bed
        snpeff_jar
        snpeff_config
        snpeff_db
        genomedb
        vep_cache
        clinvar_vcf
        clinvar_index
        

    main:
        log.info " Starting Variant Calling Workflow..."

        // **Step 1: STAR Alignment**
        star_aligned_ch = STAR_ALIGNMENT(trimmed_reads_ch, star_index)
		
		star_aligned_ch.view { "STAR Alignment Output: $it" }


        // **Step 2: Sort BAM files**
        sorted_bams = SAMTOOLS_SORT_INDEX(star_aligned_ch.map { tuple(it[0], it[1], it[5]) })
		
		sorted_bams.view { "Sorted BAMs Output: $it" }


        // **Step 3: Filter orphan reads**
        filtered_bams = SAMTOOLS_FILTER_ORPHANS(sorted_bams)
		
		filtered_bams.view { "Filtered BAMs Output: $it" }


        // **Step 4: Generate alignment statistics**
        alignment_stats = SAMTOOLS_FLAGSTAT(filtered_bams)
		
		alignment_stats.view { "Alignment stats Output: $it" }


        // **Step 5: Mark duplicates using GATK**
        marked_bams = GATK_MARK_DUPLICATES(filtered_bams)
		
        marked_bams.view { "Mark Duplicates Output: $it" }

        // **Step 6: Split N CIGAR Reads**
        split_bams = SPLIT_NCIGAR_READS(marked_bams.map { tuple(it[0], it[1], it[2], it[3]) }, reference_genome, genome_index, genome_dict)
		
		split_bams.view {"SplitNcigar Output : $it"  }


        // **Step 7: SAMTOOLS CALMD**
        calmd_bams = SAMTOOLS_CALMD(split_bams.map { tuple(it[0], it[1], it[2], it[3]) }, reference_genome, genome_index)

        // **Step 8: Base Quality Score Recalibration (BQSR)**
        recalibrated_bams = GATK_RECALIBRATION(calmd_bams.map { tuple(it[0], it[1], it[2], it[3]) }, reference_genome, genome_index, genome_dict, merged_vcf, merged_vcf_index)
		
		recalibrated_bams.view{"Recalibrated Bams Output: $it" }


        // **Step 9: Convert BED to Interval List**
        interval_list_ch = BED_TO_INTERVAL_LIST(denylist_bed, reference_genome, genome_dict)

        // **Step 10: Scatter the Interval List**
        scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, genome_dict)

        // **Step 11: Perform Variant Calling using GATK HaplotypeCaller**
        vcf_output = GATK_HAPLOTYPE_CALLER(recalibrated_bams.map { tuple(it[0], it[1], it[2], it[3]) }, reference_genome, genome_index, genome_dict, scattered_intervals_ch)

        vcf_output.view { "Raw VCF output: $it" }

        // **Step 12: Generate variant statistics**
        bcftools_stats_ch = BCFTOOLS_STATS(vcf_output)

        // **Step 13: Apply GATK Variant Filtering**
        filtered_individual_vcfs = GATK_VARIANT_FILTER(vcf_output, reference_genome, genome_index, genome_dict)
		
		filtered_individual_vcfs.view{" Filtered_vcf_output: $it"  }


        // **Step 14: Provide Stats**
        filtered_vcf_stats = BCFTOOLS_QUERY(filtered_individual_vcfs)

        if (params.merge_vcf) {
            log.info " Merging VCF files..."

            //  Collect filtered VCF paths into a channel
            filtered_vcf_list_ch = filtered_individual_vcfs.map { it[1] }.collect()

            //  Merge VCFs
            merged_filtered_vcfs = BCFTOOLS_MERGE(filtered_vcf_list_ch)

            println "Merging and annotating VCF files completed."

            // **Step 15: Annotate merged VCF with SnpEff**
            annotated_merged_vcf = ANNOTATE_VARIANTS(merged_filtered_vcfs, snpeff_jar, snpeff_config, snpeff_db, genomedb)

            // **Step 16: Annotate merged VCF with Ensembl VEP**
            annotated_merged_vcf_vep = ANNOTATEVARIANTS_VEP(merged_filtered_vcfs, vep_cache, clinvar_vcf, clinvar_index)
			
			// Step 17: Create a table from the annotated merged VCF
			extracted_csv = EXTRACT_VCF(annotated_merged_vcf)

            final_vcf_output = annotated_merged_vcf
        } else {
            log.info " Keeping individual VCFs..."

            // **Step 18: Annotate individual VCFs with SnpEff**
            annotated_individual_vcfs = ANNOTATE_INDIVIDUAL_VARIANTS(filtered_individual_vcfs, snpeff_jar, snpeff_config, snpeff_db, genomedb)

            // **Step 19: Annotate individual VCFs with Ensembl VEP**
            annotated_individual_vcf_vep = ANNOTATE_INDIVIDUAL_VARIANTS_VEP(filtered_individual_vcfs, vep_cache, clinvar_vcf, clinvar_index)

            // **Step 20: Convert Individual VCF to CSV**
            extracted_csv = EXTRACT_individual_VCF(annotated_individual_vcfs)




            final_vcf_output = annotated_individual_vcfs

        }

        log.info " Variant Annotation Completed."

    emit:
	
        final_vcf = final_vcf_output.annotated_vcf
		final_vcf_html = final_vcf_output.summary_html
        annotation_reports = extracted_csv
        star_logs = star_aligned_ch.map { [it[2], it[3], it[4]] }.flatten().collect()
        samtools_flagstat = alignment_stats.map { it[1] }.flatten().collect()
        gatk_metrics = marked_bams.map { it[4] }.flatten().collect()
        bcftools_stats = bcftools_stats_ch.map { it[2] }.flatten().collect()
        filtered_vcf_stats = filtered_vcf_stats.map { it[2] }.flatten().collect()
}
