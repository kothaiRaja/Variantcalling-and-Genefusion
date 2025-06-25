include { MultiQC } from '../modules/multiqc_quality/main.nf'

workflow MULTIQC_WRAPPER {

    take:
        fastqc_results
        fastp_reports
        star_logs
        samtools_flagstat
        gatk_metrics
        bcftools_stats
        filtered_vcf_stats
        annotation_reports
        fusion_visuals
        maf_reports
        version_yamls

    main:
        all_reports = Channel.empty()
            .mix(fastqc_results.ifEmpty([]))
            .mix(fastp_reports.ifEmpty([]))
            .mix(star_logs.ifEmpty([]))
            .mix(samtools_flagstat.ifEmpty([]))
            .mix(gatk_metrics.ifEmpty([]))
            .mix(bcftools_stats.ifEmpty([]))
            .mix(filtered_vcf_stats.ifEmpty([]))
            .mix(annotation_reports.ifEmpty([]))
            .mix(fusion_visuals.ifEmpty([]))
            .mix(maf_reports.ifEmpty([]))
            .mix(version_yamls.ifEmpty([]))
            .collect()
			
			

       

        
      reports_mqc = MultiQC(all_reports)
	  reports_ch = MultiQC.out.report
	  version_ch = MultiQC.out.versions
	  

    emit:
        reports = reports_ch
		versions = version_ch
}
