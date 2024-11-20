nextflow.enable.dsl = 2

// Include the QC Subworkflow
include { fastqc_raw } from './modules/QC/fastqc_raw.nf'
include { fastp } from './modules/QC/fastp.nf'
include { fastqc_trimmed } from './modules/QC/fastqc_trimmed.nf'
include { PREPARE_GENOME_PICARD } from './modules/prepare/dict_genome.nf'
include { PREPARE_GENOME_SAMTOOLS } from './modules/prepare/index_genome.nf'
include { PREPARE_VCF_FILE } from './modules/prepare/prepare_vcf.nf'
include { PREPARE_STAR_GENOME_INDEX } from './modules/prepare/starindex_genome.nf'
include { RNASEQ_MAPPING_STAR } from './modules/mapping/mapping_star.nf'
include { SAMTOOLS_FLAGSTAT  } from './modules/mapping/flagstat.nf'
include { SAMTOOLS_FILTER_INDEX  } from './modules/mapping/filter_index.nf'
include { ADD_READ_GROUP  } from './modules/mapping/addrg.nf'
include { RNASEQ_GATK_SPLITNCIGAR  } from './modules/split/splitncigar.nf'
include { SAMTOOLS_INDEX_SPLIT_BAM  } from './modules/split/index_splitbam.nf'
include { RNASEQ_GATK_RECALIBRATE  } from './modules/recalibrate/recalibrate.nf'
include { SAMTOOLS_INDEX_BAM  } from './modules/recalibrate/index_bam.nf'
include { RNASEQ_CALL_VARIANTS  } from './modules/variant_calling/variant.nf'
include { ANNOTATE_VARIANTS  } from './modules/annotations/snpeff.nf'
include { multiqc  } from './modules/multiqc/multiqc.nf'







workflow {
    // Step 1: Define the data directory dynamically
    def data_dir = params.data_dir ?: 'data/test' // Use test data by default if not specified

    // Step 2: Create a channel for paired-end reads
    reads_channel = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, reads -> tuple(sample_id, reads.sort()) }

    reads_channel.view()

    // Step 3: Run FastQC on raw reads
    fastqc_raw_results = fastqc_raw(reads_channel)

    // Step 4: Run Fastp for trimming
    fastp_results = fastp(reads_channel)

    // Step 5: Extract and regroup trimmed reads
    fastp_trimmed_reads = fastp_results.trimmed_reads
        .map { sample_id, r1, r2 -> tuple(sample_id, [r1, r2]) }

    // Step 6: Run FastQC on trimmed reads
    fastqc_trimmed_results = fastqc_trimmed(fastp_trimmed_reads)

    // Step 7: Data Preparation - Prepare genome indices and filtered VCF
    def genome_index_samtools = PREPARE_GENOME_SAMTOOLS(params.genome)
	def genome_dict = PREPARE_GENOME_PICARD(params.genome)
	def star_genome_dir = PREPARE_STAR_GENOME_INDEX(params.genome)
	def filtered_vcf = PREPARE_VCF_FILE(params.variants, params.denylist)


    // Step 8: STAR RNA-Seq Mapping
    star_aligned_bam = RNASEQ_MAPPING_STAR(star_genome_dir, fastp_trimmed_reads)

    // Step 9: Generate alignment statistics
    alignment_stats = SAMTOOLS_FLAGSTAT(star_aligned_bam)

    // Step 10: Process mapped reads with Samtools, Add Read Group, and SplitNCigar
    filtered_bam = SAMTOOLS_FILTER_INDEX(star_aligned_bam)
    bam_with_rg = ADD_READ_GROUP(filtered_bam)

    // Step 11: Base Recalibration and Indexing
    split_bam = RNASEQ_GATK_SPLITNCIGAR(params.genome, genome_index_samtools, genome_dict, bam_with_rg)
    indexed_split_bam = SAMTOOLS_INDEX_SPLIT_BAM(split_bam)

    // Step 12: Variant Calling 
    variant_vcf = RNASEQ_CALL_VARIANTS(params.genome, genome_index_samtools, genome_dict, indexed_split_bam)
 
	//Step 13: Annotations
	annotated_vcf = ANNOTATE_VARIANTS(variant_vcf)
    
	// Step 14: Run MultiQC to aggregate all results
    multiqc_results = multiqc(Channel.fromPath("${params.outdir}"))
}
