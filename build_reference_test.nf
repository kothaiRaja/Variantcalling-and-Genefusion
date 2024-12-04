nextflow.enable.dsl = 2

include { DOWNLOAD_TEST_GENOME } from './modules/ref_download/ref_genome_test.nf'
include { DOWNLOAD_TEST_VARIANTS } from './modules/ref_download/known_variants_test.nf'
include { DOWNLOAD_TEST_DENYLIST } from './modules/ref_download/denylist_test.nf'
include { DOWNLOAD_TEST_GTF } from './modules/ref_download/gft_file_test.nf'
//include { DOWNLOAD_TEST_FUSION } from './modules/ref_download/known_fusion_test.nf'
include { create_fasta_index } from './modules/prepare_refs/fasta_index.nf'
include { create_genome_dict } from './modules/prepare_refs/genome_dict.nf'
include { create_star_index } from './modules/prepare_refs/star_index.nf'
include { PREPARE_VCF_FILE } from './modules/prepare_refs/filtered_vcf.nf'
include { CHECK_JAVA } from './modules/Tools/check_java.nf'
include { DOWNLOAD_SNPEFF_TOOL } from './modules/Tools/SnpEff_tool.nf'
include { DOWNLOAD_SNPEFF_DB } from './modules/Tools/SnpEff_DB.nf'
include { fastqc_raw } from './modules/Quality_control/fastqc.nf'
include { trim_reads } from './modules/Quality_control/fastp.nf'
//include { multiqc  } from './modules/Quality_control/multiqc.nf'


workflow {
    // Step 1: Download reference files
    def genome = DOWNLOAD_TEST_GENOME()
    def variants = DOWNLOAD_TEST_VARIANTS()
    def denylist = DOWNLOAD_TEST_DENYLIST()
    def genome_gtf = DOWNLOAD_TEST_GTF()
	
	// Step 2: Create genome index files
    fasta_index = create_fasta_index(genome)
    
	// Step 3: Genome dictionary depends on genome download
    genome_dict = create_genome_dict(genome)
    
	// Step 4: STAR index depends on genome and GTF downloads
    star_index = create_star_index(genome, genome_gtf)
	
	// Step 5: Prepare VCF file
    prepared_vcf_ch = PREPARE_VCF_FILE(variants, denylist)

    // Step 6: Check Java installation
    CHECK_JAVA()

    // Step 7: Download and configure SnpEff
    def snpeff_jar_and_config = DOWNLOAD_SNPEFF_TOOL()
    def snpeff_db = DOWNLOAD_SNPEFF_DB(params.genomedb, snpeff_jar_and_config[0])

    // Step 8: Load and parse sample metadata
    samples_channel = Channel.fromPath(params.csv_file)
        .splitCsv(header: true)
        .map { row ->  tuple(row.sample_id, row.fastq_1, row.fastq_2) }
		
    // Step 9: Download raw reads
    //reads_channel = download_reads(samples_channel)

    // Step 10: Perform FastQC on raw reads (depends on downloaded reads)
    fastqc_raw_results = fastqc_raw(samples_channel)

    // Step 11: Trim reads with fastp (depends on FastQC)
    trimmed_reads = trim_reads(samples_channel)

    
}

