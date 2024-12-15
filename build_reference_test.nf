nextflow.enable.dsl = 2

include { DOWNLOAD_TEST_GENOME } from './modules/REFERENCES_DOWNLOAD/ref_genome_test.nf'
include { DOWNLOAD_TEST_VARIANTS_SNP } from './modules/REFERENCES_DOWNLOAD/known_variants_snps.nf'
include { DOWNLOAD_TEST_VARIANTS_INDELS } from './modules/REFERENCES_DOWNLOAD/known_variants_indels.nf'
include { DOWNLOAD_TEST_DENYLIST } from './modules/REFERENCES_DOWNLOAD/denylist_test.nf'
include { DOWNLOAD_TEST_GTF } from './modules/REFERENCES_DOWNLOAD/gft_file_test.nf'
include { CREATE_FASTA_INDEX } from './modules/PREPARE_REFs/fasta_index.nf'
include { CREATE_GENOME_DICT } from './modules/PREPARE_REFs/genome_dict.nf'
include { CREATE_STAR_INDEX } from './modules/PREPARE_REFs/star_index.nf'
include { PREPARE_VCF_FILE } from './modules/PREPARE_REFs/filtered_vcf.nf'
include { CHECK_JAVA } from './modules/TOOLS and DBs/check_java.nf'
include { DOWNLOAD_SNPEFF_TOOL } from './modules/TOOLS and DBs/SnpEff_tool.nf'
include { DOWNLOAD_SNPEFF_DB } from './modules/TOOLS and DBs/SnpEff_DB.nf'
include { DOWNLOAD_ARRIBA } from './modules/TOOLS and DBs/Arriba_tool.nf'
include { FASTQC_RAW } from './modules/QUALITY_CONTROL/fastqc.nf'
include { TRIM_READS } from './modules/QUALITY_CONTROL/fastp.nf'



workflow {
    //=======Step 1: Download reference files==========//
    def genome = DOWNLOAD_TEST_GENOME()
    def variants_snp = DOWNLOAD_TEST_VARIANTS_SNP()
	def variants_indels = DOWNLOAD_TEST_VARIANTS_INDELS()
    def denylist = DOWNLOAD_TEST_DENYLIST()
    def genome_gtf = DOWNLOAD_TEST_GTF()
	
	//=======Step 2: Create genome index files===========//
    fasta_index = CREATE_FASTA_INDEX(genome)
    
	//=======Step 3: Genome dictionary depends on genome download=========//
    genome_dict = CREATE_GENOME_DICT(genome)
    
	//=======Step 4: STAR index depends on genome and GTF downloads=======//
    star_index = CREATE_STAR_INDEX(genome, genome_gtf)
	
	//=======Step 5: Prepare VCF file================//
    prepared_vcf_ch = PREPARE_VCF_FILE(variants_snp, variants_indels, denylist)

    //========Step 6: Check Java installation=========//
    CHECK_JAVA()

    //=========Step 7: Download and configure SnpEff============//
    def snpeff_jar_and_config = DOWNLOAD_SNPEFF_TOOL()
    def snpeff_db = DOWNLOAD_SNPEFF_DB(params.genomedb, snpeff_jar_and_config[0])
	Arriba_Toolsetup = DOWNLOAD_ARRIBA()

    //=========Step 8: Load and parse sample metadata==============//
    samples_channel = Channel.fromPath(params.csv_file)
        .splitCsv(header: true)
        .map { row ->  tuple(row.sample_id, row.fastq_1, row.fastq_2) }
		
    //=========Step 9: Download raw reads===============//
    //reads_channel = download_reads(samples_channel)

    //=========Step 10: Perform FastQC on raw reads=========//
    fastqc_raw_results = FASTQC_RAW(samples_channel)

    // Step 11: Trim reads with fastp (depends on FastQC)
    trimmed_reads = TRIM_READS(samples_channel)

    
}

