nextflow.enable.dsl = 2

include { DOWNLOAD_REF_GENOME } from './modules/REFERENCES_DOWNLOAD/ref_genome_test.nf'
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
include { DOWNLOAD_VEP_CACHE  } from './modules/VEP_ANNOTATIONS/vep_prepare.nf'
include { DOWNLOAD_CLINVAR  } from './modules/VEP_ANNOTATIONS/clinvar.nf'
include { FASTQC_RAW } from './modules/QUALITY_CONTROL/fastqc.nf'
include { TRIM_READS } from './modules/QUALITY_CONTROL/fastp.nf'



workflow {

	// Define paths for all required files
    def genome_path = "${params.test_data_dir}/reference/genome.fa"
    def variants_snp_path = "${params.test_data_dir}/reference/variants_snp.vcf"
    def variants_indels_path = "${params.test_data_dir}/reference/variants_indels.vcf"
    def denylist_path = "${params.test_data_dir}/reference/denylist.bed"
    def genome_gtf_path = "${params.test_data_dir}/reference/annotations.gtf"
	def fasta_index_path = "${params.test_data_dir}/genome.fa.fai"
    def genome_dict_path = "${params.test_data_dir}/genome.dict"
    def star_index_path = "${params.test_data_dir}/STAR_index"
    def prepared_vcf_path = "${params.test_data_dir}/merged.filtered.recode.vcf.gz"
	
     // Load and parse sample metadata
    def samples_channel = Channel.fromPath(params.test_csv_file)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.fastq_1, row.fastq_2) }

    if (params.only_fastqc_fastp) {
        log.info "Running only FastQC and Fastp as 'params.only_fastqc_fastp' is set to true"

        // Execute only FastQC and Fastp
        def fastqc_raw_results = fastqc_raw(samples_channel)
        def trimmed_reads = trim_reads(samples_channel)
    } else {
        log.info "Running full pipeline setup and all processes"

        // Step 1: Download reference files (if they do not exist)
        def genome = file(genome_path).exists() ? file(genome_path) : DOWNLOAD_TEST_GENOME()
        def variants_snp = file(variants_snp_path).exists() ? file(variants_snp_path) : DOWNLOAD_TEST_VARIANTS_SNP()
        def variants_indels = file(variants_indels_path).exists() ? file(variants_indels_path) : DOWNLOAD_TEST_VARIANTS_INDELS()
        def denylist = file(denylist_path).exists() ? file(denylist_path) : DOWNLOAD_TEST_DENYLIST()
        def genome_gtf = file(genome_gtf_path).exists() ? file(genome_gtf_path) : DOWNLOAD_TEST_GTF()

        // Step 2: Prepare reference files (if they do not exist)
        def fasta_index = file(fasta_index_path).exists() ? file(fasta_index_path) : create_fasta_index(genome)
        def genome_dict = file(genome_dict_path).exists() ? file(genome_dict_path) : create_genome_dict(genome)
        def star_index = file(star_index_path).exists() ? file(star_index_path) : create_star_index(genome, genome_gtf)

        def prepared_vcf = file(prepared_vcf_path).exists() ? file(prepared_vcf_path) : PREPARE_VCF_FILE(variants_snp, variants_indels, denylist)

        
		// Step 3: Setup tools (if they do not exist)
		log.info "Checking Java installation..."
		def java_check = CHECK_JAVA()

		log.info "Checking SnpEff tool existence..."
		def snpeff_jar_path = "${params.snpeff_jar_dir}/snpEff.jar"
		def snpeff_config_path = "${params.snpeff_jar_dir}/snpEff.config"
		// Execute or skip SnpEff tool download
		def snpeff_tool = (!file(snpeff_jar_path).exists() || !file(snpeff_config_path).exists()) ?
			DOWNLOAD_SNPEFF_TOOL() :
			[file(snpeff_jar_path), file(snpeff_config_path)]

		log.info "SnpEff JAR exists: ${file(snpeff_jar_path).exists()}"
		log.info "SnpEff Config exists: ${file(snpeff_config_path).exists()}"

		// Step 2: Check for SnpEff database existence
		log.info "Checking SnpEff database existence..."
		def snpeff_db_path = "${params.snpeff_db_dir}/${params.genomedb}"
		def snpeff_db = !file(snpeff_db_path).exists() ?
			DOWNLOAD_SNPEFF_DB(params.genomedb, snpeff_tool[0]) :
			file(snpeff_db_path)

		log.info "SnpEff DB exists: ${file(snpeff_db_path).exists()}"


		log.info "Checking Arriba tool existence..."
		def arriba_dir_path = "${params.test_data_dir}/ARRIBA/arriba_v2.4.0"
		def arriba_tool = file(arriba_dir_path).exists() ? 
			file(arriba_dir_path) : 
			DOWNLOAD_ARRIBA()

		log.info "Checking VEP cache existence..."
		def vep_cache_path = "${params.test_data_dir}/vep_cache/homo_sapiens/110_GRCh38"
		def vep_cache = file(vep_cache_path).exists() ? 
			file(vep_cache_path) : 
			DOWNLOAD_VEP_CACHE()

		log.info "Checking ClinVar files existence..."
		def clinvar_vcf_path = "${params.test_data_dir}/VEP/clinvar.vcf.gz"
		def clinvar_tbi_path = "${params.test_data_dir}/VEP/clinvar.vcf.gz.tbi"
		def clinvar_files = (file(clinvar_vcf_path).exists() && file(clinvar_tbi_path).exists()) ?
			[file(clinvar_vcf_path), file(clinvar_tbi_path)] : 
			DOWNLOAD_CLINVAR()

		// Log final existence checks
		log.info "Arriba Tool: ${file(arriba_dir_path).exists()}"
		log.info "VEP Cache: ${file(vep_cache_path).exists()}"
		log.info "ClinVar VCF: ${file(clinvar_vcf_path).exists()}, TBI: ${file(clinvar_tbi_path).exists()}"


        // Step 4: Preprocessing (FastQC and Fastp)
        def fastqc_raw_results = FASTQC_RAW(samples_channel)
        def trimmed_reads = TRIM_READS(samples_channel)
    }
}


