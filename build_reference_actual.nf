nextflow.enable.dsl = 2

include { DOWNLOAD_REF_GENOME } from './modules/REFERENCES_DOWNLOAD/ref_genome_actual.nf'
include { DOWNLOAD_VARIANTS_SNP } from './modules/REFERENCES_DOWNLOAD/known_variants_snps_actual.nf'
include { DOWNLOAD_VARIANTS_INDELS } from './modules/REFERENCES_DOWNLOAD/known_variants_indels_actual.nf'
include { DOWNLOAD_DENYLIST } from './modules/REFERENCES_DOWNLOAD/denylist_actual.nf'
include { DOWNLOAD_GTF } from './modules/REFERENCES_DOWNLOAD/gtf_actual.nf'
include { CREATE_FASTA_INDEX } from './modules/PREPARE_REFs/fasta_index_actual.nf'
include { CREATE_GENOME_DICT } from './modules/PREPARE_REFs/genome_dict_actual.nf'
include { CREATE_STAR_INDEX } from './modules/PREPARE_REFs/star_index_actual.nf'
include { PREPARE_VCF_FILE } from './modules/PREPARE_REFs/filtered_vcf_actual.nf'
include { CHECK_JAVA } from './modules/TOOLS and DBs/check_java_actual.nf'
include { DOWNLOAD_SNPEFF_TOOL } from './modules/TOOLS and DBs/SnpEff_tool_actual.nf'
include { DOWNLOAD_SNPEFF_DB } from './modules/TOOLS and DBs/SnpEff_DB_actual.nf'
include { DOWNLOAD_ARRIBA } from './modules/TOOLS and DBs/Arriba_tool_actual.nf'
include { DOWNLOAD_VEP_CACHE  } from './modules/VEP_ANNOTATIONS/vep_prepare_actual.nf'
include { DOWNLOAD_CLINVAR  } from './modules/VEP_ANNOTATIONS/clinvar_actual.nf'
include { FASTQC_RAW } from './modules/QUALITY_CONTROL/fastqc_actual.nf'
include { TRIM_READS } from './modules/QUALITY_CONTROL/fastp_actual.nf'



workflow {
    //=======Step 1: Download reference files==========//
    def genome = DOWNLOAD_REF_GENOME()
    def variants_snp = DOWNLOAD_VARIANTS_SNP()
	def variants_indels = DOWNLOAD_VARIANTS_INDELS()
    def denylist = DOWNLOAD_DENYLIST()
    def genome_gtf = DOWNLOAD_GTF()
	
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
	
	//========Step 8: Download ensembl_vep dependencies==========//
	vep_cache_dir = DOWNLOAD_VEP_CACHE()
    clinvar_files = DOWNLOAD_CLINVAR()

    //=========Step 8: Load and parse sample metadata==============//
    samples_channel = Channel.fromPath(params.actual_csv_file)
        .splitCsv(header: true)
        .map { row ->  tuple(row.sample_id, row.fastq_1, row.fastq_2) }
		
    //=========Step 9: Download raw reads===============//
    //reads_channel = download_reads(samples_channel)

    //=========Step 10: Perform FastQC on raw reads=========//
    fastqc_raw_results = FASTQC_RAW(samples_channel)

    // Step 11: Trim reads with fastp (depends on FastQC)
    trimmed_reads = TRIM_READS(samples_channel)

    
}

