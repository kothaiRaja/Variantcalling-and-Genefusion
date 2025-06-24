// ============================
// BUILD_REFERENCES Subworkflow
// ============================

include { DOWNLOAD_REFERENCE_GENOME }      from './prepare_references/download_genome.nf'
include { DOWNLOAD_GTF_ANNOTATION }        from './prepare_references/download_gtf.nf'
include { DOWNLOAD_AND_PREPARE_VARIANT_VCFS } from './prepare_references/prepare_vcf.nf'
include { BUILD_STAR_INDEX }               from './prepare_references/star_index.nf'
include { SNPEFF_SETUP }                   from './prepare_references/snpeff_setup.nf'
include { ARRIBA_SETUP }                   from './prepare_references/arriba.nf'
include { VEP_SETUP }                      from './prepare_references/VEP_setup.nf'

workflow BUILD_REFERENCES {

    take:
    
	samplesheet
	genome_id
	

    main:

	// Step 1: Read and validate the samplesheet
    samples_ch = Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample_id || !row.fastq_1 || !row.fastq_2) {
                error " Missing required fields in samplesheet. Each row must have: sample_id, fastq_1, fastq_2"
            }
            def strandedness = row.strandedness ?: "unstranded"
            tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)], strandedness)
        }





    
    // Download genome + indexes
    genome_refs = DOWNLOAD_REFERENCE_GENOME()
    reference_genome_ch       = genome_refs.genome
    reference_genome_index_ch = genome_refs.genome_index
    reference_genome_dict_ch  = genome_refs.genome_dict

    // Download GTF + exons BED
    gtf_outputs = DOWNLOAD_GTF_ANNOTATION()
    gtf_ch       = gtf_outputs.gtf
    exons_bed_ch = gtf_outputs.exons_bed

    // Prepare known variants
    prepare_variants_output = DOWNLOAD_AND_PREPARE_VARIANT_VCFS()
    known_variants_ch       = prepare_variants_output.merged_vcf
    known_variants_index_ch = prepare_variants_output.merged_vcf_index

    // Build STAR index
    star_out = BUILD_STAR_INDEX(reference_genome_ch, gtf_ch)
    star_index_ch = star_out.star_index

    // Set up snpEff
    snpeff = SNPEFF_SETUP(genome_id)
    snpeff_jar_ch    = snpeff.snpeff_jar
    snpeff_config_ch = snpeff.snpeff_config
    snpeff_db_dir_ch = snpeff.snpeff_db_dir

    // Setup arriba
    arriba = ARRIBA_SETUP()
    arriba_dir_ch = arriba.arriba_dir

    // Setup VEP
    vep = VEP_SETUP()
    vep_cache_ch   = vep.vep_cache
    vep_plugins_ch = vep.vep_plugins

    emit:
    reference_genome        = reference_genome_ch
    reference_genome_index  = reference_genome_index_ch
    reference_genome_dict   = reference_genome_dict_ch
    gtf_annotation          = gtf_ch
    exons_BED               = exons_bed_ch
    known_variants          = known_variants_ch
    known_variants_index    = known_variants_index_ch
    star_genome_index       = star_index_ch
    snpeff_jar              = snpeff_jar_ch
    snpeff_config           = snpeff_config_ch
    snpeff_db_dir           = snpeff_db_dir_ch
    arriba_dir              = arriba_dir_ch
    vep_cache               = vep_cache_ch
    vep_plugins             = vep_plugins_ch
	validated_reads			= samples_ch
}
