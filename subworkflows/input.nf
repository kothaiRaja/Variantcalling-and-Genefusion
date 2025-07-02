workflow INPUT_PAIRED_READS {
    take:
    samplesheet

    main:
    Channel
        .from(samplesheet)
        .splitCsv(header: true)
        .map { row -> create_paired_meta(row) }
        .view { "Parsed sample: $it" }  
        .set { paired_reads }


    emit:
    paired_reads
}

// Function to create [meta, [fastq_1, fastq_2]]
def create_paired_meta(LinkedHashMap row) {
    def meta = [
        id: row.sample_id,
        strandedness: row.strandedness
    ]

    // Check existence of fastq_1 and fastq_2
    def fq1 = file(row.fastq_1)
    def fq2 = file(row.fastq_2)

    if (!fq1.exists()) {
        exit 1, "ERROR: fastq_1 file not found for ${row.sample_id}: ${row.fastq_1}"
    }
    if (!fq2.exists()) {
        exit 1, "ERROR: fastq_2 file not found for ${row.sample_id}: ${row.fastq_2}"
    }

    return [ meta, [ fq1, fq2 ] ]
}
