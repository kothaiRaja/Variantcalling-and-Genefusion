nextflow.enable.dsl = 2

// Process to download the test genome
process DOWNLOAD_TEST_GENOME {
    tag "Download test genome"

    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "genome.fa"

    script:
    """
    wget -q -O genome.fa ${params.test_data_url}/genome.fa
    """
}

// Process to download the test variants VCF
process DOWNLOAD_TEST_VARIANTS {
    tag "Download test variants VCF"

    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "variants.vcf"

    script:
    """
    wget -q -O variants.vcf ${params.test_data_url}/subset_fixed.vcf.gz
    """
}

// Process to download the test denylist BED
process DOWNLOAD_TEST_DENYLIST {
    tag "Download test denylist BED"

    publishDir "${params.test_data_dir}", mode: 'copy'

    output:
    path "denylist.bed"

    script:
    """
    wget -q -O denylist.bed ${params.test_data_url}/denylist.bed
    """
}

// Process to download and extract test reads
process DOWNLOAD_TEST_READS {
    tag "Download and extract test reads"

    publishDir "${params.test_data_dir}/reads", mode: 'copy'

    output:
    path "reads"

    script:
    """
    mkdir -p reads
    wget -q -O test_reads.tar.gz ${params.test_data_url}/reads/reads.tar.gz
    tar -xzf test_reads.tar.gz -C reads
    rm test_reads.tar.gz
    """
}

process CHECK_JAVA {
    tag "Check Java"

    output:
    path "java_check.txt"

    script:
    """
    if ! java -version &>/dev/null; then
        echo "Java is not installed or not in PATH." > java_check.txt
        exit 1
    else
        echo "Java is available." > java_check.txt
    fi
    """
}


process DOWNLOAD_SNPEFF_TOOL {
    tag "Download SnpEff Tool"
	publishDir "${params.test_data_dir}", mode: 'copy'
    output:
    path "${params.snpeff_jar_dir}/snpEff.jar"
	path "${params.snpeff_jar_dir}/snpEff.config"

    script:
    """
    mkdir -p ${params.snpeff_jar_dir}
    wget -q -O snpEff_latest_core.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip -j snpEff_latest_core.zip -d ${params.snpeff_jar_dir}
    rm snpEff_latest_core.zip
    """
}

process DOWNLOAD_SNPEFF_DB {
    tag "Download SnpEff Database"
	publishDir "${params.test_data_dir}/snpEff", mode: 'copy'
    input:
    val genome
    path snpeff_jar_path

    output:
    path "${params.snpeff_db_dir}/${genome}"

    script:
    """
    # Ensure the output directory exists first
    mkdir -p ${params.snpeff_db_dir}

    # Use an absolute path for the data directory
    data_dir=\$(realpath ${params.snpeff_db_dir})

    # Download the database
    java -Xmx4g -Xms2g -jar ${snpeff_jar_path} download ${genome} -dataDir \$data_dir -v
    """
}





workflow {
    // Step 1: Download test references
    def genome = DOWNLOAD_TEST_GENOME()
    def variants = DOWNLOAD_TEST_VARIANTS()
    def denylist = DOWNLOAD_TEST_DENYLIST()
    def reads = DOWNLOAD_TEST_READS()

    // Check for Java installation
    CHECK_JAVA()

    // Download the SnpEff tool
    def snpeff_jar_and_config = DOWNLOAD_SNPEFF_TOOL()

    // Download the SnpEff database
    DOWNLOAD_SNPEFF_DB(params.genomedb, snpeff_jar_and_config[0])
    
}
