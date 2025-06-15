// ========================== Download GTF Annotation File ========================== //
process CHECK_OR_DOWNLOAD_GTF {
    tag "Check or Download GTF File"
    container null
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    output: 
    path "annotations.gtf", emit: gtf

    script:
    """
    echo "Downloading GTF annotation from: ${params.gtf_download_url}"
    wget -q -O annotations.gtf.gz ${params.gtf_download_url}

    echo "Unzipping GTF file..."
    if file annotations.gtf.gz | grep -q 'gzip'; then
        gunzip annotations.gtf.gz
    else
        mv annotations.gtf.gz annotations.gtf
    fi

    echo "Checking Chromosome Naming Format in annotations.gtf..."
    FIRST_CHR=\$(awk '\$1 !~ /^#/ {print \$1; exit}' annotations.gtf)
    echo "  → Detected contig name in GTF: '\$FIRST_CHR'"
    
    """
}

