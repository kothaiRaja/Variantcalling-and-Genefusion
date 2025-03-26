// ========================== Download GTF Annotation File ========================== //
process CHECK_OR_DOWNLOAD_GTF {
    tag "Check or Download GTF File"
    container null
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    output: 
    path "annotations.gtf", emit: gtf

    script:
    """
    wget -q -O annotations.gtf.gz ${params.gtf_download_url}

    # Check if the file is gzipped and unzip if necessary
    if file annotations.gtf.gz | grep -q 'gzip'; then
        gunzip annotations.gtf.gz
    else
        mv annotations.gtf.gz annotations.gtf
    fi

    echo "Checking Chromosome Naming in GTF File..."

    # Extract the first non-comment chromosome name
    FIRST_CHR=\$(awk '\$1 !~ /^#/ {print \$1; exit}' annotations.gtf)

    if [[ "\$FIRST_CHR" == chr* ]]; then
        echo " GTF file uses 'chr' prefix. Converting to numeric format..."

        # Convert "chr1" -> "1", "chrX" -> "X" (entire file)
        awk '{ if (\$1 ~ /^chr/) sub(/^chr/, "", \$1); print }' annotations.gtf > annotations_numeric.gtf

        # Ensure the conversion was successful before replacing the file
        if [[ -s annotations_numeric.gtf ]]; then
            mv annotations_numeric.gtf annotations.gtf
            echo " GTF Chromosomes converted to numeric format."
        else
            echo " Error: GTF conversion failed, keeping original file."
        fi
    else
        echo " GTF file is already in numeric format."
    fi
    """
}

