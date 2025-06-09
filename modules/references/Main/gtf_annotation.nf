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
	
	# Decompress or rename if not gzipped
if file annotations.gtf.gz | grep -q 'gzip'; then
    gunzip annotations.gtf.gz
else
    mv annotations.gtf.gz annotations.gtf
fi

# Safety: Ensure output exists before proceeding
if [ ! -f annotations.gtf ]; then
    echo "Error: GTF file missing after download/decompression." >&2
    exit 1
fi

# Detect and convert chr-prefixed chromosomes
FIRST_CHR=\$(awk '\$1 !~ /^#/ {print \$1; exit}' annotations.gtf)

if [[ "\$FIRST_CHR" == chr* ]]; then
    echo "GTF file uses 'chr' prefix. Converting to numeric format..."
    awk '{ if (\$1 ~ /^chr/) sub(/^chr/, "", \$1); print }' annotations.gtf > annotations_numeric.gtf

    if [[ -s annotations_numeric.gtf ]]; then
        mv annotations_numeric.gtf annotations.gtf
        echo "GTF chromosome names converted to numeric format."
    else
        echo "Error: GTF conversion failed — keeping original." >&2
    fi
else
    echo "GTF chromosome names already in numeric format."
fi

    
    """
}

