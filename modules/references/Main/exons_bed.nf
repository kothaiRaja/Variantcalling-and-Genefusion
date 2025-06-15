process GENERATEEXONS_BED {
    tag "Generate Exons BED File"
    label 'process_medium'
    
    container "https://depot.galaxyproject.org/singularity/bedtools%3A2.30.0--h7d7f7ad_1"
    publishDir "${params.main_data_dir}/reference", mode: 'copy'

    input:
    path(annotation_gtf)

    output:
    path("exons.bed"), emit: exons_bed

    script:
    """
    echo "Extracting exon regions from GTF file to BED format..."

    awk '\$3 == "exon" {print \$1 "\\t" \$4-1 "\\t" \$5 "\\t" \$9}' $annotation_gtf | sort -k1,1 -k2,2n > exons.bed
    echo "Exons BED file generated: exons.bed"

    echo "Checking Chromosome Naming in BED file..."
    FIRST_CHR=\$(awk '\$1 !~ /^#/ {print \$1; exit}' exons.bed)
    echo "  → Detected contig name in BED: '\$FIRST_CHR'"
    """
}
