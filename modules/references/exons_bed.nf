process GENERATEEXONS_BED {
    tag "Generate Exons BED File"
    label 'process_medium'
    
    container "https://depot.galaxyproject.org/singularity/bedtools%3A2.30.0--h7d7f7ad_1"
    publishDir "${params.ref_base}/reference", mode: 'copy'

    input:
    path(annotation_gtf)

    output:
    path("exons.bed"), emit: exons_bed

    when:
    !params.exons_bed && !file("${params.ref_base}/reference/exons.bed").exists()

	script:
"""
echo "Extracting exon regions from GTF file to BED format..."

awk '\$3 == "exon" {
    match(\$9, /gene_id "([^"]+)"/, gene_id)
    match(\$9, /transcript_id "([^"]+)"/, transcript_id)
    match(\$9, /exon_number "([^"]+)"/, exon_number)
    name = gene_id[1] "|" transcript_id[1] "|" exon_number[1]
    print \$1 "\\t" \$4-1 "\\t" \$5 "\\t" name
}' ${annotation_gtf} | sort -k1,1 -k2,2n > exons.bed

echo "Exons BED file generated: exons.bed"

FIRST_CHR=\$(awk '\$1 !~ /^#/ {print \$1; exit}' exons.bed)
echo "  → Detected contig name in BED: '\$FIRST_CHR'"
"""

    
}
