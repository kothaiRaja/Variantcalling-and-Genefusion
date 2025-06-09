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

    awk '\$3 == "exon" {
        split(\$9, a, ";");
        gene_id = a[1]; gsub(/"|gene_id /, "", gene_id);
        print \$1 "\\t" \$4-1 "\\t" \$5 "\\t" gene_id
    }' ${annotation_gtf} | sort -k1,1 -k2,2n > exons.bed

    echo "Exons BED file generated: exons.bed"
    """
}
