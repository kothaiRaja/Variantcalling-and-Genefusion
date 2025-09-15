process GENERATEEXONS_BED {
  tag "Generate Exons BED File"
  label 'process_medium'

  container "https://depot.galaxyproject.org/singularity/bedtools%3A2.30.0--h7d7f7ad_1"
  publishDir "${params.ref_base}/reference", mode: 'copy'

  input:
  path annotation_gtf        
               

  output:
  path "exons.bed", emit: exons_bed

  when:
  !params.exons_bed && !file("${params.ref_base}/reference/exons.bed").exists()

  script:
  """
  set -euo pipefail
  awk 'BEGIN{OFS="\\t"} \$0!~/^#/ && \$3=="exon"{ print \$1,\$4-1,\$5 }' "${annotation_gtf}" > exons.tmp.bed
  bedtools sort -i exons.tmp.bed > exons.bed
  rm -f exons.tmp.bed
  """
}
