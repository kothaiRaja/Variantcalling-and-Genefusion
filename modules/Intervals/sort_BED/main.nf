process SORT_BED_BY_FAIDX {
  tag "Sort BED by faidx"
  label 'process_low'
  container "https://depot.galaxyproject.org/singularity/bedtools%3A2.30.0--h7d7f7ad_1"
  publishDir params.bed_to_interval_outdir, mode: "copy"

  input:
  path bed_in
  path faidx    

  output:
  path "exons.sorted.bed", emit: bed_sorted

  script:
  """
  set -euo pipefail
  bedtools sort -faidx "${faidx}" -i "${bed_in}" > exons.sorted.bed
  """
}
