import csv
from cyvcf2 import VCF

# Input and output file paths
input_vcf = "annotated.vcf"
output_csv = "filtered_variants.csv"

# Open VCF and CSV
vcf = VCF(input_vcf)
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["CHROM", "POS", "REF", "ALT", "FILTER", "DP", "Impact", "Gene", "Variant"])

    for record in vcf:
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = ",".join(record.ALT)
        filter_status = record.FILTER or "PASS"
        dp = record.INFO.get("DP", 0)
        ann = record.INFO.get("ANN", "").split(",")[0].split("|")

        if len(ann) > 2:
            writer.writerow([chrom, pos, ref, alt, filter_status, dp, ann[2], ann[3], ann[1]])
