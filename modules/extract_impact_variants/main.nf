process EXTRACT_VARIANT_IMPACT {
    tag { sample_id }
    publishDir "${params.outdir}/annotations", mode: 'copy'
    container params.container

    cpus 2
    memory '2GB'
    time '30m'

    input:
    tuple val(sample_id), path(annotated_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_variants.csv")

    script:
    """
    echo "Running Variant Extraction..."
   
    python3 <<EOF
import pandas as pd
import vcfpy

# File paths
vcf_file = "${annotated_vcf}"
output_csv = "${sample_id}_variants.csv"

# Initialize list to store parsed VCF data
variant_data = []

# Parse the VCF file
reader = vcfpy.Reader.from_path(vcf_file)

# Extract sample names
sample_names = reader.header.samples.names

# Extract data from VCF
for record in reader:
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt_list = [str(a) for a in record.ALT] if record.ALT else ["NA"]
    filt = ",".join(record.FILTER) if record.FILTER else "PASS"
    dp = record.INFO.get('DP', 'NA')  

    # Determine variant type (SNV, Insertion, Deletion)
    for alt in alt_list:
        if len(ref) == len(alt) == 1:
            variant_type = "SNV"
        elif len(alt) > len(ref):
            variant_type = "Insertion"
        elif len(ref) > len(alt):
            variant_type = "Deletion"
        else:
            variant_type = "Complex"

        # Extract ANN field (annotation)
        ann_field = record.INFO.get('ANN', [])
        if ann_field:
            for annotation in ann_field:
                ann_details = annotation.split('|')
                gene = ann_details[3] if len(ann_details) > 3 else "NA"
                impact = ann_details[2] if len(ann_details) > 2 else "NA"
                effect = ann_details[1] if len(ann_details) > 1 else "NA"

                # Iterate over samples to include sample-specific data
                for sample in sample_names:
                    sample_data = record.call_for_sample.get(sample, None)
                    genotype = sample_data.data.get('GT', './.') if sample_data else './.'

                    variant_data.append({
                        "Sample_ID": sample,
                        "Chromosome": chrom,
                        "Position": pos,
                        "Reference": ref,
                        "Alternate": alt,
                        "Variant_Type": variant_type,
                        "Filter": filt,
                        "DP": dp,
                        "Impact": impact,
                        "Effect": effect,
                        "Gene": gene,
                        "Genotype": genotype
                    })

# Convert the VCF data to a DataFrame
vcf_df = pd.DataFrame(variant_data)

# Remove exact duplicate rows
vcf_df = vcf_df.drop_duplicates()

# Save to CSV
vcf_df.to_csv(output_csv, index=False)

print(f"Extracted variant data saved to: {output_csv}")
EOF
    """
}
