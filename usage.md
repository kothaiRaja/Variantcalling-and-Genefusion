# RNA-Seq Variant Calling and RNA Fusion Detection Pipeline

## Overview

This pipeline is designed for RNA-seq variant calling and RNA fusion detection. It is modularized to ensure scalability, reproducibility, and ease of maintenance. Key functionalities include:

- Quality control
- Alignment
- Variant identification and annotation
- RNA fusion detection

Integrated tools include **FastQC**, **STAR**, **GATK**, **SnpEff**, and **Arriba**, ensuring high accuracy and robustness.

## Objectives

- **Genetic Variant Identification**: Detect genetic variants from RNA-seq data.
- **RNA Fusion Detection**: Identify RNA fusion events in transcriptomic data.
- **Comprehensive Reporting**: Provide automated results with clear reporting.

## Usage

This pipeline leverages **Nextflow DSL 2** for modular design, enabling flexibility and maintainability. Below are the key modules and their functionalities:

### References Download

- `DOWNLOAD_REF_GENOME`: Downloads the reference genome.
- `DOWNLOAD_VARIANTS_SNP` & `DOWNLOAD_VARIANTS_INDELS`: Downloads SNP and indel variant files.
- `DOWNLOAD_DENYLIST`: Downloads blacklist regions.
- `DOWNLOAD_GTF`: Fetches genome annotation files.

### Reference Preparation

- `CREATE_FASTA_INDEX`: Generates an index file for the reference genome.
- `CREATE_GENOME_DICT`: Creates a dictionary for the genome.
- `CREATE_STAR_INDEX`: Builds STAR indices for alignment.
- `PREPARE_VCF_FILE`: Merges and filters VCF files for known variants.

### Tools and Databases

- `CHECK_JAVA`: Verifies the Java environment.
- `DOWNLOAD_SNPEFF_TOOL`: Downloads the SnpEff tool.
- `DOWNLOAD_SNPEFF_DB`: Fetches SnpEff database.
- `DOWNLOAD_ARRIBA`: Downloads the Arriba tool for RNA fusion detection.
- `DOWNLOAD_VEP_CACHE`: Prepares the VEP cache for annotation.
- `DOWNLOAD_CLINVAR`: Fetches ClinVar files for variant annotation.

### Quality Control

- `FASTQC_RAW`: Performs quality checks on raw FASTQ files.
- `TRIM_READS`: Trims low-quality reads using Fastp.

## Workflow Steps

### 1. Pipeline Initialization

- Parses sample metadata and loads file paths for required references.
- Conditional execution based on user-defined parameters (e.g., `params.only_fastqc_fastp`).

### 2. Reference and Tool Setup

- **Initial Checks**: The processes check for the existence of required reference files, tools, databases, and prepared files. If a file is unavailable, the process downloads it from the respective URL; otherwise, the process is skipped.

- **Renaming Reference Files**: If already available, reference files are renamed as follows:

  1. Whole genome GRCh38: `genome.fa`
  2. GTF files: `annotations.gtf`
  3. Blacklist: `denylist.bed`
  4. Variants indels: `variants_indels.vcf`
  5. Variants SNP: `variants_snp.vcf`

- Downloads and prepares reference genome, annotation, and variant files.

- Builds necessary indices and checks tool availability (e.g., SnpEff, Arriba).

### 3. Data Preprocessing

- Executes quality control using FastQC and Fastp.
- Optionally runs only FastQC and Fastp if `params.only_fastqc_fastp` is set.

### 4. Integration with Tools

- Annotates variants using SnpEff.
- Detects RNA fusion events using Arriba.
- Utilizes VEP for additional annotations (e.g., ClinVar data).

### 5. Conditional Logic

- Steps are executed only if required files are unavailable, ensuring efficiency and avoiding redundancy.

### 6. Output

- Results include:
  - Trimmed reads
  - Quality reports
  - Alignment files
  - Variant calls
  - Fusion detections
  - Final summary reports
  - **Final Output**: All intermediate files and final annotated files are located in the `output` folder.

## Sample Sheet Format

The pipeline is designed for **paired-end RNA-seq data**. Below is an example of how the sample sheet should look:

| sample\_id | fastq\_1                        | fastq\_2                        |
| ---------- | ------------------------------- | ------------------------------- |
| Sample\_01 | /path/to/sample\_1\_R1.fastq.gz | /path/to/sample\_1\_R2.fastq.gz |
| Sample\_02 | /path/to/sample\_2\_R1.fastq.gz | /path/to/sample\_2\_R2.fastq.gz |

- **sample\_id**: Unique identifier for each sample.
- **fastq\_1**: Path to the forward reads file.
- **fastq\_2**: Path to the reverse reads file.

## Command to Run

To test the pipeline with default test datasets:

```bash
nextflow run build_reference_test.nf -c nextflow_ref.config --mode test --only_fastqc_fastp false --mode test -profile singularity
```

To run only FastQC and Fastp with all available reference files:

```bash
nextflow run build_reference_test.nf -c nextflow_ref.config --mode test --only_fastqc_fastp true --mode test -profile singularity
```

Once the reference files and input files are ready, run:

```bash
nextflow run main.nf -c nextflow_main.config --mode test --merge_vcf true -profile singularity
```

### Running the Pipeline on Actual Data

To test the pipeline with actual datasets:

```bash
nextflow run build_reference_actual.nf -c nextflow_ref.config --mode actual --only_fastqc_fastp false --mode test -profile singularity
```

To execute FastQC and Fastp separately in actual mode:

```bash
nextflow run build_reference_actual.nf -c nextflow_ref.config --mode actual --only_fastqc_fastp true --mode test -profile singularity
```

To run the main pipeline on actual data:

```bash
nextflow run main.nf -c nextflow_main.config --mode actual --merge_vcf true -profile singularity
```

To annotate individual VCF files without merging:

```bash
nextflow run main.nf -c nextflow_main.config --mode actual --merge_vcf false -profile singularity
```

### Adjusting Parameters

Users can adjust the default filtration parameters by editing the `params.config` file. Additionally, users can select the desired SnpEff database for annotation.

### Data Handling

- The files are fetched from the respective folders (`data<test/actual>`).

## Notes

- The pipeline dynamically checks for necessary files and tools, downloading or preparing them only when missing.
- Modular design ensures flexibility and ease of adding new steps or tools.

---

Feel free to reach out for support or contributions!

