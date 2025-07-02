# RNA-Seq Variant Calling & Fusion Detection Pipeline

**A scalable, reproducible, and modular Nextflow pipeline for RNA-seq variant calling and fusion detection.**

**GitHub Repository:** [Variantcalling-and-Genefusion](https://github.com/kothaiRaja/Variantcalling-and-Genefusion.git)

---

## RNA-seq Variant Calling and Gene Fusion Workflow

This flowchart represents the RNA-seq analysis workflow, including variant calling and gene fusion detection.

![RNA-seq Pipeline Workflow](Documentation/complete%20pipeline.png)


---

## Overview

This pipeline processes RNA-seq data to:

* **Identify genetic variants (SNPs & Indels)** from transcriptomic data
* **Detect RNA fusion events** critical in cancer research
* **Perform quality control, alignment, annotation, and reporting** in an automated workflow

### Built With

* **Nextflow DSL2** for modular, scalable design
* **Singularity** for containerized, reproducible runs
* Tools: **FastQC**, **Fastp**, **STAR**, **GATK**, **SnpEff**, **VEP**, **Arriba**, **vcf2maf**, **maftools**, **MultiQC**

### Key Features

* Preprocessing: FastQC, Fastp, MultiQC
* Variant Calling: STAR, GATK HaplotypeCaller, filtering, annotation with SnpEff & VEP
* Fusion Detection: Arriba + visualization
* Modular Design: Easy extension and reuse of workflows and subworkflows
* Optional MAF reporting via vcf2maf and maftools
* Flexible configuration using custom and profile-based settings

---

## Workflow Structure

### 1. Preprocessing

* **FastQC & Fastp**: Raw read quality check and trimming
* **MultiQC**: Aggregates QC metrics for an overview

### 2. STAR Alignment

* **Two-pass STAR alignment** with splice junctions
* Filters orphan reads and computes **alignment statistics**

### 3. BAM Processing

* **Sorting**, **marking duplicates**, and **SplitNCigarReads**
* **Merge BAMs** (if split by intervals)
* Reset read groups and run **samtools calmd** for MD/NR tags

### 4. Base Recalibration

* **GATK BaseRecalibrator** on known variants
* **ApplyBQSR** to adjust base quality scores

### 5. Variant Calling

* **GATK HaplotypeCaller** across scattered intervals
* Merge VCFs and filter for high-confidence SNPs/Indels
* Compress, index, and generate variant statistics

### 6. Variant Annotation

* **SnpEff** and/or **VEP** (with plugins: LoF, CADD, REVEL, etc.)
* Convert VCFs to **MAF** using `vcf2maf`
* Generate mutation summary plots via **maftools**

### 7. Fusion Detection (Optional)

* **Arriba** identifies gene fusions from STAR-aligned BAMs
* Automatically produces **visual plots** of fusions

### 8. Reports

* **MultiQC** for preprocessing and alignment metrics
* **Arriba fusion plots** and **MAF visualizations**
* Software versions and logs for reproducibility

---

## Installation & Setup

### 1. Install Nextflow & Dependencies

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

Make sure **Singularity** or **Docker** is installed.

### 2. Clone the Repository

```bash
git clone https://github.com/kothaiRaja/Variantcalling-and-Genefusion.git
cd Variantcalling-and-Genefusion
```

---

## Running the Pipeline

### Step 1: Test Configuration(default)

 `test.config` inside the `test_data/` folder with parameters like:

```groovy
params {
  vep_enable = false
  annotation_tools = ['snpeff']
  samplesheet = "${baseDir}/test_data/sample_sheet_new.csv"
  arriba_input_dir = "${baseDir}/test_data/sample_ARRIBA"
  run_fusion = true
  maftools = false
  ref_base = "${baseDir}/test_data/new_test"
}
```

### Step 2: Run with Test Data

```bash
nextflow run main.nf \
  -c test_data/test.config \
  --ref_base /path/to/test_data/new_test \
  --resultsdir /path/to/results \
  -profile test,singularity
```



### Step 3: Run with Real Data

Before running the pipeline, create your own `custom.config` file.  
This file should include paths to available reference files or URLs, required databases, `workDir`, and `ref_base`.

```bash
nextflow run main.nf \
  -c path/to/custom.config \
  --resultsdir path/to/final_results_directory \
  -profile singularity
  
```

### Step 3a: Run with Real Data with reference files saving (also can be given in CLI)
```bash
nextflow run main.nf \
  -c path/to/custom.config \
  --ref_base /path/to/test_data/folder \
  --samplesheet path/to/sample_sheet.csv \
  --resultsdir path/to/final_results_directory \
  -profile singularity
```

---

## Input Format: Sample Sheet

Prepare a `samplesheet.csv`:

| sample\_id | fastq\_1                        | fastq\_2                        | strandedness |
| ---------- | ------------------------------- | ------------------------------- | ------------ |
| Sample\_01 | /path/to/sample\_1\_R1.fastq.gz | /path/to/sample\_1\_R2.fastq.gz | forward      |
| Sample\_02 | /path/to/sample\_2\_R1.fastq.gz | /path/to/sample\_2\_R2.fastq.gz | reverse      |

---

## Output Summary

Detailed output directories configured in `params.config`:

### Quality Control

* `results/multiqc_input/`: FastQC, Fastp, MultiQC inputs
* `results/multiqc_quality/`: MultiQC summary HTML report

### BAM Processing & Alignment

* `cache/sorted_bam/`, `filtered_bam/`, `merged_bam/`, `split_ncigar/`, `calmd/`, `recalibrated_bams/`: intermediate BAM steps
* `results/multiqc_input/`: STAR logs, flagstats, duplication metrics

### Variant Calling

* `cache/haplotype_caller/`, `merged_vcf/`, `variant_filter/`, `selected_variants*/`: VCF outputs
* `results/multiqc_input/`: BCFtools stats/query
* `results/vcf2table/`: Table of variants

### Variant Annotation

* `results/multiqc_input/`: Annotated VCFs (SnpEff and/or VEP)

### Fusion Detection

* `results/ARRIBA/`: Fusions TSV + discarded TSV
* `results/ARRIBA_visualisation/`: Fusion gene plots

### MAF Reporting

* `results/maftools/vcf2maf/`: MAF conversion output
* `results/maftools/visualisation/`: Mutation plot summary

### Pipeline Reports

* `results/timeline.html`, `dag.png`, `trace.txt`, `report.html`: Nextflow logs and runtime tracking

---

## Configuration Parameters

| Parameter          | Description                                            |
| ------------------ | ------------------------------------------------------ |
| `samplesheet`      | Path to CSV with sample information                    |
| `resultsdir`       | Directory for result summaries and final outputs       |
| `ref_base`         | Base path for data, cache, and reference lookup        |
| `run_fusion`       | Run Arriba fusion detection along with variant calling |
| `maftools`         | Enable MAF visual reporting using vcf2maf + maftools   |
| `annotation_tools` | List of tools to use: \['snpeff'], \['vep'], or both   |

---


