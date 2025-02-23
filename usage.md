# **RNA-Seq Variant Calling and RNA Fusion Detection Pipeline**  

## **Overview**  

This pipeline is designed for **RNA-seq variant calling** and **RNA fusion detection**, providing a comprehensive and automated solution for transcriptomic analysis. By leveraging **Nextflow DSL 2**, the pipeline is modular, scalable, and reproducible, making it easy to maintain and adapt to different datasets.  

### **Key Features**  
**Quality Control**: Ensures high-quality sequencing data through FastQC and Fastp.  
**Alignment & Processing**: Uses STAR for alignment and GATK for variant processing.  
**Variant Calling & Annotation**: Detects genetic variants and annotates them with SnpEff and VEP.  
**RNA Fusion Detection**: Identifies gene fusions using Arriba, crucial for cancer genomics.  
**Comprehensive Reporting**: Generates MultiQC reports and structured variant/fusion outputs.  

## **Pipeline Objectives**  

- **Genetic Variant Identification**: Detect and annotate genetic variants from RNA-seq data.  
- **RNA Fusion Detection**: Identify gene fusion events critical in cancer and other diseases.  
- **Automated & Structured Reporting**: Generate quality control metrics, structured variant reports, and fusion visualizations.  

---

# **Pipeline Components**  

This pipeline is divided into key **subworkflows**, each responsible for a specific aspect of data processing.  

## **Preprocessing**  

The **Preprocessing subworkflow** ensures that raw sequencing reads are cleaned, quality-checked, and prepared for downstream analysis.  

### **Steps:**  
1**Read Sample Metadata**: Extracts sample IDs, FASTQ paths, and strandedness from the user-provided `samplesheet.csv`.  
2️**FASTQ File Concatenation**: If enabled (`params.concatenate = true`), multiple FASTQ files per sample are merged.  
3️**Quality Control**:  
   - **FastQC** generates quality reports on raw sequencing reads.  
   - **Fastp** trims low-quality bases and removes adapter sequences.  
4️**MultiQC Aggregation**: Compiles FastQC and Fastp reports into a single quality summary.  

---

## **Variant Calling**  

The **VARIANT_CALLING** subworkflow is responsible for detecting genetic variants (SNPs & Indels) from aligned RNA-seq data.  

### **Steps:**  
1️**Alignment & BAM Processing**:  
   - Reads are aligned using **STAR**.  
   - BAM files are sorted and orphan reads are removed.  
   - Duplicate reads are marked using **GATK MarkDuplicates**.  
   - Base Quality Score Recalibration (BQSR) enhances accuracy.  

2️**Variant Calling & Filtering**:  
   - Variants are detected using **GATK HaplotypeCaller**.  
   - Filters are applied using **GATK Variant Filter**.  
   - Variant statistics are generated with **BCFtools**.  

3️**Variant Annotation**:  
   - Variants are annotated using **SnpEff** and **VEP**.  
   - Structured **CSV reports** are created for downstream analysis.  

### **Outputs:**  
 Annotated **VCF files**  
 **CSV summary reports** with filtered variant data  
 MultiQC metrics for alignment and variant calling  

---

## **3️Gene Fusion Detection**  

The **GENE_FUSION** subworkflow identifies **gene fusion events** from RNA-seq data using **Arriba**, a leading fusion detection tool.  

### **Steps:**  
1️**STAR Fusion Alignment**: Aligns RNA-seq reads in fusion detection mode.  
2️**Fusion Detection with Arriba**: Identifies fusion events using known databases and blacklist filtering.  
3️**Fusion Visualization**: Generates fusion diagrams for result interpretation.  

### **Outputs:**  
**Fusion results** (list of detected gene fusions).  
**Discarded fusions** (filtered potential false positives).  
**Graphical fusion visualizations** for easier interpretation.  

---

## **4️Reference Preparation (`build_references` Subworkflow)**  

This pipeline supports **automated reference preparation**, ensuring all required genome and annotation files are correctly downloaded or processed.  

### **Reference Downloads:**  
**Reference Genome** (FASTA)  
**SNP & Indel Variants** (VCF)  
**Blacklist Regions** (BED)  
**Genome Annotation** (GTF)  

### **Reference Indexing:**  
**FASTA Indexing** (`samtools faidx`)  
**Genome Dictionary** (`picard CreateSequenceDictionary`)  
**STAR Indexing** (for alignment)  
**VCF Merging & Filtering** (for known variants)  

### **Tools & Databases Setup:**  
**SnpEff & Database Download** (for variant annotation)  
**Arriba Tool & Known Fusions DB**  
**VEP Cache Preparation**  
**ClinVar Database Download**  

---

# **Workflow Execution**  

## **1️Sample Sheet Format**  

The pipeline requires a **sample sheet** to specify RNA-seq input files.  

### **Example Format:**  

| sample_id  | fastq_1                          | fastq_2                          | strandedness  |  
|------------|----------------------------------|----------------------------------|--------------|  
| Sample_01  | /path/to/sample_1_R1.fastq.gz   | /path/to/sample_1_R2.fastq.gz   | forward      |  
| Sample_02  | /path/to/sample_2_R1.fastq.gz   | /path/to/sample_2_R2.fastq.gz   | reverse      |  

- **`sample_id`**: Unique identifier for each sample.  
- **`fastq_1` / `fastq_2`**: Paths to the paired-end sequencing reads.  
- **`strandedness`**: Specifies strand orientation (`forward`, `reverse`, or `unstranded`).  

---

## **2️Running the Pipeline**  

### **A. Reference Preparation (If Required)**  

Run the reference setup process **before** executing the main pipeline:  

```bash
nextflow run build_reference_main.nf -c nextflow_ref_main.config -profile singularity
```  

---

### **B. Execute the Main Pipeline**  

Once references are ready, launch the full RNA-seq analysis pipeline:  

```bash
nextflow run main.nf -c nextflow_main.config -profile singularity
```  

**Adjustable Parameters:**  
- **`merge_vcf: true`** → Merge all VCFs into a single file for annotation.  
- **`only_variant_calling: true`** → Runs only variant calling workflow.  
- **`only_fusion_detection: true`** → Runs only fusion detection workflow.  

Example command to **run only variant calling**:  

```bash
nextflow run main.nf -c nextflow_main.config -profile singularity --only_variant_calling true
```  

---

# **Final Outputs**  

**MultiQC Reports** (Quality control overview)  
**Variant Annotation CSV Reports** (Structured variant data)  
**Fusion Detection Results** (List & visualizations)  
**Final Processed VCF Files** (for downstream analysis)  

---

# **Conclusion**  

This **Nextflow-based RNA-Seq pipeline** provides a **scalable, reproducible, and modular** approach to variant calling and gene fusion detection. By integrating state-of-the-art tools like **GATK, SnpEff, and Arriba**, it ensures high accuracy and efficiency for transcriptomic studies.  

 

---
