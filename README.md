# RNA-Seq Variant Calling and Gene Fusion Detection Pipeline  

This repository contains a **Nextflow pipeline** designed for RNA-Seq analysis, including **variant calling**, **gene fusion detection**, and **variant annotation**. It consists of two key pipelines:  
1. **buildreference.nf**: Prepares references, performs quality control, and manages resources efficiently by building/downloading only if files are unavailable.
2. **main.nf**: Conducts downstream analyses such as variant calling, gene fusion detection, and variant annotation.  

![Pipeline Workflow](https://github.com/user-attachments/assets/c6e4f029-acc8-47db-aa0d-4928d50c8538)  

---

## Features  
- **Variant Calling**: Using GATK HaplotypeCaller.  
- **Gene Fusion Detection**: Utilizing ARRIBA for comprehensive fusion identification.  
- **Annotation**: Functional annotation with Ensembl VEP and SnpEff.  
- **Quality Control**: Implements FastQC and Fastp for data preprocessing.  
- **Scalability**: Modular design supporting test and actual datasets.  

---

## Setup Instructions  

### 1. Clone the Repository  
  
git clone https://github.com/kothaiRaja/Variantcalling-and-Genefusion.git  
 

### 2. Navigate to the Project Directory  
 
cd Variantcalling-and-Genefusion  
  

### 3. Run the Pipelines  

#### a. Build Reference Files  
  
nextflow run build_reference_<test/actual>.nf -c nextflow_ref.config --only_fastqc_fastp <false/true> -profile singularity  
 

#### b. Main Analysis  
 
nextflow run main.nf -c nextflow_main.config --merge_vcf <true/false> -mode <test/actual> -profile singularity  
  

---

## Pipeline Descriptions  

### Pipeline 1: buildreference.nf  

#### Purpose  
The `buildreference.nf` pipeline prepares reference files and ensures data quality through rigorous quality control, laying the foundation for downstream analysis.  

#### Key Features  
1. **Input Handling**  
   - Requires a CSV file listing sample IDs and paths to paired-end FASTQ files.  
   - Example CSV format:  
     csv  
     sample_id,fastq_1,fastq_2  
     sample_1,/path/to/sample_1_R1.fastq.gz,/path/to/sample_1_R2.fastq.gz  
     sample_2,/path/to/sample_2_R1.fastq.gz,/path/to/sample_2_R2.fastq.gz  
       

2. **Efficient Reference Preparation**  
   - Downloads and prepares all necessary reference files only if they are not available in the specified data directory  
   - Ensures no redundant operations, saving time and computational resources.  

3. **Flexible Quality Control**  
   - **FastQC**: Assesses raw read quality and detects potential issues.
   -**Fastp**: Trims low-quality bases and removes adapter contamination.
   -Allows running only FastQC and Fastp by setting the parameter --only_fastp_fastqc true.

#### Outputs  
- Cleaned and trimmed FASTQ files.  
- Quality reports for raw and processed reads.  
- Prepared reference files for downstream analysis.  

---

### Pipeline 2: main.nf  

#### Purpose  
The `main.nf` pipeline performs the core analysis, including alignment, variant calling, gene fusion detection, and variant annotation.  

#### Key Features  
1. **Alignment and Variant Calling**  
   - Aligns cleaned reads to the reference genome using STAR.  
   - Calls variants (SNPs and indels) using GATK HaplotypeCaller.  

2. **Gene Fusion Detection**  
   - Detects gene fusions using ARRIBA, generating detailed reports for downstream analysis.  

3. **Variant Annotation**  
   - Annotates variants using SnpEff and Ensembl VEP, incorporating functional predictions and gene information.  

4. **Comprehensive Reporting**  
   - Produces annotated VCF files, HTML reports, and gene fusion details in TSV format.  

 #### Outputs  
- **Variants**: Annotated VCF file (`final.vcf`) with detailed variant information.  
- **Gene Fusion**: TSV file (`fusion.tsv`) summarizing detected gene fusions.  
- **Annotation Report**: HTML file providing detailed functional and gene information for variants.  

---

## Benefits of This Pipeline  

- **Robust Quality Control**: Ensures only high-quality reads are used, improving accuracy.  
- **Comprehensive Workflow**: From reference preparation to functional annotations, the pipeline handles all steps of RNA-Seq analysis.  
- **Customizable Modes**: Supports "test" and "actual" modes for easy testing and production runs.  
- **Scalability**: Leveraging Nextflow, the pipelines are scalable across computational environments.  
- **Standalone QC Option**: Allows running only FastQC and Fastp with --only_fastp_fastqc.

---

## Get Started  

Follow these steps to clone the repository, set up the required files, and execute the pipelines for a complete RNA-Seq analysis workflow.  

