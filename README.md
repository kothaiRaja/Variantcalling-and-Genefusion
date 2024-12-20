Pipeline Implementation
1. Working with the pipeline: 
1.1 Input CSV Files: 
This pipeline is compatible with paired-end RNA-Seq data, and the input CSV file should contain all the necessary details about the samples you are processing. The file should include the sample identifiers, paths to the raw FASTQ files (both read1 and read2). Below is the structure of the CSV file used in the pipeline.
sample_id,fastq_1,fastq_2
sample_1,/home/XXX/cq-git-sample/Internship/data/test/fastp/filtered_sample_1_R1.fastq.gz,/home/XXX/cq-git-sample/Internship/data/test/fastp/filtered_sample_1_R2.fastq.gz
sample_2,/home/XXX/cq-git-sample/Internship/data/test/fastp/filtered_sample_2_R1.fastq.gz,/home/XXX/cq-git-sample/Internship/data/test/fastp/filtered_sample_2_R2.fastq.gz

1.2 Reading the Input CSV File
In the pipeline, the input CSV file is read using Nextflow’s built-in functionality. The file path is passed to the pipeline as a parameter (params.csv_file), and the data from the CSV file is loaded into the pipeline for further processing.
1.3 Reference files:
Running build_reference.nf pipeline downloads all the reference files and prepares all the files required for downstream analysis. Path to these files are also updated in nextflow.config file.
Nextflow run build_reference_<test/actual>.nf -c nextflow_ref.config -mode <test/actual> -profile singularity
1.4 Quality control :
1. Purpose of Quality Control in buildreference.nf:
 The primary goal of incorporating quality control in buildreference.nf is to ensure that users can:
o	Assess the quality of raw sequencing reads before starting the main steps of the pipeline. 
o	Identify potential issues such as low-quality bases, adapter contamination, or other artifacts that could affect downstream analyses. 
Make adjustments based on the quality of the reads, such as trimming low-quality bases, filtering poor-quality reads, or adjusting parameters for better alignment and variant calling. 
2. Quality Control Steps in buildreference.nf :
The quality control steps are executed using two main tools: FastQC and Fastp. 
Step 1: Raw Read Quality Check with FastQC FastQC is run on the raw FASTQ files to evaluate the quality of the sequencing reads. It provides a detailed summary of several key quality metrics, including: Per base sequence quality: Measures the average quality score of bases at each position in the read. Per sequence GC content: Evaluates the GC content distribution. Per base N content: Analyzes the number of ambiguous bases (N) across the reads. Sequence duplication levels: Identifies duplicate reads that might indicate PCR bias. Adapter content: Detects the presence of adapter sequences in the reads. 
Step 2: Trimming and Quality Filtering with Fastp Fastp is used to trim the reads and filter out low-quality sequences or adapter contamination. It ensures that the reads passing through the pipeline are of high quality by: Trimming low-quality bases: Removing bases with low quality scores from the ends of the reads. Removing adapter sequences: Identifying and removing any adapter contamination from the reads. Quality filtering: Removing reads that are too short or have low overall quality scores. Fastp can also generate per-read and summary-level reports, which help users understand how much data was trimmed or filtered and the overall quality of the data after processing. 
3. Flow in buildreference.nf In the buildreference.nf pipeline, quality control is handled early in the process, before the main alignment and variant calling steps. This ensures that the quality of the input data is assessed and cleaned before any computationally intensive steps are performed. Input: Raw FASTQ files (fastq_1 and fastq_2). QC Step 1: FastQC is run to generate reports on the quality of the raw reads. QC Step 2: Fastp is applied to trim and filter the reads. Output: Cleaned and high-quality FASTQ files, ready for alignment and downstream analysis. 
4. Adjustments Based on QC Reports Once FastQC and Fastp processes have been run, the quality control reports will provide a comprehensive overview of the data quality. Based on these reports, adjustments can be made before proceeding with the main pipeline. Quality Cutoff: If a large proportion of reads have low quality, the user can decide to discard the problematic reads or adjust the trimming parameters to remove low-quality bases. Adapter Removal: If adapter contamination is detected in the FastQC reports, Fastp can be used to trim the adapters more aggressively. Read Length Adjustment: If many reads are too short after trimming, it might be necessary to filter them out or change the trimming parameters to maintain enough read length for downstream analyses. 
5. Benefits of Including QC in buildreference.nf Early Quality Check: The quality of the input reads is checked and adjusted before the main analysis, preventing errors from propagating in downstream steps. Improved Data Quality: By ensuring that only high-quality reads are used, the accuracy of alignment, variant calling, and fusion detection improves. User Control: Users can decide how to handle different quality issues, such as trimming, filtering, or discarding problematic reads. This is so big. Please make it short precise
Pipeline 2: main.nf
Purpose:
The main.nf pipeline is designed to perform downstream analyses on the prepared data. It takes the cleaned and processed files from the previous buildreference.nf pipeline and conducts variant calling, gene fusion detection, and variant annotation.
Command to Run:
Once all the required files are ready, the pipeline can be run with the following command:

nextflow run main.nf -c nextflow_main.config  -mode <test/actual> -profile singularity 

This command initiates the pipeline in "actual mode" using Singularity containers.
Output:
•	Final VCF File: The pipeline generates a final.vcf file containing the variants identified through the analysis.
•	Annotation HTML Report: An HTML file with all annotations, including the functional impacts of the variants and associated gene information.
•	Gene Fusion Output and visualisation as pdf: A fusion.tsv file containing detected gene fusions.
Steps Performed in main.nf:
1.	Alignment and Variant Calling: Using the cleaned FASTQ files, STAR is used for alignment, followed by variant calling (e.g., with GATK HaplotypeCaller).
2.	Fusion Detection: The pipeline includes steps for detecting gene fusions using tools like ARRIBA, based on the aligned data.
3.	Variant Annotation: Variants are annotated with additional information using tools like SnpEff and GATK Variant Annotator, adding functional predictions and gene associations.
4.	Final Reports: Outputs a VCF file for variants, an HTML file for annotations, and a TSV file for gene fusion details.
Key Outputs:
•	Final.vcf: Contains the identified variants (SNPs and indels).
•	Fusion.tsv: Provides a list of detected gene fusions.
•	Annotation HTML file: Contains a comprehensive set of annotations for the variants.
