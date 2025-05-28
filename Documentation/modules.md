
# Chapter 1: Trimming and Quality Control using Fastp (`TRIM_READS` Process)

## Overview

The `TRIM_READS` process is the first step in the RNA-seq preprocessing pipeline. It uses the **[fastp](https://github.com/OpenGene/fastp)** tool, a fast and versatile FASTQ preprocessor, to clean paired-end sequencing data by:

- Removing adapter sequences
- Trimming low-quality bases
- Filtering short reads
- Producing quality control reports (HTML & JSON)

Cleaning raw reads before alignment improves downstream performance and accuracy in variant calling, gene expression quantification, and fusion detection.

---

## Inputs

The process takes in the following inputs as a tuple:

```nextflow
tuple val(sample_id), path(r1), path(r2), val(strandedness)
```

| Input         | Description                                                  |
|---------------|--------------------------------------------------------------|
| `sample_id`   | A string identifier for the sample                           |
| `r1`          | Path to Read 1 (forward) FASTQ file                          |
| `r2`          | Path to Read 2 (reverse) FASTQ file                          |
| `strandedness`| Indicates library strandedness; passed downstream as-is      |

---

## Outputs

The process emits two output tuples:

1. **Trimmed Reads**

```nextflow
tuple val(sample_id), path("trimmed_${sample_id}_R1.fastq.gz"), path("trimmed_${sample_id}_R2.fastq.gz"), val(strandedness)
```

- Gzipped FASTQ files after adapter trimming and quality filtering
- Retains sample ID and strandedness

2. **fastp Reports**

```nextflow
tuple val(sample_id), path("${sample_id}_fastp.html"), path("${sample_id}_fastp.json")
```

- `.html`: Human-readable visual QC report
- `.json`: Structured QC report for parsing and MultiQC aggregation

---

## Code Block: The `fastp` Command

```bash
fastp -i "${r1}" -I "${r2}" \
		-o "trimmed_${sample_id}_R1.fastq.gz" \
		-O "trimmed_${sample_id}_R2.fastq.gz" \
		${params.fastp_extra} \
		--html "${sample_id}_fastp.html" \
		--json "${sample_id}_fastp.json"
```

---

## Parameter Breakdown and Purpose

| Parameter                               | Description                                                             |
| --------------------------------------- | ----------------------------------------------------------------------- |
| `-i "${r1}"`                            | Input FASTQ file – Read 1 (forward reads)                               |
| `-I "${r2}"`                            | Input FASTQ file – Read 2 (reverse reads)                               |
| `-o "trimmed_${sample_id}_R1.fastq.gz"` | Output file for trimmed Read 1 FASTQ                                    |
| `-O "trimmed_${sample_id}_R2.fastq.gz"` | Output file for trimmed Read 2 FASTQ                                    |
| `${params.fastp_extra}`                 | **User-configurable space to add extra parameters**; passed from config |
| `--html "${sample_id}_fastp.html"`      | Generates an interactive HTML report for QC visualization               |
| `--json "${sample_id}_fastp.json"`      | Generates a structured JSON report used by MultiQC or other tools       |


---

## Output Directory

The following directive is used to specify where outputs should be written:

```nextflow
publishDir trim_reads_outdir
```

This is configured using the pipeline parameter:  
`params.trim_reads_outdir`

All output files (trimmed reads and reports) will be saved here.

---

## How to Inspect and Interpret the Output

### ✅ Trimmed FASTQ Files

You can inspect the trimmed FASTQ files with tools like:

```bash
zcat trimmed_sample_R1.fastq.gz | head -n 8
```

Check for:

- Proper read formatting
- Read lengths matching expectations
- Absence of adapter contamination

### 📊 fastp HTML Report

Open the `.html` file in a browser to see:

- Base quality profiles (before and after)
- Adapter content
- Read length distribution
- Overrepresented sequences
- Filter summary

### 📈 fastp JSON Report

Useful for automation and summary tools:

```bash
cat sample_fastp.json | jq '.summary'   # View summary section (requires jq)
```

Can be combined with **MultiQC** for batch-level quality summaries.

---

## Summary

The `TRIM_READS` process is a critical first step to ensure your RNA-seq data is clean and ready for reliable downstream analysis. By trimming adapters, removing low-quality bases, and filtering short reads, we improve the accuracy of alignments and downstream variant detection. The use of `fastp` makes this process efficient and well-documented via its reports.

---


# Chapter 2: Read Alignment with STAR (`STAR_ALIGNMENT` Process)

## Overview

The `STAR_ALIGNMENT` process uses the **[STAR aligner](https://github.com/alexdobin/STAR)**, a highly efficient and accurate aligner tailored for RNA-seq data. It aligns paired-end reads to a reference genome, handles splicing events, and outputs sorted BAM files, splice junctions, gene counts, and chimeric alignments (for fusion detection).

This process performs a **two-pass alignment strategy**:
1. **First pass**: Detect splice junctions from the data.
2. **Second pass**: Use detected junctions and the GTF file to improve alignment quality.

---

## Inputs

```nextflow
tuple val(sample_id), path(trimmed_r1), path(trimmed_r2), val(strandedness)
path star_index_dir
path gtf_file
```

| Input         | Description |
|---------------|-------------|
| `sample_id`   | Sample identifier |
| `trimmed_r1`  | Trimmed Read 1 FASTQ (gzip) |
| `trimmed_r2`  | Trimmed Read 2 FASTQ (gzip) |
| `strandedness`| Strand-specificity of library |
| `star_index_dir` | STAR genome index directory |
| `gtf_file`    | Gene annotation file in GTF format |

---

## Outputs

This process emits multiple files, organized into tuples:

| Output Label | Description |
|--------------|-------------|
| `bam`        | Sorted BAM file of aligned reads |
| `log_final`  | STAR's final summary log (alignment statistics) |
| `log_out`    | STAR standard log |
| `log_progress` | STAR real-time progress log |
| `chimeric_sam` | Detected chimeric reads (for gene fusions, optional) |
| `junctions`  | Splice junctions detected |
| `versions`   | YAML file with STAR version used |

---

## Key STAR Parameters Explained

### General Options

| Parameter | Why It’s Used |
|----------|----------------|
| `--runThreadN` | Multithreading for speed. Uses `task.cpus`. |
| `--readFilesIn` | Input trimmed paired-end reads. |
| `--readFilesCommand zcat` | Unzips `.fastq.gz` files on-the-fly. |
| `--genomeDir` | Location of pre-built STAR genome index. |
| `--limitBAMsortRAM` | Controls memory for BAM sorting. Dynamically set from task memory. |
| `--outSAMunmapped Within` | Keeps unmapped reads in the output BAM. |

### Splice Junction & Filtering

| Parameter | Why It’s Used |
|----------|----------------|
| `--outFilterType BySJout` | Filters alignments based on known splice junctions. |
| `--alignSJoverhangMin` | Minimum overhang for junction reads; improves splicing accuracy. |
| `--alignSJDBoverhangMin` | Overhang required when using the GTF annotation junctions. |
| `--outFilterMismatchNmax` | Max mismatches per read; set high (999) to avoid over-filtering. |
| `--outFilterMatchNmin` | Minimum number of matching bases. |
| `--outFilterMatchNminOverLread` / `--outFilterScoreMinOverLread` | Normalize filtering by read length. |

### Strandness Handling

If the library is **unstranded**, STAR requires:

```bash
--outSAMstrandField intronMotif
```

This enables correct processing of strand-agnostic RNA-seq data.

### Output Format and Attributes

| Parameter | Why It’s Used |
|----------|----------------|
| `--outSAMtype BAM SortedByCoordinate` | Outputs sorted BAM ready for downstream use. |
| `--outSAMattributes NH HI AS nM MD NM` | Adds tags useful for variant calling and QC. |
| `--outSAMattrRGline` | Adds read group info (sample, platform, center). Essential for GATK compatibility. |

---

## Two-Pass Strategy

### First Pass: Splice Junction Discovery

In the first pass, STAR aligns reads **without annotation** to discover novel junctions:

```bash
STAR ... --outFileNamePrefix ${sample_id}_pass1_
```

The key output is:  
`${sample_id}_pass1_SJ.out.tab`

### Second Pass: Refined Alignment with GTF

In the second pass, STAR uses:

- The junctions from pass 1 (`--sjdbFileChrStartEnd`)
- The GTF file (`--sjdbGTFfile`)  
to produce final alignments and counts:

```bash
STAR ... --outFileNamePrefix ${sample_id}_ ...
```

---

## Chimeric Detection (Gene Fusion Support)

STAR is configured to detect potential gene fusions:

| Parameter | Purpose |
|----------|---------|
| `--chimSegmentMin` | Minimum length of each chimeric segment. |
| `--chimJunctionOverhangMin` | Minimum overhang on each side of junction. |
| `--chimScoreMin` | Minimum total score for a chimeric alignment. |
| `--chimScoreJunctionNonGTAG` | Penalty for non-canonical junctions. |
| `--chimOutType WithinBAM HardClip` | Outputs chimeric reads into BAM; helpful for Arriba. |

---

## Quantification Options

STAR can quantify gene expression simultaneously:

```bash
--quantMode TranscriptomeSAM GeneCounts
```

- `TranscriptomeSAM`: Used by tools like RSEM or Salmon for transcript quant.
- `GeneCounts`: STAR’s built-in gene-level counting.

---

## Extra Options and Version Capture

- `params.star_extra_args` allows pipeline users to pass any extra STAR arguments.
- At the end, the STAR version is recorded in a YAML file:

```yaml
"STAR_ALIGNMENT":
  STAR: "2.7.10a"
```

---

## Output Directory

Defined by:

```nextflow
publishDir params.star_outdir, mode: "copy"
```

This saves all STAR output files to the user-defined directory (`params.star_outdir`) and copies them to retain original names.

---

## How to Validate Outputs

### BAM File

```bash
samtools flagstat sample_Aligned.sortedByCoord.out.bam
```

- Confirms read mapping quality and rates.

### Logs

- `Log.final.out`: Summary statistics
- `Log.progress.out`: Real-time alignment stats (can be monitored during runtime)
- `Log.out`: General STAR log

### Splice Junctions

- Found in `${sample_id}_SJ.out.tab`
- Useful for checking novel junction discovery

### Gene Counts

- Found in `${sample_id}_ReadsPerGene.out.tab` (if `--quantMode` enabled)

---

## Summary

The `STAR_ALIGNMENT` process is a powerful, flexible step that maps high-quality RNA-seq reads to the genome while accounting for splicing, strandedness, and gene fusion events. With its two-pass alignment strategy and extensive logging, it provides high confidence alignments ready for variant calling, quantification, or downstream analyses like gene fusion detection.

---


# Chapter 3: Duplicate Marking with GATK (`GATK_MARK_DUPLICATES` Process)

## Overview

The `GATK_MARK_DUPLICATES` process uses the **GATK MarkDuplicates** tool to identify and mark duplicate reads in a sorted BAM file. Duplicate reads typically result from PCR amplification during library preparation. Removing or marking them ensures more accurate variant calling and downstream quantification.

This process also:
- Generates a BAM index file
- Outputs a duplicate metrics report
- Captures the GATK version used

---

## Inputs

```nextflow
tuple val(sample_id), val(strandedness), path(sorted_bam), path(sorted_bam_index)
```

| Input            | Description |
|------------------|-------------|
| `sample_id`      | Unique sample identifier |
| `strandedness`   | Strand information (passed forward) |
| `sorted_bam`     | Sorted BAM file (from STAR or other aligner) |
| `sorted_bam_index` | BAM index (.bai) corresponding to the sorted BAM |

---

## Outputs

| Output Tuple | Description |
|--------------|-------------|
| `marked_duplicates.bam`, `marked_duplicates.bai` | BAM file with duplicates marked and indexed |
| `dup_metrics.txt` | Duplicate marking metrics from GATK |
| `versions.yml` | GATK version metadata for reproducibility |

The outputs are copied to the user-defined directory specified in:

```nextflow
publishDir params.markduplicates_outdir, mode: "copy"
```

---

## Core Command and Parameters

```bash
gatk MarkDuplicates \
    -I ${sorted_bam} \
    -O ${sample_id}_marked_duplicates.bam \
    -M ${sample_id}_dup_metrics.txt \
    --CREATE_INDEX true \
    --REMOVE_DUPLICATES ${params.remove_duplicates ? 'true' : 'false'} \
    --VALIDATION_STRINGENCY ${params.validation_stringency ?: 'LENIENT'}
```

### Parameter Explanations

| Parameter | Purpose |
|----------|---------|
| `-I` | Input BAM file (sorted) |
| `-O` | Output BAM with marked duplicates |
| `-M` | Metrics file reporting duplicate stats |
| `--CREATE_INDEX true` | Creates `.bai` index for the output BAM |
| `--REMOVE_DUPLICATES` | Removes duplicates instead of just marking them; can be toggled with `params.remove_duplicates` |
| `--VALIDATION_STRINGENCY` | Controls error handling for malformed records. Set to `LENIENT` by default to avoid breaking on minor issues |

---

## Duplicate Metrics

The metrics file `${sample_id}_dup_metrics.txt` contains important statistics such as:

- **Number of reads examined**
- **Number of duplicates found**
- **Percent duplication**
- **Unmapped reads**

You can inspect this file manually or include it in MultiQC.

---

## Error Checking

The process includes safety checks to ensure that key output files are generated correctly:

```bash
if [ ! -s ${sample_id}_marked_duplicates.bam ]; then
    echo "Error: Marked duplicates BAM file not generated for ${sample_id}" >&2
    exit 1
fi
```

These checks prevent silent failures and ensure reliability in the workflow.

---

## Version Recording

For reproducibility, the process captures the GATK version and writes it to `versions.yml`:

```yaml
"GATK_MARK_DUPLICATES":
  gatk: "4.4.0.0"
```

This file is emitted with the `versions` label and can be collected later in the pipeline or included in final reports.

---

## How to Validate Output

### BAM File

```bash
samtools flagstat sample_marked_duplicates.bam
```

- Confirms number of reads, duplicates, and proper indexing

### Metrics File

- Look for `% duplication` and number of optical/PCR duplicates
- Useful to determine library complexity and PCR bias

---

## Summary

The `GATK_MARK_DUPLICATES` process ensures that PCR-induced read duplicates are properly handled before variant calling. By either marking or removing these duplicates and producing detailed metrics, this step contributes to a more accurate and interpretable analysis pipeline. All results are validated for completeness and include version-tracked metadata for reproducibility.


# Chapter 4: Creating and Scattering Genomic Intervals

This chapter describes two related processes used to define and manage genomic intervals for parallelized variant calling:

1. **`BED_TO_INTERVAL_LIST`** — Converts BED files into GATK-compatible `.interval_list` files.
2. **`SCATTER_INTERVAL_LIST`** — Splits interval lists into multiple chunks to enable scatter-gather parallel processing.

These steps are typically used in GATK workflows such as Base Quality Score Recalibration (BQSR), variant calling, and coverage calculations.

---

## Process: `BED_TO_INTERVAL_LIST`

### Purpose

Converts a BED file (standard for describing genomic regions) into a **GATK `.interval_list`** format. This is required by GATK tools like `BaseRecalibrator`, `HaplotypeCaller`, etc., which expect interval lists with reference contig metadata.

---

### Inputs

```nextflow
tuple val(meta_id), path(bed_file)
path(genome_fasta)
path(genome_dict)
```

| Input | Description |
|-------|-------------|
| `meta_id` | Sample or project-level ID |
| `bed_file` | BED file with regions of interest (e.g., exome capture kit) |
| `genome_fasta` | Reference genome FASTA file (not used directly here but implied) |
| `genome_dict` | Reference genome sequence dictionary `.dict` file required by GATK |

---

### Output

```nextflow
tuple val(meta_id), path("${bed_file.baseName}.interval_list")
```

- A GATK-compatible `.interval_list` file corresponding to the input BED file.

---

### Script Breakdown

```bash
gatk BedToIntervalList \
    -I ${bed_file} \
    -O ${bed_file.baseName}.interval_list \
    -SD ${genome_dict}
```

### Explanation

| Parameter | Purpose |
|----------|---------|
| `-I` | Input BED file |
| `-O` | Output interval list file |
| `-SD` | Sequence dictionary from the reference genome |

The sequence dictionary adds necessary metadata about chromosomes and contigs to the `.interval_list`.

---

## Process: `SCATTER_INTERVAL_LIST`

### Purpose

Splits a single `.interval_list` file into multiple **smaller interval files** to enable parallel processing across genomic regions — a key step in optimizing GATK tools for performance on large datasets.

---

### Inputs

```nextflow
tuple val(meta), path(interval_list)
path(genome_dict)
```

| Input | Description |
|-------|-------------|
| `meta` | Sample or task-level identifier |
| `interval_list` | GATK interval list from previous step |
| `genome_dict` | Reference genome dictionary file |

---

### Output

```nextflow
tuple val(meta), path("*.interval_list")
```

- A set of split `.interval_list` files, ready for parallelized GATK processing.

---

### Script Breakdown

```bash
gatk IntervalListTools \
    --INPUT ${interval_list} \
    --OUTPUT scattered_intervals \
    --SCATTER_COUNT ${params.scatter_count} \
    --UNIQUE true
```

### Explanation of Parameters

| Parameter | Purpose |
|----------|---------|
| `--INPUT` | Original interval list file to scatter |
| `--OUTPUT` | Output directory for scattered files |
| `--SCATTER_COUNT` | Number of intervals to create (set via pipeline param) |
| `--UNIQUE` | Ensures non-overlapping intervals in output |

The output files are then **renamed** to avoid filename collisions, and **moved back** into the working directory:

```bash
for f in scattered_intervals/*/*; do
    dir_name=$(basename $(dirname "$f"))
    file_name=$(basename "$f")
    mv "$f" "scattered_intervals/${dir_name}_${file_name}.interval_list"
done
mv scattered_intervals/*.interval_list .
rm -r scattered_intervals
```

This ensures consistent and organized output for downstream processes like `BaseRecalibrator` or `HaplotypeCaller`.

---

## Example Use Case

Let’s say you have a BED file for exome targets. You would:

1. Convert it to `.interval_list` using `BED_TO_INTERVAL_LIST`
2. Scatter it into 20 parts using `SCATTER_INTERVAL_LIST`
3. Run BQSR or variant calling on each part in parallel

This enables **faster computation** and **better resource management** in high-throughput pipelines.

---

## Summary

| Process | Goal | Output |
|--------|------|--------|
| `BED_TO_INTERVAL_LIST` | Convert BED to `.interval_list` | GATK-compatible intervals |
| `SCATTER_INTERVAL_LIST` | Split intervals into chunks | Multiple `.interval_list` files for parallelization |

These two steps provide the foundation for scalable, efficient execution of computationally intensive GATK processes.

---

# Chapter 5: Splitting Reads with N CIGARs (`SPLIT_NCIGAR_READS` Process)

## Overview

The `SPLIT_NCIGAR_READS` process prepares RNA-seq aligned BAM files for variant calling using **GATK’s `SplitNCigarReads`**. This step is critical for **handling RNA spliced reads**, where the CIGAR string contains `N` operations representing large introns.

Variant callers like `HaplotypeCaller` do not natively support reads with `N` operators in CIGAR strings. Hence, `SplitNCigarReads`:

- Splits reads at N CIGAR operations
- Adjusts mapping qualities
- Rewrites alignments to be compatible with GATK tools

---

## Inputs

```nextflow
tuple val(sample_id), val(strandedness), path(bam), path(bai), path(interval)
path genome_fasta
path index
path genome_dict
```

| Input | Description |
|-------|-------------|
| `sample_id` | Unique sample name |
| `strandedness` | Strand-specific info (passed forward) |
| `bam` | Aligned BAM file from STAR |
| `bai` | BAM index file |
| `interval` | A genomic interval file (one chunk from scattered intervals) |
| `genome_fasta` | Reference genome (FASTA format) |
| `index` | Index file for FASTA (e.g., `.fai`) |
| `genome_dict` | Sequence dictionary file for the genome |

---

## Outputs

```nextflow
tuple val(sample_id), val(strandedness), path("${sample_id}_split_${interval.baseName}.bam"), path("${sample_id}_split_${interval.baseName}.bai")
```

| Output File | Description |
|-------------|-------------|
| `.bam` | BAM file with reads split at N CIGARs |
| `.bai` | BAM index for the above file |

---

## Why This Step Is Important

RNA-seq reads span exons and introns. When mapped with STAR, the intronic regions are represented with `N` in the CIGAR string (e.g., `76M100N24M`). GATK tools like `HaplotypeCaller` are not designed to interpret these `N` operations correctly. This can lead to missed variants or errors during execution.

`SplitNCigarReads` preprocesses these reads to:

- Break spliced reads at junctions (`N`)
- Adjust MAPQ scores if needed
- Remove or reprocess supplementary alignments
- Create BAMs suitable for GATK variant calling

---

## Core Command Breakdown

```bash
gatk --java-options "-Xmx${avail_mem}g" SplitNCigarReads \
    -R ${genome_fasta} \
    -I ${bam} \
    -O ${sample_id}_split_${interval.baseName}.bam \
    --skip-mapping-quality-transform false \
    --max-mismatches-in-overhang 1 \
    --max-bases-in-overhang 50 \
    --create-output-bam-index true \
    --process-secondary-alignments true \
    ${interval_command}
```

### Java Memory

```bash
--java-options "-Xmx${avail_mem}g"
```

- Dynamically allocates available memory from the Nextflow task specification.

---

## Key GATK Parameters

| Parameter | Purpose |
|----------|---------|
| `-R` | Reference genome (FASTA) |
| `-I` | Input BAM file with RNA-seq alignments |
| `-O` | Output BAM with split reads |
| `--skip-mapping-quality-transform false` | Retains STAR's default MAPQ=255, avoids setting MAPQ=60 |
| `--max-mismatches-in-overhang` | Limits mismatches in exon overhangs during realignment |
| `--max-bases-in-overhang` | Caps bases allowed in the junction overhang region |
| `--create-output-bam-index true` | Automatically creates `.bai` index |
| `--process-secondary-alignments true` | Processes secondary alignments (e.g., multi-mapped reads) |
| `--intervals` | Optional: restricts processing to one genomic region (if set) |

---

## Parallelization with Interval Lists

This process is usually executed **per-interval**, where each interval is processed separately for speed and scalability. The `interval.baseName` is used to uniquely name each output chunk.

Example interval name: `chr1_0_1000000`  
Output: `sample1_split_chr1_0_1000000.bam`

---

## Logging and Debugging

At runtime, you’ll see messages like:

```
Running SplitNCigarReads for sample: sample1 on interval: chr1_0_1000000
...
SplitNCigarReads completed for sample: sample1 on interval: chr1_0_1000000
```

If any interval fails, logs can help isolate the issue by interval and sample.

---

## Validation and Inspection

You can validate the split BAM using:

```bash
samtools view sample_split_chr1_0_1000000.bam | head
samtools flagstat sample_split_chr1_0_1000000.bam
```

Check for:

- Reasonable read counts
- No broken headers
- Correct splitting at junctions

---

## Summary

The `SPLIT_NCIGAR_READS` step is essential for making spliced RNA-seq reads compatible with GATK. It transforms BAM files containing `N` CIGAR operations into a format that `HaplotypeCaller` and other tools can process without errors. The use of intervals allows for distributed computing, making this step scalable for large datasets.

---

# Chapter 6: Base Quality Score Recalibration (`GATK_BASERECALIBRATOR` Process)

## Overview

The `GATK_BASERECALIBRATOR` process performs **Base Quality Score Recalibration (BQSR)** using the **GATK `BaseRecalibrator`** tool. BQSR models systematic errors made by the sequencing machine when it estimates the quality of each base call. These recalibrated scores improve the accuracy of downstream **variant calling**.

This step generates a recalibration table (`.recal_data.table`) that will later be applied to adjust the quality scores in the BAM file.

---

## Why BQSR Is Important

Sequencing machines tend to produce **systematic errors** in quality scores, which can mislead variant callers. GATK BQSR uses a **machine learning model** to detect and correct these errors by comparing observed base calls to **known, trusted variant sites** (e.g., dbSNP).

Corrected quality scores lead to:
- Fewer false positives
- Higher precision and recall in variant calling
- Better calibration for downstream filtering

---

## Inputs

```nextflow
tuple val(sample_id), val(strandedness), path(bam), path(bai), path(interval)
path(genome_fasta)
path(index)
path(dict)
path(known_variants)
path(known_variants_index)
```

| Input | Description |
|-------|-------------|
| `sample_id` | Sample identifier |
| `strandedness` | Strand-specificity info |
| `bam`, `bai` | Split BAM and index from `SplitNCigarReads` |
| `interval` | Genomic interval for regional recalibration |
| `genome_fasta`, `index`, `dict` | Reference genome, index (`.fai`), and dictionary (`.dict`) |
| `known_variants` | VCF file of known variant sites (e.g., dbSNP, Mills) |
| `known_variants_index` | Index file for the known variants VCF |

---

## Outputs

```nextflow
tuple val(sample_id), val(strandedness), path("${sample_id}_recal_data.table")
```

- A text file containing recalibration data for the BAM file
- Used in the `ApplyBQSR` step to adjust base quality scores

---

## Core Command

```bash
gatk BaseRecalibrator \
    -R ${genome_fasta} \
    -I ${bam} \
    --known-sites ${known_variants} \
    -O ${sample_id}_recal_data.table \
    ${interval_command}
```

### Parameter Explanation

| Parameter | Purpose |
|----------|---------|
| `-R` | Reference genome in FASTA format |
| `-I` | Input BAM file with split CIGARs |
| `--known-sites` | Trusted variant sites used to model sequencing errors |
| `-O` | Output recalibration table |
| `--intervals` | Optional: limit recalibration to a specific interval |

> Intervals allow recalibration to be run in parallel chunks for performance.

---

## Quality Control Check

The process includes a file check after the run:

```bash
if [ ! -s ${sample_id}_recal_data.table ]; then
    echo "Error: Recalibration table not generated for ${sample_id}" >&2
    exit 1
fi
```

This ensures only successful recalibration tables are passed forward.

---

## Optional Intervals

When the process is part of a scatter/gather workflow, it supports per-interval recalibration via:

```bash
--intervals ${interval}
```

This allows the full BAM to be processed in smaller genomic chunks in parallel — speeding up workflows significantly.

---

## Validation Tips

You can inspect the output `.table` file for calibration metrics:

```bash
head sample_recal_data.table
```

Look for:
- Recalibration parameters
- Covariate bins
- Quality score adjustments

GATK also supports plotting these results (`AnalyzeCovariates`), though it's typically used in research or benchmarking, not production workflows.

---

## Summary

The `GATK_BASERECALIBRATOR` process generates a recalibration model that corrects systematic errors in base quality scores using a trusted set of known variants. This improves the sensitivity and specificity of downstream variant calling. It is a required step in best-practices GATK RNA-seq workflows and is fully parallelizable via intervals for performance.

---

# Chapter 7: Applying Base Quality Score Recalibration (`GATK_APPLYBQSR` Process)

## Overview

The `GATK_APPLYBQSR` process uses the **GATK `ApplyBQSR`** tool to apply base quality score corrections to RNA-seq aligned BAM files, using a recalibration table previously generated by `BaseRecalibrator`.

This step adjusts the quality scores in the BAM file based on the statistical model learned in the previous step, improving the accuracy of downstream **variant calling**.

---

## Why ApplyBQSR Is Important

After modeling quality score errors with known variants (via `BaseRecalibrator`), we need to **apply** that model to the actual sequencing reads. This ensures:

- Corrected base quality scores
- Reduced systematic bias
- Improved performance of GATK tools like `HaplotypeCaller`

---

## Inputs

```nextflow
tuple val(sample_id), val(strandedness), path(bam), path(bai), path(recal_table), path(interval)
path(genome_fasta)
path(index)
path(dict)
```

| Input | Description |
|-------|-------------|
| `sample_id` | Sample name |
| `strandedness` | Strand information (passed through) |
| `bam`, `bai` | Input BAM and index file from `SplitNCigarReads` |
| `recal_table` | Recalibration table generated by `BaseRecalibrator` |
| `interval` | Genomic interval (for scatter-parallelization) |
| `genome_fasta` | Reference genome FASTA |
| `index` | FASTA index (.fai) |
| `dict` | Sequence dictionary file for reference |

---

## Outputs

```nextflow
tuple val(sample_id), val(strandedness),
      path("${sample_id}_${interval.baseName}_recalibrated.bam"),
      path("${sample_id}_${interval.baseName}_recalibrated.bai")
```

- A **recalibrated BAM file** per interval
- A **corresponding BAM index** file

These are saved to the directory set in:

```nextflow
publishDir params.recalibrated_bams_outdir, mode: "copy"
```

---

## Core Command

```bash
gatk ApplyBQSR \
    -R ${genome_fasta} \
    -I ${bam} \
    --bqsr-recal-file ${recal_table} \
    -O ${sample_id}_${interval.baseName}_recalibrated.bam \
    --intervals ${interval}
```

### Parameter Explanations

| Parameter | Purpose |
|----------|---------|
| `-R` | Reference genome in FASTA format |
| `-I` | Input BAM file to recalibrate |
| `--bqsr-recal-file` | Recalibration table from `BaseRecalibrator` |
| `-O` | Output recalibrated BAM |
| `--intervals` | Restricts recalibration to a genomic region (if applicable) |

---

## Additional Step: BAM Indexing

After recalibration, the BAM file is indexed using:

```bash
gatk BuildBamIndex -I recalibrated.bam
```

This ensures the recalibrated BAM is queryable and ready for downstream analysis like variant calling.

---

## Error Checking

The process includes safety checks:

```bash
if [ ! -s ${sample_id}_${interval.baseName}_recalibrated.bam ]; then
    echo "Error: Recalibrated BAM not generated for ${sample_id}" >&2
    exit 1
fi
```

This prevents silent failures from affecting the rest of the pipeline.

---

## Validation and QC

### Confirm BAM Generation

```bash
samtools quickcheck -v sample_recalibrated.bam
```

### Check BAM Statistics

```bash
samtools flagstat sample_recalibrated.bam
```

This helps confirm read counts and successful processing.

---

## Interval Support

This step supports interval-based parallelism, just like `SplitNCigarReads` and `BaseRecalibrator`, enabling scatter/gather execution across genomic regions.

Each output BAM is named using the interval, e.g.:

```
sample1_chr1_0_1000000_recalibrated.bam
```

---

## Summary

The `GATK_APPLYBQSR` step applies the learned recalibration model to correct base quality scores in RNA-seq alignments. This ensures that downstream tools like `HaplotypeCaller` work with more accurate base quality data, reducing bias and improving the confidence of variant calls.

---

# Chapter 8: Variant Calling with GATK HaplotypeCaller (`GATK_HAPLOTYPE_CALLER` Process)

## Overview

The `GATK_HAPLOTYPE_CALLER` process is the final step in the variant detection portion of your RNA-seq pipeline. It uses **GATK's `HaplotypeCaller`**, configured with RNA-seq specific flags, to call SNPs and indels from **recalibrated RNA-seq BAM files**.

This process generates a **compressed VCF file** (`.vcf.gz`) and its index (`.vcf.gz.tbi`), containing variant calls for a given sample.

---

## Why Use HaplotypeCaller for RNA-seq?

Although designed primarily for DNA sequencing, `HaplotypeCaller` supports RNA-seq variant calling with a few key adaptations:

- Spliced reads (handled via `SplitNCigarReads`)
- Strand-specific base recalibration
- Avoidance of soft-clipped bases
- Disabling duplicate read filtering (handled earlier)

These adjustments allow for **high-confidence variant detection** even with the complexity of RNA data.

---

## Inputs

```nextflow
tuple val(sample_id), val(strandedness), path(bam), path(bai)
path(genome)
path(genome_index)
path(genome_dict)
path(known_sites_vcf)
path(known_sites_vcf_index)
```

| Input | Description |
|-------|-------------|
| `sample_id` | Sample name |
| `strandedness` | Library strand info (carried forward) |
| `bam`, `bai` | Recalibrated BAM file and index |
| `genome`, `genome_index`, `genome_dict` | Reference genome FASTA, `.fai`, and `.dict` |
| `known_sites_vcf`, `known_sites_vcf_index` | dbSNP or other known variants for annotation |

---

## Outputs

```nextflow
tuple val(sample_id), val(strandedness), path("output_${bam.baseName}.vcf.gz"), path("output_${bam.baseName}.vcf.gz.tbi")
```

| Output File | Description |
|-------------|-------------|
| `.vcf.gz` | Gzipped variant call file |
| `.vcf.gz.tbi` | Tabix index for VCF (used for downstream tools and browsing) |

---

## Core Command

```bash
gatk HaplotypeCaller \
    --native-pair-hmm-threads ${task.cpus} \
    --reference ${genome} \
    --output output_${bam.baseName}.vcf.gz \
    -I ${bam} \
    --standard-min-confidence-threshold-for-calling 10.0 \
    --min-base-quality-score 10 \
    --output-mode EMIT_VARIANTS_ONLY \
    --dont-use-soft-clipped-bases true \
    --disable-read-filter NotDuplicateReadFilter \
    --dbsnp ${known_sites_vcf} \
    --verbosity INFO
```

### RNA-seq Specific Flags

| Parameter | Purpose |
|----------|---------|
| `--dont-use-soft-clipped-bases true` | Ignores bases soft-clipped by the aligner (can be misleading) |
| `--disable-read-filter NotDuplicateReadFilter` | Avoids double-filtering already marked duplicates |
| `--min-base-quality-score 10` | Minimum base quality to consider in variant calling |
| `--standard-min-confidence-threshold-for-calling 10.0` | Minimum confidence to emit a variant |
| `--output-mode EMIT_VARIANTS_ONLY` | Outputs only variant sites (not reference-matching positions) |
| `--dbsnp` | Annotates known variants for reference and comparison |

### Performance Settings

| Parameter | Purpose |
|----------|---------|
| `--native-pair-hmm-threads` | Enables multithreaded variant calling using the native HMM implementation |

---

## Error Checking

The script includes checks to ensure that critical files exist before proceeding:

```bash
if [ ! -s ${bam} ]; then
    echo "Error: BAM file not found or empty." >&2
    exit 1
fi
if [ ! -s ${genome} ]; then
    echo "Error: Reference genome not found or empty." >&2
    exit 1
fi
```

---

## Output File Naming

The output files are named based on the BAM file, e.g.:

```
input BAM: sample1_chr1_0_1000000_recalibrated.bam
output VCF: output_sample1_chr1_0_1000000_recalibrated.vcf.gz
```

This ensures easy traceability between input BAM chunks and resulting variant calls.

---

## Validation and QC

### Validate VCF File

```bash
bcftools stats output_sample.vcf.gz
```

Check for:
- Number of SNPs and indels
- Transition/transversion ratio (Ti/Tv)
- QUAL score distributions

### Inspect First Few Lines

```bash
zcat output_sample.vcf.gz | head -n 20
```

Make sure VCF headers are present and variants are properly annotated.

---

## Summary

The `GATK_HAPLOTYPE_CALLER` step completes the RNA-seq variant calling workflow by identifying high-confidence variants in spliced transcriptomic data. With RNA-specific tuning and integration of known variants, this process provides high-quality VCFs suitable for downstream analysis, annotation, or filtering.

---

# Chapter 9: Variant Filtering with GATK (`GATK_VARIANT_FILTER` Process)

## Overview

The `GATK_VARIANT_FILTER` process uses **GATK's `VariantFiltration`** tool to apply hard filters to variants based on their quality metrics. These filters help remove low-confidence calls and technical artifacts, improving the reliability of the final VCF.

This step is **crucial before downstream interpretation**, annotation, or biological analysis.

---

## Inputs

```nextflow
tuple val(sample_id), path(vcf_file), path(vcf_index)
path genome
path genome_index
path genome_dict
```

| Input | Description |
|-------|-------------|
| `sample_id` | Unique sample identifier |
| `vcf_file`, `vcf_index` | Raw VCF and its Tabix index from `HaplotypeCaller` |
| `genome`, `genome_index`, `genome_dict` | Reference FASTA, `.fai` index, and `.dict` |

---

## Outputs

```nextflow
tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi")
```

| Output File | Description |
|-------------|-------------|
| `.vcf.gz` | Gzipped filtered VCF file |
| `.vcf.gz.tbi` | Tabix index for the filtered VCF |

All outputs are saved to the directory defined by:

```nextflow
publishDir params.variant_filter_outdir, mode: "copy"
```

---

## Core Command

```bash
gatk VariantFiltration \
    -R ${genome} \
    -V ${vcf_file} \
    --cluster-window-size ${params.gatk_vf_window_size} \
    --cluster-size ${params.gatk_vf_cluster_size} \
    --filter-name "LowQual" --filter-expression "QUAL < ${params.gatk_vf_qual_filter}" \
    --filter-name "LowQD" --filter-expression "QD < ${params.gatk_vf_qd_filter}" \
    --filter-name "HighFS" --filter-expression "FS > ${params.gatk_vf_fs_filter}" \
    --filter-name "LowMQ" --filter-expression "MQ < ${params.gatk_vf_mq_filter}" \
    --filter-name "HighSOR" --filter-expression "SOR > ${params.gatk_vf_sor_filter}" \
    --filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum < ${params.gatk_vf_read_pos_filter}" \
    --filter-name "LowBaseQRankSum" --filter-expression "BaseQRankSum < ${params.gatk_vf_baseq_filter}" \
    -O ${sample_id}_filtered.vcf.gz
```

---

## Parameter Explanations

| Filter | Expression | Purpose |
|--------|------------|---------|
| `LowQual` | `QUAL < ...` | Removes variants with overall low quality score |
| `LowQD` | `QD < ...` | QD = Quality / Depth; low QD suggests low signal-to-noise |
| `HighFS` | `FS > ...` | FS = Fisher Strand Bias; high FS indicates strand-specific errors |
| `LowMQ` | `MQ < ...` | Low mapping quality can indicate alignment errors |
| `HighSOR` | `SOR > ...` | SOR = Strand Odds Ratio; high values mean strand bias |
| `LowReadPosRankSum` | `< ...` | Position bias in reads; low values indicate edge-biased calls |
| `LowBaseQRankSum` | `< ...` | Tests if base quality is lower in variant-supporting reads |

These thresholds are **fully configurable** via `params.gatk_vf_*` variables in the pipeline configuration file.

---

## Additional Options

| Option | Purpose |
|--------|---------|
| `--cluster-window-size` | Number of base pairs to define a cluster of nearby variants |
| `--cluster-size` | Max number of clustered variants allowed within the window |

> This helps identify artifact-prone regions where too many variants are called together.

---

## Output Validation

The script includes built-in checks to ensure the output VCF and its index are correctly generated:

```bash
if [ ! -s ${sample_id}_filtered.vcf.gz ] || [ ! -s ${sample_id}_filtered.vcf.gz.tbi ]; then
    echo "Error: Filtered VCF or index is empty for ${sample_id}" >&2
    exit 1
fi
```

---

## How to Inspect Results

### View Filtered Variants

```bash
bcftools view -f "PASS" sample_filtered.vcf.gz | less
```

### Count Total vs. Passed Variants

```bash
bcftools stats sample_filtered.vcf.gz | grep -E "number of SNPs|number of indels|number of records"
```

Check how many variants passed the filters (`PASS`) vs. were flagged.

---

## Summary

The `GATK_VARIANT_FILTER` step is a configurable, hard-filtering step that marks or removes low-quality variant calls. It ensures that the final VCF file contains **only high-confidence calls** suitable for annotation, biological interpretation, or clinical applications.

With clearly named filters and full control via pipeline parameters, this step forms the last major quality gate in your variant calling workflow.

---

# Chapter 10: Variant Annotation with SnpEff (`ANNOTATE_VARIANTS` Process)

## Overview

The `ANNOTATE_VARIANTS` process uses **[SnpEff](http://snpeff.sourceforge.net/)** to annotate genetic variants in the VCF file with biological and functional information. SnpEff predicts the **effect of each variant** on genes, transcripts, and proteins — e.g., whether it's a synonymous mutation, missense variant, or a frameshift.

This step is critical for making your variant calls interpretable for research or clinical applications.

---

## Inputs

```nextflow
tuple val(sample_id), path(vcf), path(tbi)
path(snpEffJar)
path(snpEffConfig)
path(snpEffDbDir)
val(genomedb)
```

| Input | Description |
|-------|-------------|
| `sample_id` | Sample name |
| `vcf`, `tbi` | Filtered VCF and its Tabix index |
| `snpEffJar` | Path to the SnpEff JAR executable |
| `snpEffConfig` | Path to the SnpEff configuration file |
| `snpEffDbDir` | Directory containing SnpEff genome databases |
| `genomedb` | Genome database name (e.g. `GRCh38.86`) used for annotation |

---

## Outputs

```nextflow
tuple val(sample_id), path("annotated_${sample_id}.vcf")
path("annotated_${sample_id}.summary.html")
```

| Output File | Description |
|-------------|-------------|
| `annotated_${sample_id}.vcf` | VCF file enriched with SnpEff annotations |
| `annotated_${sample_id}.summary.html` | HTML report summarizing annotation statistics |

All output files are copied to the directory defined by:

```nextflow
publishDir params.annotate_outdir, mode: 'copy'
```

---

## Core Annotation Command

```bash
java -Xmx16G -jar ${snpEffJar} \
    -c ${snpEffConfig} \
    -v ${genomedb} \
    -dataDir ${snpEffDbDir} \
    ${vcf} > annotated_${sample_id}.vcf
```

### Parameter Explanations

| Parameter | Purpose |
|----------|---------|
| `-Xmx16G` | Allocates 16GB of memory for SnpEff |
| `-jar ${snpEffJar}` | Executes the SnpEff JAR file |
| `-c ${snpEffConfig}` | Configuration file specifying genome info, cache path, and databases |
| `-v ${genomedb}` | Name of the SnpEff genome database (e.g., GRCh38.86) |
| `-dataDir ${snpEffDbDir}` | Location of SnpEff genome databases |
| `${vcf}` | Input VCF to annotate |

This command reads each variant in the VCF and appends annotation fields like:

- `ANN=...` with consequence, impact, gene name, biotype, transcript ID
- Effects ranked by predicted severity

---

## HTML Summary Report

A second SnpEff command generates a **visual report** showing:

```bash
java -Xmx16G -jar ${snpEffJar} \
    -c ${snpEffConfig} \
    -v ${genomedb} \
    -dataDir ${snpEffDbDir} \
    -stats annotated_${sample_id}.summary.html \
    ${vcf} > /dev/null
```

This outputs:

- Summary of all variant types (SNP, InDel, StopGained, etc.)
- Bar plots and tables of functional effects
- Total count of variants per class and impact

---

## How to Inspect Results

### Annotated VCF

```bash
zcat annotated_sample.vcf | grep -v "^#" | head
```

- Look for the `ANN=` field in the INFO column
- Confirm that variants are annotated with gene names and effects

### HTML Report

Open in a browser:

```
open annotated_sample.summary.html
```

Check for:

- Proportion of high/moderate/low impact variants
- Most frequently affected genes and regions
- Types of mutation consequences (e.g. missense, frameshift)

---

## Summary

The `ANNOTATE_VARIANTS` process enhances filtered variant calls with rich biological context using **SnpEff**. This enables researchers and clinicians to interpret the **functional impact** of each variant — a crucial step before downstream filtering, prioritization, or reporting.

It produces both an annotated VCF and a graphical summary report, making it easier to communicate and analyze the results.

---

# Chapter 11: Variant Annotation with Ensembl VEP (`ANNOTATEVARIANTS_VEP` Process)

## Overview

The `ANNOTATEVARIANTS_VEP` process uses the **[Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)** to annotate genetic variants with information about their predicted effects on genes, transcripts, proteins, and regulatory elements.

VEP is one of the most comprehensive variant annotation tools and is especially useful for clinical-grade annotations, providing HGNC symbols, transcript IDs, functional consequence terms, population allele frequencies, and links to known variant databases.

---

## Inputs

```nextflow
tuple val(sample_id), path(input_vcf), path(input_vcf_tbi)
path vep_cache    
val genome_assembly
val cache_version
val species
```

| Input | Description |
|-------|-------------|
| `sample_id` | Unique sample identifier |
| `input_vcf`, `input_vcf_tbi` | Input filtered VCF and Tabix index |
| `vep_cache` | Path to local Ensembl VEP cache directory |
| `genome_assembly` | Genome assembly name (e.g., GRCh38 or GRCh37) |
| `cache_version` | VEP cache version (e.g., 108) |
| `species` | Species identifier (e.g., `homo_sapiens`) |

---

## Outputs

```nextflow
tuple val(sample_id), path("vep_annotated_${sample_id}.vcf")
path("vep_annotated_${sample_id}.html")
path("versions.yml")
```

| Output File | Description |
|-------------|-------------|
| `vep_annotated_${sample_id}.vcf` | VCF annotated by VEP with INFO fields and transcript consequences |
| `vep_annotated_${sample_id}.html` | Summary statistics and charts in HTML |
| `versions.yml` | YAML file containing VEP version metadata |

Output files are written to the directory defined by:

```nextflow
publishDir params.annotate_vep_outdir, mode: 'copy'
```

---

## Core Command

```bash
vep \
  --input_file "${input_vcf}" \
  --output_file "vep_annotated_${sample_id}.vcf" \
  --stats_file "vep_annotated_${sample_id}.html" \
  --cache \
  --dir_cache "${vep_cache}" \
  --species "${species}" \
  --assembly "${genome_assembly}" \
  --cache_version ${cache_version} \
  --format vcf \
  --vcf \
  --symbol \
  --protein \
  --check_existing \
  --everything \
  --filter_common \
  --per_gene \
  --total_length \
  --force_overwrite \
  --offline
```

---

## Important VEP Flags and Their Purpose

| Flag | Description |
|------|-------------|
| `--cache` | Use local VEP cache for speed and offline use |
| `--offline` | Prevents online API calls; requires full local cache |
| `--dir_cache` | Location of pre-downloaded VEP cache |
| `--species` | Species identifier (e.g. `homo_sapiens`) |
| `--assembly` | Genome build used for annotation (e.g. GRCh38) |
| `--cache_version` | Specifies which cache version to use (e.g. 108) |
| `--symbol` | Adds HGNC gene symbols |
| `--protein` | Includes protein sequence information |
| `--check_existing` | Flags variants already known in public databases |
| `--everything` | Enables all available standard annotations |
| `--filter_common` | Filters out variants common in public populations |
| `--per_gene` | Reports only one consequence per gene |
| `--force_overwrite` | Overwrites any existing output files |
| `--stats_file` | Generates an interactive HTML summary of annotated variants |

---

## Version Capture

VEP version is logged for reproducibility using:

```bash
vep --help | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*$//'
```

Saved to:

```yaml
versions.yml
```

---

## Output Example

### Annotated VCF

Inspect it with:

```bash
zcat vep_annotated_sample.vcf | grep -v "^#" | head
```

Look for the `CSQ=` field in the INFO column with transcript and protein-level annotations.

### HTML Summary

Open in a browser:

```
open vep_annotated_sample.html
```

Includes:
- Charts on variant effects, consequence types
- Number of affected genes and transcripts
- Most common variant impacts

---

## When to Use VEP vs. SnpEff

| Tool | Best For |
|------|----------|
| **SnpEff** | Fast, simple annotation for large datasets |
| **VEP** | Clinical-grade, rich functional and population annotations |

Some pipelines use both to cross-reference annotations or for specific downstream needs.

---

## Summary

The `ANNOTATEVARIANTS_VEP` process uses **Ensembl VEP** to enrich your final VCF with detailed, biologically and clinically relevant annotations. It provides powerful transcript-level and population-level insights into each variant and generates both a structured VCF and an interactive HTML report.

This is the final interpretation layer in your variant calling pipeline and forms the foundation for filtering, prioritization, and reporting.

---

# Chapter 12: Fusion Gene Detection with Arriba (`ARRIBA` Process)

## Overview

The `ARRIBA` process uses the **[Arriba](https://github.com/suhrig/arriba)** tool to detect gene fusions from RNA-seq alignments. Gene fusions are critical in cancer biology, as they can create abnormal, oncogenic transcripts. This step analyzes STAR-aligned RNA-seq BAMs and identifies high-confidence fusion events.

---

## Inputs

```nextflow
tuple val(sample_id), path(bam), path(bai)
path fasta
path gtf
path blacklist
path known_fusions
```

| Input | Description |
|-------|-------------|
| `sample_id` | Sample identifier |
| `bam`, `bai` | BAM file and index from STAR alignment |
| `fasta` | Reference genome FASTA file |
| `gtf` | Gene annotation file (e.g., GENCODE, Ensembl) |
| `blacklist` | List of known false positives (provided by Arriba) |
| `known_fusions` | Database of validated known fusions (optional for prioritization) |

---

## Outputs

```nextflow
tuple val(sample_id), path("*.fusions.tsv"), emit: fusions
tuple val(sample_id), path("*.fusions.discarded.tsv"), emit: fusions_discarded
path("versions.yml"), emit: versions
```

| Output | Description |
|--------|-------------|
| `*.fusions.tsv` | High-confidence gene fusions detected |
| `*.fusions.discarded.tsv` | Low-confidence fusions that were filtered out |
| `versions.yml` | Arriba version information for reproducibility |

All outputs are copied to the path defined by:

```nextflow
publishDir params.arriba_outdir, mode: 'copy'
```

---

## Core Command

```bash
arriba \
  -x "${bam}" \
  -a "${fasta}" \
  -g "${gtf}" \
  -b "${blacklist}" \
  -k "${known_fusions}" \
  -o "${sample_id}.fusions.tsv" \
  -O "${sample_id}.fusions.discarded.tsv"
```

### Parameter Explanations

| Flag | Purpose |
|------|---------|
| `-x` | Input BAM file aligned with STAR |
| `-a` | Reference genome FASTA |
| `-g` | Gene annotation GTF |
| `-b` | Blacklist file to remove false positives (e.g., readthroughs, artifacts) |
| `-k` | Known fusion database to help prioritize clinically relevant fusions |
| `-o` | Output: detected fusion genes |
| `-O` | Output: filtered/discarded fusions |

---

## Version Logging

The version of Arriba is automatically captured using:

```bash
arriba -h | grep 'Version:'
```

Saved in:

```yaml
versions.yml
```

---

## Output Inspection

### Fusions File (`*.fusions.tsv`)

This file contains high-confidence gene fusions with fields like:

- `gene1`, `gene2`: Names of fused genes
- `breakpoint1`, `breakpoint2`: Chromosomal positions of fusion
- `confidence`: `high`, `medium`, `low`
- `type`: Fusion mechanism (e.g., duplication, inversion, read-through)
- `reads1`, `reads2`: Read support for fusion junctions

### Discarded Fusions (`*.fusions.discarded.tsv`)

- Contains fusions filtered out based on confidence, blacklist, or low read support
- Still useful for manual review or exploratory analysis

---

## Visualization

You can visualize fusion results using Arriba’s companion R scripts, IGV snapshots, or tools like **FusionInspector** or **FusionViewer**.

---

## Summary

The `ARRIBA` process performs fusion gene discovery using RNA-seq data and STAR-aligned BAMs. It distinguishes real fusion events from technical artifacts using a combination of:

- Reference genome and gene annotations
- Known fusion databases
- Blacklists of common false positives

This step is essential in cancer studies and disease transcriptomics, where fusion genes can be drivers of pathology or therapeutic targets.

---

# Chapter 13: Fusion Gene Visualization with Arriba (`ARRIBA_VISUALIZATION` Process)

## Overview

The `ARRIBA_VISUALIZATION` process uses the **Arriba `draw_fusions.R`** script to generate publication-ready **PDF plots of gene fusions** detected in RNA-seq data. These plots help visualize the structure of gene fusions, their breakpoints, protein domains, and associated transcripts.

This step supports both real fusion events and placeholder plots when no fusions are found.

---

## Inputs

```nextflow
tuple val(sample_id), path(bam), path(bai), path(fusions_tsv)
path gtf
```

| Input | Description |
|-------|-------------|
| `sample_id` | Sample name |
| `bam`, `bai` | STAR-aligned RNA-seq BAM and index file |
| `fusions_tsv` | Gene fusions identified by `arriba` (`*.fusions.tsv`) |
| `gtf` | Gene annotation GTF used for plotting transcript structures |

Optional parameters via `params`:

| Parameter | Description |
|----------|-------------|
| `params.cytobands` | Cytoband file (optional, for chromosomal band annotation) |
| `params.protein_domains` | Protein domain file (optional, for domain annotation in the plot) |

---

## Outputs

```nextflow
tuple val(sample_id), path("*.pdf"), emit: fusion_plot
path "versions.yml", emit: versions
```

| Output | Description |
|--------|-------------|
| `*.pdf` | A single PDF file visualizing detected fusions for the sample |
| `versions.yml` | Arriba version metadata for reproducibility |

---

## Core Workflow Logic

This process automatically handles both:
- **Samples with detected fusions** → generates a detailed plot
- **Samples without fusions** → generates a placeholder "no fusions found" PDF

---

## Fusion Plot Command

```bash
draw_fusions.R \
  --fusions=${fusions_tsv} \
  --alignments=${bam} \
  --output=${sample_id}.fusion_plot.pdf \
  --annotation=${gtf} \
  --cytobands=${params.cytobands} \
  --proteinDomains=${params.protein_domains}
```

### Optional Parameters

| Flag | Description |
|------|-------------|
| `--cytobands` | Adds chromosomal band information to fusion plots |
| `--proteinDomains` | Adds protein domain information (Pfam/SMART) |
| `--annotation` | Gene model used to draw exon/intron structures |
| `--fusions` | Input TSV from Arriba |
| `--alignments` | BAM file used to show read support for fusions |

---

## No Fusion Case Handling

If the `*.fusions.tsv` file is empty, a placeholder PDF is generated:

```bash
echo "No fusions found for ${sample_id}." \
  | convert -background white -fill black -font Helvetica -pointsize 20 text:- ${sample_id}.fusion_plot.pdf
```

This ensures that the workflow always emits a plot file, regardless of fusion detection outcome.

---

## Version Logging

The Arriba version is recorded using:

```bash
arriba -h | grep 'Version:'
```

Output is stored in:

```yaml
versions.yml
```

---

## How to View Fusion Plot

Open the PDF using a viewer:

```
open sample1.fusion_plot.pdf
```

Each page of the PDF includes:

- Fused genes and breakpoints
- Exon structures from the GTF
- Optional protein domains
- Read-level support from the alignment

---

## Summary

The `ARRIBA_VISUALIZATION` process produces clear, visual representations of fusion events detected in RNA-seq data. It enhances the interpretability of fusion calls by showing the structural impact of each fusion on gene models and protein domains. Even when no fusions are found, a placeholder is emitted to ensure consistent outputs for downstream reporting.

---
