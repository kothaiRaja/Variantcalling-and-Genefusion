

# Chapter 1: Preprocessing of Raw Reads (`PREPROCESSING` Subworkflow)

## Overview

The `PREPROCESSING` subworkflow performs **read quality control and trimming** before alignment. It prepares high-quality input for downstream analysis by:

- Optionally concatenating multiple FASTQ files per sample
- Running quality control on raw reads using **FastQC**
- Trimming adapters and low-quality bases using **Fastp**
- Aggregating all QC results using **MultiQC**
- Dumping software versions used for reproducibility

---

## Inputs

```nextflow
take:
    samplesheet
    dump_script
```

| Input | Description |
|-------|-------------|
| `samplesheet` | CSV file with columns: `sample_id`, `fastq_1`, `fastq_2`, `strandedness` |
| `dump_script` | Script required by `nf-core/software_versions` to collect tool versions |

---

## Workflow Logic

### 1. **Read the Sample Sheet**

Parses the samplesheet using `.splitCsv()` and maps each row to a tuple:

```groovy
tuple(sample_id, [file(fastq_1), file(fastq_2)], strandedness)
```

### 2. **Concatenate FASTQ Files (Optional)**

If `params.concatenate = true`, the workflow calls:

```nextflow
CONCAT_FASTQ(samples_ch)
```

Otherwise, it directly uses the input FASTQ files for each sample.

---

### 3. **Run FastQC on Raw Reads**

FastQC is applied on the untrimmed FASTQ files to assess:

- Per base sequence quality
- Adapter content
- GC content
- Sequence duplication levels

The results are mixed into the report and version tracking channels.

---

### 4. **Trim Reads with Fastp**

Fastp performs:

- Adapter detection and trimming
- Sliding window base quality filtering
- Per-read quality assessment

Outputs include:

- Trimmed `R1/R2` FASTQs
- Fastp HTML and JSON reports

Reports are logged and added to MultiQC input.

---

### 5. **Dump Software Versions**

Uses the `CUSTOM_DUMPSOFTWAREVERSIONS` module to:

- Collect versions of all tools used
- Generate a `.yml` file for MultiQC

This enables full reproducibility tracking.

---

### 6. **Generate Combined QC Report with MultiQC**

The reports from FastQC, Fastp, and version logs are:

- Collected using `.collect()`
- Flattened and passed to `MultiQC`
- A final interactive HTML report is generated summarizing all samples

---

## Emitted Channels

```nextflow
emit:
    qc_results     = qc_results_ch
    fastp_reports  = fastp_reports_ch
    trimmed_reads  = trimmed_reads_ch
    reports        = reports_ch
    multiqc        = multiqc_quality.report
    versions       = ch_versions
```

| Output | Description |
|--------|-------------|
| `qc_results` | FastQC results for raw reads |
| `fastp_reports` | HTML and JSON reports from Fastp |
| `trimmed_reads` | Paired-end trimmed FASTQs |
| `reports` | All collected reports (FastQC, Fastp, versions) |
| `multiqc` | Final HTML report aggregating all QC metrics |
| `versions` | Combined version metadata for all tools used |

---

## Key Modules Used

| Module | Purpose |
|--------|---------|
| `CONCAT_FASTQ` | Concatenate FASTQ files from multiple lanes or replicates |
| `FASTQC_RAW` | Run FastQC on raw reads |
| `TRIM_READS` | Trim reads using Fastp |
| `CUSTOM_DUMPSOFTWAREVERSIONS` | Dump software versions for traceability |
| `MultiQC_quality` | Generate summary report using MultiQC |

---

## Highlights

-  Handles both raw FASTQ and lane-merged samples
-  Comprehensive QC for both pre- and post-trimming reads
-  MultiQC summary is easy to visualize and export
-  Debug-friendly `.view` calls let users trace what’s passed to each step

---

## Summary

The `PREPROCESSING` subworkflow ensures your raw sequencing data is clean, consistent, and ready for alignment. By trimming, checking quality, and generating comprehensive reports, this step lays a solid foundation for trustworthy downstream analysis in any RNA-seq or variant-calling pipeline.

---

# Chapter 2:Subworkflow: STAR_ALIGN

## Purpose

The `STAR_ALIGN` subworkflow performs **RNA-seq read alignment and preprocessing** using the STAR aligner and several post-alignment tools from SAMtools. It is designed to:

- Align paired-end reads using STAR (if pre-aligned BAMs are not provided)
- Sort and index aligned BAMs
- Generate alignment statistics and flagstat reports
- Filter orphaned reads
- Track versions for reproducibility

---

## Inputs

| Input | Type | Description |
|-------|------|-------------|
| `trimmed_reads_ch` | `tuple(val(sample_id), path(r1), path(r2), val(strandedness))` | Trimmed paired-end FASTQ reads with sample ID and strand info |
| `star_index` | `path` | Prebuilt STAR genome index directory |
| `gtf_file` | `path` | Gene annotation file (GTF) used by STAR |
| `aligned_bam_samplesheet` | `path` (optional) | CSV with sample ID, strandedness, and BAM path |
| `aligned_bam_folder` | `path` (optional) | Folder with pre-aligned BAMs |

---

## Logic Flow

###  Conditional Branching

- If **pre-aligned BAMs are *not* provided**, it:
  - Runs the `STAR_ALIGNMENT` module
  - Extracts BAMs, chimeric reads, STAR logs, and version info

- If **pre-aligned BAMs *are* provided**, it:
  - Reads from a CSV (`aligned_bam_samplesheet`) or folder (`aligned_bam_folder`)
  - Assumes strandedness as `"unstranded"` if missing
  - Skips STAR alignment and emits empty STAR-specific outputs

---

## Modules Included

| Module | Purpose |
|--------|---------|
| `STAR_ALIGNMENT` | Align reads with STAR using 2-pass alignment |
| `SAMTOOLS_SORT_INDEX` | Sort and index aligned BAMs |
| `SAMTOOLS_STATS` | Generate alignment summary stats |
| `SAMTOOLS_FILTER_ORPHANS` | Remove orphaned reads |
| `SAMTOOLS_FLAGSTAT` | Run flagstat for alignment quality metrics |

---

## Intermediate Channels

| Channel | Description |
|---------|-------------|
| `star_bam_ch` | BAMs from STAR or pre-aligned inputs |
| `chimeric_reads_ch` | STAR chimeric read output for fusion detection |
| `flagstats_ch` | BAM-level alignment quality metrics |
| `align_stats_ch` | Alignment summary from SAMtools stats |
| `filtered_bams_ch` | Orphan-filtered BAMs for downstream analysis |
| `ch_versions` | Version metadata for STAR and SAMtools modules |

---

## Emitted Outputs

```nextflow
emit:
    bam_sorted     = sorted_bams_ch
    chimeric_reads = chimeric_reads_ch
    flagstats      = flagstats_ch
    align_stats    = align_stats_ch
    star_logs      = star_logs_ch
    filtered_bams  = filtered_bams_ch
    versions       = ch_versions
```

| Output | Type | Description |
|--------|------|-------------|
| `bam_sorted` | `tuple(...)` | Sorted and indexed BAMs from STAR or pre-aligned sources |
| `chimeric_reads` | `tuple(...)` | Chimeric reads from STAR (used in fusion calling) |
| `flagstats` | `tuple(...)` | Flagstat output for filtered BAMs |
| `align_stats` | `tuple(...)` | Alignment statistics (e.g., % mapped) |
| `star_logs` | `tuple(...)` | STAR alignment logs (`Log.final.out`) |
| `filtered_bams` | `tuple(...)` | BAMs after orphan filtering |
| `versions` | `path("versions.yml")` | Combined version info from all modules |

---

## Highlights

-  **Flexibility**: Works with raw FASTQs *or* pre-aligned BAMs
-  **Multi-level QC**: Combines STAR logs, SAMtools stats, and flagstat outputs
-  **Fusion-ready**: Emits chimeric reads from STAR for tools like Arriba
-  **Reproducible**: Collects and emits versions of every module used

---

## Example Use

```nextflow
workflow {
    call STAR_ALIGN {
        take:
        trimmed_reads_ch = my_trimmed_fastqs
        star_index = params.star_index
        gtf_file = file(params.gtf)
    }
}
```

---

## Tips

- Use `aligned_bam_samplesheet` if you already have STAR-aligned BAMs
- Make sure BAM file names follow STAR naming conventions if using `aligned_bam_folder`
- Use emitted `filtered_bams` for downstream variant calling or quantification
- Combine `chimeric_reads` with a fusion detection module like `ARRIBA`

---

## Summary

The `STAR_ALIGN` subworkflow handles RNA-seq alignment in a robust and modular way, supporting both fresh alignment with STAR and reuse of existing BAMs. It combines best-practice preprocessing and multiple levels of QC for optimal downstream performance.

---

# Chapter 3: Interval Processing (`INTERVAL_PROCESSING` Subworkflow)

## Overview

The `INTERVAL_PROCESSING` subworkflow handles the preparation of genomic intervals for downstream tools like **BaseRecalibrator** and **HaplotypeCaller**. It supports both:

- Converting a **BED file** to a GATK-compatible `.interval_list`
- **Scattering** the interval list into smaller chunks for parallel execution

This step ensures compatibility and performance optimization in workflows that rely on genomic intervals.

---

## Inputs

```nextflow
take:
    bed_file_ch
    reference_genome
    reference_genome_dict
```

| Input | Type | Description |
|-------|------|-------------|
| `bed_file_ch` | `tuple(val(meta_id), path(bed))` | BED file with regions of interest (e.g., exome targets) |
| `reference_genome` | `path` | Reference genome FASTA file |
| `reference_genome_dict` | `path` | Reference genome dictionary (`.dict`), required by GATK tools |

---

## Main Workflow Logic

###  Step 1: Convert BED → Interval List

```nextflow
interval_list = BED_TO_INTERVAL_LIST(bed_file_ch, reference_genome, reference_genome_dict)
```

- Uses the GATK `BedToIntervalList` tool
- Ensures that intervals are compatible with downstream GATK processes
- Emits `.interval_list` for the given BED

Version is recorded in:

```nextflow
BED_TO_INTERVAL_LIST.out.versions
```

---

###  Step 2: Scatter Intervals (Optional)

Controlled via:

```nextflow
params.scatterintervals = true
```

If enabled:

```nextflow
scattered_intervals = SCATTER_INTERVAL_LIST(interval_list_ch, params.reference_genome_dict)
```

- Uses `IntervalListTools` from GATK to split the interval list
- Enables scatter/gather parallelization in downstream processes (e.g., BQSR)
- Each interval chunk is flattened into a separate file path

Each scattered `.interval_list` file is displayed with `.view {}` for debugging:

```
Scattered interval: chr1_0_1000000.interval_list
```

If not enabled, the full unscattered list is passed as-is.

---

## Outputs

```nextflow
emit:
    intervals = scattered_intervals_ch
    versions  = ch_versions
```

| Output | Description |
|--------|-------------|
| `intervals` | Scattered or full `.interval_list` files (depending on config) |
| `versions` | Combined version info from the interval tools used |

---

## Modules Used

| Module | Purpose |
|--------|---------|
| `BED_TO_INTERVAL_LIST` | Converts BED → `.interval_list` (GATK-compatible) |
| `SCATTER_INTERVAL_LIST` | Splits interval list for parallel processing (optional) |

---

## Highlights

-  Compatible with **GATK Best Practices**
-  Enables **parallelization** through interval scattering
-  Can be turned off to use a **single interval list**
-  Emission of versions supports reproducibility

---

## Example Use Case

If you're calling variants only on exonic regions:

```nextflow
workflow {
    call INTERVAL_PROCESSING {
        take:
        bed_file_ch = Channel.fromPath("targets.bed").map { tuple("target", it) }
        reference_genome = file(params.genome_fasta)
        reference_genome_dict = file(params.genome_dict)
    }
}
```

---

## Summary

The `INTERVAL_PROCESSING` subworkflow ensures that your BED file is transformed into a GATK-friendly format and can optionally be scattered for performance. This makes it essential for any pipeline involving **BQSR**, **HaplotypeCaller**, or other GATK tools that operate on intervals.

---

# Chapter 4: Duplicate Marking (`MARK_DUPLICATES` Subworkflow)

## Overview

The `MARK_DUPLICATES` subworkflow identifies and marks **PCR duplicates** in aligned BAM files using **GATK's `MarkDuplicates`** tool. This step is essential to avoid overcounting reads that arise from amplification artifacts during sequencing, particularly important for accurate variant calling.

---

## Inputs

```nextflow
take:
    bam_input_ch
```

| Input | Type | Description |
|-------|------|-------------|
| `bam_input_ch` | `tuple(val(sample_id), val(strandedness), path(bam), path(bai))` | Sorted BAM and index from alignment step or STAR subworkflow |

---

## Main Workflow Logic

```nextflow
dedup_bams = GATK_MARK_DUPLICATES(bam_input_ch)
```

- Calls the `GATK_MARK_DUPLICATES` module
- Uses GATK’s Picard-style MarkDuplicates to:
  - Detect and mark duplicate reads
  - Create indexed BAMs
  - Generate duplication metrics

---

## Outputs

```nextflow
emit:
    marked_bams_bai        = dedup_bams_ch
    marked_bams_bai_metrics = dedup_metrics_ch
    versions               = ch_versions
```

| Output | Type | Description |
|--------|------|-------------|
| `marked_bams_bai` | Tuple with BAM and BAI after duplicates are marked |
| `marked_bams_bai_metrics` | Duplicate metrics file (`.dup_metrics.txt`) |
| `versions` | Collected tool version info from `GATK_MARK_DUPLICATES` |

The `.collect { it[2] }` operation extracts just the metrics files from the 3-element tuple (sample_id, strandedness, metrics).

---

## Tool: GATK `MarkDuplicates`

### Key Parameters (from module)

| Parameter | Description |
|-----------|-------------|
| `-I` | Input BAM |
| `-O` | Output BAM with marked duplicates |
| `-M` | Duplication metrics file |
| `--CREATE_INDEX true` | Ensures a BAM index is also created |
| `--REMOVE_DUPLICATES` | Optional: remove instead of marking (via `params`) |
| `--VALIDATION_STRINGENCY LENIENT` | Avoids hard errors for minor formatting issues |

---

## Highlights

-  Detects PCR duplicates using coordinate and read group
-  Emits detailed duplication metrics
-  Fully parallelizable per sample
-  Version tracking included

---

## Example Use

If you're running this step standalone:

```nextflow
workflow {
    call MARK_DUPLICATES {
        take:
        bam_input_ch = some_sorted_bams
    }
}
```

---

## Summary

The `MARK_DUPLICATES` subworkflow wraps GATK’s `MarkDuplicates` process into a scalable, version-aware unit that ensures high-quality input for recalibration and variant calling. By identifying PCR artifacts, it enhances the **accuracy and reliability** of all downstream analysis.

---

# Chapter 5: Split and Merge BAMs (`SPLIT_MERGE_BAMS` Subworkflow)

## Overview

The `SPLIT_MERGE_BAMS` subworkflow processes RNA-seq aligned BAM files to make them compatible with **GATK's variant calling** tools. It performs:

1. **SplitNCigarReads** — for handling spliced reads
2. **Merging BAMs** — across genomic intervals
3. **Applying `calmd`** — for adding MD tags required by some tools

This ensures BAMs are cleaned, standardized, and properly annotated across scattered genomic intervals.

---

## Inputs

```nextflow
take:
    bam_input_ch           // BAMs with duplicates marked
    intervals_ch           // List of genomic intervals (scattered or full)
    reference_genome
    reference_genome_index
    reference_genome_dict
```

| Input | Description |
|-------|-------------|
| `bam_input_ch` | Tuples: `(sample_id, strandedness, bam, bai)` |
| `intervals_ch` | `.interval_list` files, possibly scattered |
| `reference_genome` | Reference FASTA file |
| `reference_genome_index` | `.fai` index for the FASTA |
| `reference_genome_dict` | Sequence dictionary file (`.dict`) |

---

## Main Workflow Steps

###  Step 1: Split BAMs using `SplitNCigarReads`

RNA-seq aligners like STAR use `N` CIGAR operations for spliced alignments. GATK tools like `HaplotypeCaller` require **unsplit reads**, so this step:

```nextflow
split_bams = SPLIT_NCIGAR_READS(ch_splitncigar_bam_bai_interval, reference_genome, reference_genome_index, reference_genome_dict)
```

- Runs once per `(sample, interval)` tuple
- Uses interval scattering to parallelize
- Emits split BAMs + index files

Output: `split_bams_ch`

---

###  Step 2: Merge Split BAMs Per Sample

All split BAMs for each sample are grouped and merged using:

```nextflow
merged_bams = MERGE_BAMS(ch_merged_bams)
```

This creates a single BAM per sample across all processed intervals.

Output: `merged_bams_ch`

---

###  Step 3: Apply `samtools calmd`

```nextflow
calmd_bams = SAMTOOLS_CALMD(merged_bams_ch, reference_genome, reference_genome_index)
```

- Adds **MD tags** and **NM fields** to BAM alignments
- These are essential for accurate base recalibration and variant calling

Output: `calmd_bams_ch`

---

## Key Transformations

- `.combine(intervals_ch)` pairs each BAM with every interval
- `.groupTuple()` groups results by sample for merging
- `.flatten()` ensures multi-interval files are joined into unified lists

These steps enable scalable, interval-based processing followed by per-sample aggregation.

---

## Outputs

```nextflow
emit:
    split_bams         = split_bams_ch
    merged_bams        = merged_bams_ch
    merged_calmd_bams  = calmd_bams_ch
    versions           = ch_versions
```

| Output | Description |
|--------|-------------|
| `split_bams` | Split BAMs by interval, suitable for BQSR and variant calling |
| `merged_bams` | Sample-level merged BAMs from scattered splits |
| `merged_calmd_bams` | Merged BAMs with MD/NM tags added |
| `versions` | Version info from all modules (SplitNCigar, Merge, Calmd) |

---

## Modules Used

| Module | Purpose |
|--------|---------|
| `SPLIT_NCIGAR_READS` | Fix spliced reads for GATK compatibility |
| `MERGE_BAMS` | Merge split BAMs back into per-sample BAMs |
| `SAMTOOLS_CALMD` | Annotate BAMs with base-level MD/NM tags |

---

## Summary

The `SPLIT_MERGE_BAMS` subworkflow enables parallelized RNA-seq BAM processing over intervals and ensures compatibility with GATK variant calling tools. With interval-aware splitting, merging, and tag restoration, it forms a critical bridge between alignment and recalibration/calling.

---

# Chapter 6: Base Quality Recalibration (`BASE_RECALIBRATION` Subworkflow)

## Overview

The `BASE_RECALIBRATION` subworkflow corrects systematic errors in base quality scores using known variant sites. It wraps two key GATK tools:

1. `BaseRecalibrator` – Builds a model of covariates based on known variants.
2. `ApplyBQSR` – Applies the learned model to adjust quality scores in BAM files.

This process improves the accuracy of variant calling by modeling technical biases in the sequencer’s quality estimations.

---

## Inputs

```nextflow
take:
    bam_input_ch           // BAMs post-calmd/merge
    intervals_ch           // Scattered interval list
    reference_genome
    reference_genome_index
    reference_genome_dict
    merged_vcf             // Known variants (e.g., dbSNP)
    merged_vcf_index       // Index for known variants
```

| Input | Description |
|-------|-------------|
| `bam_input_ch` | BAMs from prior step (tuple: sample_id, strandedness, BAM, BAI) |
| `intervals_ch` | Scattered `.interval_list` files |
| `reference_genome`, `index`, `dict` | FASTA + support files for reference |
| `merged_vcf`, `merged_vcf_index` | Known sites of variation (e.g., dbSNP, Mills) |

---

## Workflow Logic

###  Step 1: Base Recalibrator (BQSR table generation)

Each input BAM is paired with each interval:

```nextflow
bam_input_ch.combine(intervals_ch)
```

Then run through `GATK_BASERECALIBRATOR`:

```nextflow
GATK_BASERECALIBRATOR(
    sample_id, strandedness, bam, bai, interval, 
    reference_genome, index, dict, known_variants
)
```

**Output**: A recalibration table (`*.recal_data.table`) per sample and interval.

---

###  Step 2: Join BAMs and Recal Tables

The BAMs and recal tables are joined by `sample_id`, producing tuples like:

```groovy
(sample_id, strandedness, bam, bai, recal_table)
```

This ensures each recalibrated table is correctly paired with the originating BAM.

---

###  Step 3: Apply Recalibration

The combined BAM + recal_table tuples are again paired with intervals:

```groovy
(sample_id, strandedness, bam, bai, recal_table, interval)
```

Then passed to `GATK_APPLYBQSR`, which outputs:

```bash
sampleID_chr1_0_1000000_recalibrated.bam
```

Each recalibrated BAM is still per-interval at this stage.

---

## Outputs

```nextflow
emit:
    recalibrated_bams = bams_base_recalibrated_ch
    versions = ch_versions
```

| Output | Description |
|--------|-------------|
| `recalibrated_bams` | Interval-level BAMs with adjusted quality scores |
| `versions` | Version metadata for BQSR tools used |

---

## Modules Used

| Module | Tool | Purpose |
|--------|------|---------|
| `GATK_BASERECALIBRATOR` | `BaseRecalibrator` | Builds the recalibration model from known variant sites |
| `GATK_APPLYBQSR` | `ApplyBQSR` | Applies the recalibration model to adjust base quality scores |

---

## Highlights

-  Ensures GATK best practices for RNA-seq data
-  Uses known variants to correct sequencing bias
-  Interval-parallelizable for performance
-  Clean channel joins using `.combine()` and `.join()`
-  Emits both results and versions for reproducibility

---

## Example Join Strategy (for ApplyBQSR)

Joining recalibration table with original BAMs:

```groovy
.join(recalibrated_bams_table_ch.map { sample_id, _, table -> tuple(sample_id, table) }, by: 0)
```

Ensures tables align with the correct BAMs by `sample_id`.

---

## Summary

The `BASE_RECALIBRATION` subworkflow performs base quality recalibration in two modular steps: modeling and application. With clean channel logic and interval-aware parallelism, it improves the reliability of variant calls by correcting for systemic biases in sequencing quality estimation.

---

# Chapter 7: Variant Calling (`VARIANT_CALLING` Subworkflow)

## Overview

The `VARIANT_CALLING` subworkflow performs **variant discovery, merging, filtering, and summary reporting** using a combination of:

- **GATK HaplotypeCaller**
- **GATK VariantFiltration**
- **BCFtools (stats and query)**
- **GATK SelectVariants (SNPs and INDELs)**

This step transforms interval-wise recalibrated BAMs into final, filtered VCFs ready for annotation or interpretation.

---

## Inputs

```nextflow
take:
    recalibrated_bams_ch
    reference_genome
    reference_genome_index
    reference_genome_dict
    known_variants
    known_variants_index
```

| Input | Description |
|-------|-------------|
| `recalibrated_bams_ch` | Recalibrated BAMs across intervals |
| `reference_*` | Reference FASTA and its index and dictionary |
| `known_variants` | Known VCF (e.g., dbSNP, Mills) for annotation |
| `known_variants_index` | Tabix index for the known VCF |

---

## Main Workflow Steps

###  Step 1: Variant Calling (GATK HaplotypeCaller)

```nextflow
GATK_HAPLOTYPE_CALLER(recalibrated_bams_ch, ...)
```

- Performs RNA-seq specific variant calling per interval
- Outputs VCFs and `.tbi` index files

These are grouped per sample using `.groupTuple()` and flattened:

```groovy
tuple(sample_id, strandedness, [vcfs...], [tbis...])
```

---

###  Step 2: Merge VCFs (GATK MergeVcfs)

```nextflow
GATK_MERGEVCFS(sample_id, vcf_list, tbi_list)
```

- Combines all interval-level VCFs into a single per-sample VCF
- Emits `.vcf.gz` and `.vcf.gz.tbi`

---

###  Step 3: Variant Statistics (BCFtools stats)

```nextflow
BCFTOOLS_STATS(merged_vcf_ch)
```

- Generates summary statistics such as:
  - SNP/INDEL counts
  - Transition/Transversion ratios
  - Quality distribution

---

###  Step 4: Variant Filtering (GATK VariantFiltration)

```nextflow
GATK_VARIANT_FILTER(merged_vcf_ch, ...)
```

- Applies hard filtering based on thresholds for:
  - QD, FS, MQ, SOR, QUAL, etc.
- Outputs final filtered `.vcf.gz` and index

---

###  Step 5: BCFtools Query (optional summary)

```nextflow
BCFTOOLS_QUERY(filtered_vcf)
```

- Extracts tabular representation (e.g., CHROM, POS, REF, ALT, INFO) for quick inspection or report generation

---

###  Step 6 (Optional): Select SNPs and INDELs Separately

If `params.select_variants = true`, the workflow:

```nextflow
SELECT_SNPs(filtered_vcf)
SELECT_INDELs(filtered_vcf)
```

- Separates SNPs and INDELs for downstream analysis
- Joins them for annotation in a combined channel

---

## Outputs

```nextflow
emit:
    final_variants    = ch_variant_filter
    selected_snps     = selected_snps_ch
    selected_indels   = selected_indels_ch
    selected_variants = ch_selected_variants
    bcftools_stats    = bcftools_stats_ch
    bcftools_query    = bcftools_query_ch
    versions          = ch_versions
```

| Output | Description |
|--------|-------------|
| `final_variants` | Fully filtered `.vcf.gz` file |
| `selected_snps` | Optional SNP-only VCF |
| `selected_indels` | Optional INDEL-only VCF |
| `selected_variants` | Paired SNPs + INDELs (for combined annotation) |
| `bcftools_stats` | Summary statistics from `bcftools stats` |
| `bcftools_query` | Tabular variant summary |
| `versions` | Version metadata for reproducibility |

---

## Modules Used

| Module | Description |
|--------|-------------|
| `GATK_HAPLOTYPE_CALLER` | Variant calling per BAM/interval |
| `GATK_MERGEVCFS` | Merges interval VCFs |
| `GATK_VARIANT_FILTER` | Applies hard filters |
| `BCFTOOLS_STATS` | Variant statistics |
| `BCFTOOLS_QUERY` | Extracts INFO fields for summary |
| `SELECT_SNPs`, `SELECT_INDELs` | Optional separation of variant types |

---

## Highlights

-  RNA-seq friendly calling via `SplitNCigarReads`-compatible BAMs
-  Variant filtering and decomposition using GATK + BCFtools
-  Tabular querying and stats for QA and reporting
-  Modular SNP/INDEL output for flexible downstream use

---

## Summary

The `VARIANT_CALLING` subworkflow is the analytical core of your RNA-seq pipeline. It performs per-interval variant calling, merging, filtering, and optional stratification into SNPs and INDELs — all designed for high sensitivity and specificity in transcriptomic variant discovery.

---

# Chapter 8: Variant Annotation (`ANNOTATE` Subworkflow)

## Overview

The `ANNOTATE` subworkflow enriches the called variants (VCFs) with biological, functional, and clinical annotations using the following tools:

- **SnpEff**: Local annotation tool that predicts effects based on transcript models.
- **VEP**: Ensembl's Variant Effect Predictor, offering rich annotations with regulatory, frequency, and impact data.
- **Combined Mode**: SnpEff followed by VEP to enhance completeness.

This subworkflow is **modular and conditional** — you can run any combination of tools by setting `params.tools` to include `'snpeff'`, `'vep'`, or `'combine'`.

---

## Inputs

```nextflow
take:
    vcf             // Tuple: [ val(meta_id), path(vcf) ]
    tools           // List: ["snpeff", "vep", or "combine"]
    snpeff_jar
    snpeff_config
    snpeff_db
    genomedb
    genome_assembly
    species
    vep_version
    vep_cache
```

| Input | Description |
|-------|-------------|
| `vcf` | Channel of called variants (VCF + meta) |
| `tools` | List of annotation tools to run (`["snpeff"]`, `["vep"]`, `["combine"]`) |
| `snpeff_jar`, `snpeff_config`, `snpeff_db`, `genomedb` | Inputs for SnpEff |
| `vep_cache`, `vep_version`, `genome_assembly`, `species` | Inputs for VEP |

---

## Workflow Logic

###  Case 1: SnpEff Annotation

```nextflow
if (tools.contains('snpeff') || tools.contains('combine')) {
    VARIANT_ANNOTATION(vcf, snpeff_jar, snpeff_config, snpeff_db, genomedb)
}
```

- Annotates variants using the SnpEff database for the specified genome
- Outputs include:
  - Annotated VCF (`annotated_sample.vcf`)
  - HTML summary (`annotated_sample.summary.html`)
- Emitted to:
  - `final_vcf_annotated`
  - `reports_html`

---

###  Case 2: Combined (SnpEff + VEP)

```nextflow
if (tools.contains('combine')) {
    COMBINED_ANNOTATE(...)
}
```

- Takes SnpEff-annotated VCF and runs Ensembl VEP
- Enhances annotations with:
  - Regulatory regions
  - Frequency in populations (gnomAD, 1000 Genomes)
  - Canonical transcript effects
- Combines outputs from both tools for rich annotation

---

###  Case 3: VEP Only

```nextflow
if (tools.contains('vep')) {
    VEP_ANNOTATION_WORKFLOW(...)
}
```

- Directly runs VEP on the filtered VCF
- Produces annotated VCF and summary report

---

## Outputs

```nextflow
emit:
    final_vcf_annotated     = annotated_vcfs
    reports_html            = annotation_reports
    versions                = ch_versions
```

| Output | Description |
|--------|-------------|
| `final_vcf_annotated` | Final annotated VCF(s), merged from all tools selected |
| `reports_html` | HTML summary reports from SnpEff and/or VEP |
| `versions` | Tool version metadata for reproducibility |

---

## Tools Used

| Tool | Description |
|------|-------------|
| `SnpEff` | Local transcript-aware effect prediction |
| `VEP` | Ensembl variant effect predictor for rich annotations |
| `Combined` | Chains SnpEff + VEP for comprehensive annotation coverage |

---

## Highlights

-  **Flexible configuration** via `params.tools`
-  Supports both **transcript** and **regulatory** annotation
-  **Multi-format outputs** (VCF + HTML)
-  Enables downstream variant filtering, prioritization, and reporting

---

## Example Configuration

In your `nextflow.config`:

```groovy
params.tools = ['combine']  // Run both snpEff and VEP
```

Or:

```groovy
params.tools = ['vep']      // VEP only
```

---

## Summary

The `ANNOTATE` subworkflow is a powerful, flexible, and modular final stage in your variant pipeline. It enhances variant calls with functional, regulatory, and clinical context — enabling biologically meaningful interpretation and downstream prioritization.

---

# Chapter 9: Gene Fusion Detection (`GENE_FUSION` Subworkflow)

## Overview

The `GENE_FUSION` subworkflow detects and visualizes **gene fusions** in RNA-seq data using the Arriba tool. Gene fusions can represent critical **oncogenic drivers** or **disease markers** in cancer and other pathologies.

This workflow consists of:
1. **Fusion Detection** using the `ARRIBA` process
2. **Fusion Visualization** using `ARRIBA_VISUALIZATION`

---

## Inputs

```nextflow
take:
    STAR_bam_output
    reference_genome
    gtf_annotation
    arriba_blacklist
    arriba_known_fusions
```

| Input | Type | Description |
|-------|------|-------------|
| `STAR_bam_output` | Tuple: `(sample_id, strandedness, bam, bai)` |
| `reference_genome` | FASTA file used for alignment |
| `gtf_annotation` | Gene model (GTF) used for fusion interpretation |
| `arriba_blacklist` | List of known false positives (Arriba-provided) |
| `arriba_known_fusions` | Known fusion reference database (e.g., Mitelman DB) |

---

## Workflow Logic

###  Step 1: Prepare ARRIBA Input

The input channel is simplified by dropping the `strandedness` field:

```nextflow
arriba_input_bam_ch = STAR_bam_output.map { sample_id, _, bam, bai -> tuple(sample_id, bam, bai) }
```

---

###  Step 2: Detect Fusions with ARRIBA

```nextflow
ARRIBA(arriba_input_bam_ch, reference_genome, gtf_annotation, arriba_blacklist, arriba_known_fusions)
```

- Alignments are scanned for fusion events
- Known artifacts are filtered using the blacklist
- Known fusions are labeled from a reference set

**Output**:
- `*.fusions.tsv`: High-confidence fusion calls
- `*.fusions.discarded.tsv`: Filtered/low-confidence calls

---

###  Step 3: Fusion Visualization

Input BAMs are joined to fusions using `sample_id`:

```nextflow
fusion_viz_input_ch = arriba_input_bam_ch.join(ARRIBA.out.fusions, by: 0)
```

This allows **visual context** to be added via:

```nextflow
ARRIBA_VISUALIZATION(fusion_viz_input_ch, gtf_annotation)
```

**Output**:
- `.pdf` files with gene diagrams, breakpoints, and read-level support

---

## Outputs

```nextflow
emit:
    fusion_results        = ARRIBA.out.fusions
    discarded_results     = ARRIBA.out.fusions_discarded
    fusion_visualizations = fusion_visual_ch
    versions              = ch_versions
```

| Output | Description |
|--------|-------------|
| `fusion_results` | Detected gene fusions (`.fusions.tsv`) |
| `discarded_results` | Low-confidence or filtered-out fusions |
| `fusion_visualizations` | PDF visual reports for each fusion |
| `versions` | Combined version info from ARRIBA and visualization modules |

---

## Tools Used

| Tool | Purpose |
|------|---------|
| `ARRIBA` | Detects gene fusions from STAR-aligned RNA-seq BAMs |
| `ARRIBA_VISUALIZATION` | Creates detailed PDF plots of detected fusions |

---

## Highlights

-  Sensitive detection of fusion events with STAR + ARRIBA
-  Optional visualization of fusion architecture and breakpoints
-  Filters known artifacts while preserving novel candidates
-  Emits both raw and visual results for downstream interpretation

---

## Example Visualization

PDF output includes:

- Fused genes and exons
- Genomic breakpoint positions
- Supporting read evidence
- Optional: protein domains (if configured)

---

