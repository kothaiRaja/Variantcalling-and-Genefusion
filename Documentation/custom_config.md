
# Appendix A: Understanding the `custom.config` File

## Purpose

This `custom.config` file is used to **centralize user-defined parameters** for the Nextflow pipeline. It allows you to:

- Customize input/output paths
- Configure tool settings
- Set computational resource allocations
- Enable/disable specific pipeline modules

By separating these parameters from the main workflow, it makes the pipeline:
- Easier to reuse
- More portable
- More user-friendly

---

## 1. `params {}` Block

### Base Directories

```groovy
outdir = "$baseDir/output"
resultsdir = "$baseDir/results"
base_dir = "$baseDir"
```

Standardized output paths for all results and files.

---

### Dataset Inputs

```groovy
samplesheet = "path/to/samples.csv"
reads = "path/to/sample_*_Aligned.sortedByCoord.out.bam"
strandedness_map = [
    "sample_1": "forward",
    "sample_2": "reverse",
    "sample_3": "unstranded"
]
```

Defines your sample sheet, raw read locations, and the strandness of each sample.

---

### Reference Files

```groovy
reference_genome_path = ""
reference_genome_index_path = ""
reference_genome_dict_path = ""
reference_genome_gtf = ""
reference_denylist_path = ""
```

Reference genome files used across alignment and variant calling steps.

---

### Variant Calling Reference Files

```groovy
variants_snp_path = ""
variants_indels_path = ""
```

VCFs of known SNPs and INDELs used by GATK for base recalibration and variant filtering.

---

### Tool Paths

```groovy
snpeff_jar_path = "path/to/snpEff.jar"
snpeff_config_path = "path/to/snpEff.config"
arriba_known_fusions = "path/to/known_fusions.tsv"
arriba_blacklist = "path/to/blacklist.tsv"
vep_cache_dir = "path/to/vep/cache"
```

Paths to annotation tools like SnpEff, VEP, and Arriba.

---

### Pipeline Behavior Flags

```groovy
concatenate = false
only_qc = false
only_star = false
skip_star = false
run_fusion = true
```

Flags to enable or skip specific modules:
- `only_qc`: Runs only QC steps and produces a multiqc report
- `skip_star`: Skips STAR alignment if BAMs are already aligned and bams sample sheet or star output folder is given in "aligned_bam_samplesheet = path/to/bam or aligned_bam_folder = path/to/star_folder(in previous run)       = null 
- `run_fusion`: Enables fusion detection

---

###  Fastp Trimming Configuration

```groovy
trim_reads_length_required = 75
trim_reads_cut_window_size = 6
trim_reads_cut_mean_quality = 25
```

Fastp parameters for trimming low-quality bases and short reads.

---

###  STAR Configuration

```groovy
star_alignSJoverhangMin = 12
star_outFilterMismatchNmax = 800
star_chimScoreDropMax = 40
```

Tuning STAR alignment sensitivity and chimeric read filtering.

---

### Interval List Processing

```groovy
scatterintervals = true
scatter_count = 3
```

If enabled, genomic intervals are split to allow parallel GATK processing.

---

###  GATK MarkDuplicates Settings

```groovy
remove_duplicates = false
validation_stringency = "STRICT"
```

Controls behavior of duplicate marking.

---

###  GATK Variant Filtering Thresholds

```groovy
gatk_vf_qual_filter = 30.0
gatk_vf_qd_filter = 2.0
gatk_vf_fs_filter = 60.0
gatk_vf_mq_filter = 30.0
gatk_vf_sor_filter = 4.0
gatk_vf_read_pos_filter = -2.0
gatk_vf_baseq_filter = -2.0
```

Sets hard thresholds for filtering out low-confidence variants.

---

###  Annotation Parameters

```groovy
genomedb = 'GRCh38.86'
ensembl_release = "113"
genome_assembly = "GRCh38"
species = "homo_sapiens"
```

Configures SnpEff and VEP for appropriate genome version and species.

---

###  Arriba Visualization

```groovy
protein_domains = ""
cytobands = ""
```

Optional extra data files for better fusion gene visualizations.

---

## 2. `process {}`: Resource Profiles by Label

Each process in the pipeline can be tagged with a `label`, which maps to specific resource profiles:

```groovy
withLabel: 'process_low' {
    cpus = 2
    memory = '8GB'
    time = '4h'
}
withLabel: 'process_medium' {
    cpus = 6
    memory = '36GB'
    time = '8h'
}
withLabel: 'process_high' {
    cpus = 12
    memory = '40GB'
    time = '16h'
}
withLabel: 'process_very_high' {
    cpus = 16
    memory = '200GB'
    time = '24h'
}
```

This helps in tailoring compute usage for different steps like alignment, variant calling, and annotation.

---

## 3. `process.withName`: Custom Tool Flags

For individual processes, extra arguments can be injected:

```groovy
withName: 'ANNOTATEVARIANTS_VEP' {
    ext.args = [
        '--fork 4',
        '--everything'
    ].join(' ').trim()
}
```

This sets `--fork 4` and `--everything` for the VEP annotation process.

---

##  Summary

| Feature | Description |
|--------|-------------|
| `params {}` | Central place for user-supplied paths, behavior flags, and tool configs |
| `process.withLabel` | Defines CPU, memory, and runtime for tagged steps |
| `process.withName` | Adds custom options for specific tools like VEP |
| Reproducibility | All changes are external to pipeline code for easier version control |

---

This file makes your pipeline **flexible**, **configurable**, and **production-ready** across HPC, cloud, or local environments.


# ⚠️ Appendix B: Chromosome Naming Convention – Important Caution

## Why This Matters

This pipeline expects all reference files to follow a **consistent chromosome naming convention**. Specifically, the pipeline is configured to work with **numeric chromosome names** such as:

```
22, X, Y, MT
```

Not UCSC-style names like:

```
chr22, chrX, chrY, chrM
```

---

## ❗ Potential Issues If Names Are Inconsistent

If the naming convention is inconsistent across files, it may result in:

-  **No variants being called**
-  **Empty VCF outputs**
-  **Errors during annotation or alignment**
-  **Contig mismatch errors** from tools like GATK, STAR, or bcftools

---

##  What You Should Check

Ensure all of the following files use **numeric chromosome names only**:

| File Type | Description |
|-----------|-------------|
| `*.fa / *.fasta` | Reference genome |
| `*.dict` | Reference dictionary |
| `*.gtf` | Gene annotation |
| `*.vcf` | Known variant files (e.g., dbSNP, Mills) |
| `*.bed / *.interval_list` | Target regions or intervals |

---

##  How to Convert Chromosome Names

You can use simple tools like:

```bash
# Example: Convert 'chr22' to '22' in a VCF file
sed 's/^chr//' input.vcf > output_cleaned.vcf

# For GTF files
sed 's/^chr//' input.gtf > cleaned.gtf

# For FASTA headers
sed -i 's/^>chr/>/' reference.fasta
```

Alternatively, tools like `bcftools annotate`, `picard LiftoverVcf`, or `GATK LiftoverIntervalList` may be required for more complex transformations.

---

##  Final Reminder

Please **ensure consistent chromosome naming before starting the workflow**. This is a **common cause of silent failures** or misinterpreted results.

If you're unsure, you can inspect chromosome headers with:

```bash
grep '^>' reference.fasta        # For FASTA
zcat input.vcf.gz | grep -v '^#' | cut -f1 | sort | uniq
head -n 5 annotation.gtf
```

---

Stay consistent, and the pipeline will thank you 
```

---
