
#  RNA-Seq Variant Calling & Gene Fusion Detection Pipeline

This document describes the step-by-step setup and validation required before executing the `buildreference.nf` workflow. Proper preparation of reference files ensures smooth alignment, variant calling, and annotation.

---

## ✅ Step 1: Input Reference Files Required for `buildreference.nf`

Before running the pipeline, make sure the following reference files and folders are available and correctly formatted.

These files are essential for building genome indexes, annotations, and known variant databases used by downstream processes.

###  Required Files & Their Purpose

| File Name                   | Description                                                                 |
|-----------------------------|-----------------------------------------------------------------------------|
| `genome.fa`                 | Reference genome in FASTA format                                            |
| `genome.fa.fai`             | FASTA index file ir it is  generated using `samtools faidx`                 |
| `genome.dict`               | Sequence dictionary generated using GATK `CreateSequenceDictionary`         |
| `annotations.gtf`           | Gene annotation file in GTF format (used for STAR and SnpEff)               |
| `dbsnp.vcf.gz`              | dbSNP VCF file containing known SNPs                                        |
| `dbsnp.vcf.gz.tbi`          | Tabix index for the dbSNP VCF                                               |
| `known_indels.vcf.gz`       | VCF file of known INDELs (e.g., Mills & 1000G)                              |
| `known_indels.vcf.gz.tbi`   | Tabix index for the INDEL VCF                                               |
| `snpeff_cache/`             | Folder containing the SnpEff cache for your reference genome                |
| `vep_cache/`                | Folder containing the VEP cache database                                    |
| `star_index/`               | STAR genome index folder (optional if already generated)                    |

>  **Unzipping** compressed reference files (e.g., `.fa.gz`, `.gtf.gz`) is **recommended** for compatibility and performance:
> ```bash
> gunzip genome.fa.gz
> gunzip annotations.gtf.gz
> ```

---

##  Step 2: Ensure Consistent Chromosome Naming (Ensembl Style)

Before using these files in the pipeline, ensure they **all follow the same chromosome naming convention**.

>  Mismatched styles like `chr1` vs `1` can cause pipeline failures or misalignments.

###  This Pipeline Uses: **Ensembl-style naming**

```
1, 2, 3, ..., 22, X, Y, MT
```

>  UCSC-style (`chr1`, `chrX`, etc.) is not supported unless manually converted to Ensembl style.

---

###  Files to Check

Check naming consistency in:

| File Type              | Example File               |
|------------------------|----------------------------|
| Reference genome FASTA | `genome.fa`                |
| Genome dictionary      | `genome.dict`              |
| Genome index (`.fai`)  | `genome.fa.fai`            |
| GTF annotation file    | `annotations.gtf`          |

>  SNP and INDEL VCF files do **not** need to be checked — the pipeline will auto-convert chromosome naming when merging them.

---

###  How to Check Naming

####  1. FASTA
```bash
zgrep '^>' genome.fa.gz | cut -d' ' -f1
```

####  2. Genome Index
```bash
head genome.fa.fai
```

####  3. Genome Dictionary
```bash
grep '^@SQ' genome.dict
```

####  4. GTF File
```bash
zcat annotations.gtf.gz | cut -f1 | grep -v '^#' | sort -u
```

>  All must use either `1`, `2`, … `X`, `Y` or `chr1`, `chr2`, … — not mixed.

![Chromosome Naming Convention](Documentation/chromosome naming.png)


---

##  Once Steps 1 & 2 are verified, you're ready to proceed to **Step 3: Run `buildreference.nf`**.
```

---

```
#  Step 3: Run the `buildreference.nf` Pipeline

After verifying that your reference files are correctly formatted (Step 1) and consistently named (Step 2), you can now proceed to run the **reference-building workflow**.

This step generates all necessary indexes and annotation resources used throughout the pipeline.

---

##  Command to Run

```bash
nextflow run main.nf -c nextflow_ref_main.config --build_references -profile singularity
```

---

##  Explanation of Parameters

| Option/Flag                   | Description                                                                |
|------------------------------|----------------------------------------------------------------------------|
| `main.nf`                    | The main Nextflow script                                                   |
| `-c nextflow_ref_main.config`| Loads the config file with paths to reference files                       |
| `--build_references`         | Tells the pipeline to execute the reference-building workflow only         |
| `-profile singularity`       | Uses the Singularity profile for containerized execution                   |

---

##  Output

Once complete, this workflow will produce:

- Genome indexes (FASTA index, dict, STAR index)
- GTF/GFF checks
- Known variant files
- SnpEff database verification
- VEP cache confirmation

All outputs will be stored in the folders defined in your `nextflow_ref_main.config`.


---

##  Step 4: Run the Full Pipeline with Samples

Once all references are built successfully, you're ready to run the **complete RNA-seq variant calling and gene fusion workflow** on your sample data.

### ▶️ Command to Execute

```bash
nextflow run main.nf -c nextflow_main.config -profile singularity
```

---

###  Explanation

| Option/Flag                   | Description                                                         |
|------------------------------|---------------------------------------------------------------------|
| `main.nf`                    | The main Nextflow pipeline script                                   |
| `-c nextflow_main.config`    | Configuration file containing paths to input samples and references |
| `-profile singularity`       | Enables Singularity containers for reproducible environments         |

---

###  What Happens

This runs the **full workflow**, including:
- Quality control (FastQC, fastp)
- Alignment (STAR)
- Variant calling (GATK)
- Variant annotation (SnpEff, VEP)
- Fusion detection (Arriba)
- MultiQC reporting

---

 Outputs will be stored in structured directories as defined in your config, with logs, results, and version tracking.


```
