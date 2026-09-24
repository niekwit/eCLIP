# Snakemake workflow: `eCLIP`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow for eCLIP sequencing data analysis, based on the [ENCODE eCLIP pipeline](https://www.encodeproject.org/pipelines/ENCPL357ADL/) (eCLIP-seq Processing Pipeline v2.2, Yeo lab, UCSD). Supports single-end and paired-end (inline barcodes) reads.

---

## Workflow overview

```
FASTQ reads (single-end or paired-end, auto-detected)
    │
    ▼
FastQC
    │
    ├── single-end: UMI extraction (umi_tools extract)
    └── paired-end: inline barcode demultiplexing + UMI extraction (eclipdemux)
    │
    ▼
Adapter trimming, 2 rounds (cutadapt; round 2 removes double ligation events) + FastQC
    │
    ▼
fastq-sort → STAR mapping to repeat elements → keep unmapped reads → fastq-sort
    │
    ▼
STAR mapping to genome (unique mappers only)
    │
    ├── single-end: umi_tools dedup
    └── paired-end: barcodecollapsepe.py, merge of both barcodes, keep read 2 only
    │
    ├── RPM normalised strand specific bigWig files
    │
    ▼
CLIPper peak calling (IP samples)
    │
    ▼
Input normalisation (size-matched input: Fisher exact/chi-square test, log2 fold change)
    │
    ▼
Peak compression → blacklist filtering → narrowPeak / bigBed / total entropy
    │
    ▼
Replicates of the same condition: entropy ranked IDR → reproducible peaks (narrowPeak / bigBed)
    │
    ▼
MultiQC (FastQC, cutadapt, STAR)
```

Differences with the ENCODE pipeline:

* The repeat element filter uses a free substitute for RepBase (which is not freely available): the curated Dfam consensus sequences of the species plus the rDNA repeating unit. The repeat-mapped BAM files are not kept.
* Current versions of the tools (STAR, cutadapt, umi_tools) are used instead of the versions in the SOP. IDR 2.0.2 can not be installed anymore, so IDR 2.0.4.2 is used. Yeo lab scripts that are not on conda (CLIPper, eclipdemux, and the perl/python scripts for input normalisation and IDR) are used at fixed commits.
* Genome and annotation are from GENCODE (hg38: v29, hg19: v19, mm10: vM25), which is what CLIPper's built-in annotations are based on.

---

## Deployment

### 1. Install Snakemake

Install Snakemake (≥ 8.25.5) using conda/mamba:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake=8.25.5
conda activate snakemake
```

### 2. Clone the workflow

```bash
git clone https://github.com/<owner>/<repo>.git
cd <repo>
```

### 3. Place FASTQ files

Copy or symlink your FASTQ files into the `reads/` directory (see [config/README.md](config/README.md)).

**Single-end**:

```
reads/{sample}.fastq.gz
```

**Paired-end**:

```
reads/{sample}_R1_001.fastq.gz
reads/{sample}_R2_001.fastq.gz
```

> **Note:** IP sample names must contain only alphanumeric characters and underscores, and end with `_` followed by the replicate number (e.g. `RBFOX2_1`). Size-matched input libraries must be present in `reads/` as well.

### 4. Configure

Describe all libraries in `config/samples.csv` and adjust `config/config.yaml`, see [config/README.md](config/README.md).

### 5. Run

```bash
snakemake --use-conda --cores 32
```

All reference files (genome, GENCODE annotation, blacklist, repeat elements) are downloaded and the STAR indices are built on the first run.

---

## Output

| Directory                                   | Contents                                                                                  |
| ------------------------------------------- | ----------------------------------------------------------------------------------------- |
| `results/qc/multiqc.html`                   | MultiQC report                                                                            |
| `results/mapped/{sample}.bam`               | PCR duplicate removed, uniquely mapped reads (paired-end: read 2 only)                    |
| `results/bigwig/`                           | RPM normalised strand specific read density (`.norm.pos.bw`, `.norm.neg.bw`)              |
| `results/clipper/`                          | CLIPper peak clusters                                                                     |
| `results/peaks/{sample}/`                   | Input normalised peaks: `.normed.bed(.full)`, `.normed.compressed.bed`, final blacklist filtered `.peaks.bed`, `.peaks.narrowPeak`, `.peaks.bb`, `.total_entropy.txt` |
| `results/idr/{sample1}_vs_{sample2}/`       | IDR analysis and reproducible peaks (`.reproducible_peaks.bed`, `.custombed`, `.narrowPeak`, `.bb`) |

Columns of the peak BED files: chromosome, start, end, -log10(p-value), log2 fold change over input, strand.

---

## Reference

Van Nostrand, E.L., Pratt, G.A., Shishkin, A.A. et al. Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP). *Nat Methods* 13, 508–514 (2016).

Yeo lab eCLIP pipeline: https://github.com/YeoLab/eclip
