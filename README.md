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

## Usage

### 1. Install Snakemake

Install Snakemake (≥ 8.25.5) using conda/mamba. A recent conda (≥ 24.7) is needed for `--use-conda`:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake=8.25.5
conda activate snakemake
```

### 2. Get the workflow

```bash
git clone https://github.com/niekwit/eCLIP.git
cd eCLIP
```

Run the workflow from this directory (or from a project directory that has the same `config/` and `reads/` layout; Snakemake finds `workflow/Snakefile` automatically).

### 3. Add your reads

Copy or symlink the FASTQ files of **all** libraries (IP and size-matched input) into `reads/`. The layout is auto-detected from the file names.

| Layout       | File names                                                             |
| ------------ | ---------------------------------------------------------------------- |
| Single-end   | `reads/{sample}.fastq.gz`                                              |
| Paired-end   | `reads/{sample}_R1_001.fastq.gz` and `reads/{sample}_R2_001.fastq.gz`  |

Paired-end means ENCODE-style eCLIP where read 1 starts with an inline barcode and read 2 starts with the UMI. Single-end means seCLIP, where the first 10 nt of the read are the UMI.

> **FASTQ files from the ENCODE portal are already demultiplexed** (the inline barcodes are removed from read 1 and the UMI is the first part of the read name, e.g. `@CAAAA:HWI-D00611:...`). For these files, set `demultiplexed: True` in `config/config.yaml`: the demultiplexing step is skipped and each library is processed as a whole. The barcode IDs in `samples.csv` are still needed, as they define the adapters that are trimmed. Leave it `False` for raw reads.

### 4. Describe your samples

Edit `config/samples.csv`: one row per library, with the matching size-matched input for every IP sample.

Single-end example:

```csv
sample,control,adapter
RBFOX2_1,RBFOX2_input_1,InvRil19
RBFOX2_2,RBFOX2_input_2,InvRil19
RBFOX2_input_1,,InvRil19
RBFOX2_input_2,,InvRil19
```

Paired-end example (extra columns with the inline barcode IDs; `NIL` for barcode-less input):

```csv
sample,control,barcode_a,barcode_b
RBFOX2_1,RBFOX2_input_1,A01,B06
RBFOX2_2,RBFOX2_input_1,C01,D8f
RBFOX2_input_1,,NIL,NIL
```

Rules for sample names (checked when the workflow starts):

* Only letters, numbers and `_`.
* IP samples end with `_` and the replicate number (`RBFOX2_1`, `RBFOX2_2`). The rest of the name is the condition: replicates of the same condition are compared with IDR.
* Sample names can not contain `_vs_`.

All columns are explained in [config/README.md](config/README.md).

### 5. Configure

Set the genome in `config/config.yaml` (`hg38`, `hg19` or `mm10`). The other settings default to the ENCODE pipeline values and usually do not need to be changed. Also adjust the `resources` section (CPUs and time in minutes per step) to your machine.

### 6. Check and run

Do a dry-run first. This checks the samples file, the read files, and shows the jobs that will be run:

```bash
snakemake --use-conda --cores 32 -n
```

Then run the workflow:

```bash
snakemake --use-conda --cores 32
```

On a cluster, use your Snakemake profile or an executor plugin (here `snakemake-executor-plugin-slurm`), for example:

```bash
snakemake --use-conda --executor slurm --jobs 50 --default-resources slurm_account=<account> slurm_partition=<partition>
```

Generate a report with the results afterwards:

```bash
snakemake --report report.zip
```

What to expect on the first run:

* Conda environments are created (a few minutes each; CLIPper is built from source).
* The reference files are downloaded and indexed: GENCODE genome and annotation, ENCODE blacklist, the Dfam repeat elements, and the Yeo lab scripts. Two STAR indices (genome and repeat elements) are built; the genome index needs roughly 30-40 GB of RAM for human/mouse (the workflow requests 40 GB). This only happens once, in `resources/`.
* CLIPper is the slowest step: it processes the whole genome annotation and takes hours per IP sample for a full dataset. Give it enough CPUs (`resources: clipper: cpu` in the config).
* A warning is printed when a condition has only one replicate: peaks are still called for its samples, but no IDR analysis is done.

---

## Expected output

For the example samples above, a successful run gives (single-end and paired-end have the same structure):

```
results/
├── qc/
│   ├── multiqc.html                      # MultiQC report: FastQC, cutadapt, STAR
│   ├── fastqc/  cutadapt/                # per sample QC data
│   └── umi_tools/  (single-end)          # PCR duplicate stats
│       demux/  barcode_collapse/  (paired-end)
├── mapped/
│   ├── RBFOX2_1.bam(.bai)                # unique, PCR-duplicate removed reads (paired-end: read 2 only)
│   └── RBFOX2_1.readnum.txt              # number of mapped reads (used for input normalisation)
├── bigwig/
│   ├── RBFOX2_1.norm.pos.bw              # RPM normalised read density, + strand
│   └── RBFOX2_1.norm.neg.bw              # RPM normalised read density, - strand (negative values)
├── clipper/
│   └── RBFOX2_1.peakClusters.bed         # CLIPper peak clusters (IP samples)
├── peaks/RBFOX2_1/
│   ├── RBFOX2_1.normed.bed(.full)        # peaks with p-value and fold change over input
│   ├── RBFOX2_1.normed.compressed.bed    # overlapping peaks merged
│   ├── RBFOX2_1.peaks.bed                # FINAL peaks (blacklist filtered)
│   ├── RBFOX2_1.peaks.narrowPeak         # same, narrowPeak format
│   ├── RBFOX2_1.peaks.bb                 # same, bigBed format
│   └── RBFOX2_1.total_entropy.txt        # total relative information content of significant peaks
├── idr/RBFOX2_1_vs_RBFOX2_2/
│   ├── RBFOX2_1_vs_RBFOX2_2.reproducible_peaks.bed         # FINAL reproducible peaks
│   ├── RBFOX2_1_vs_RBFOX2_2.reproducible_peaks.custombed   # per replicate values of each peak
│   ├── RBFOX2_1_vs_RBFOX2_2.reproducible_peaks.narrowPeak
│   ├── RBFOX2_1_vs_RBFOX2_2.reproducible_peaks.bb
│   ├── RBFOX2_1_vs_RBFOX2_2.idr.out(.png)                  # IDR output and diagnostic plot
│   └── entropy/, *.idr_peaks.normed.bed(.full)              # intermediate files
└── star/                                 # STAR logs (repeat element and genome mapping)
resources/                                # reference files and STAR indices (created on first run)
logs/                                     # log file of every job
```

Files you will normally use:

| File                                                  | What it is                                                                         |
| ----------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `results/peaks/{sample}/{sample}.peaks.bed`           | Input normalised, blacklist filtered peaks of one replicate                        |
| `results/idr/{s1}_vs_{s2}/*.reproducible_peaks.bed`   | Reproducible peaks between two replicates (use these as the final binding sites)   |
| `results/bigwig/*.bw`                                 | Tracks for a genome browser (+ strand positive, - strand negative)                 |
| `results/qc/multiqc.html`                             | Quality control (adapter trimming, repeat and genome mapping rates)                |

Columns of `.peaks.bed` and `.reproducible_peaks.bed`: chromosome, start, end, -log10(p-value) over input, log2 fold change over input, strand. In the `.narrowPeak` files the score is 1000 for significant peaks (`peaks: l10p` and `l2fc` in the config, both 3 by default) and 200 for other peaks. In `.reproducible_peaks.bed` the p-value is the minimum of the two replicates and the fold change the geometric mean.

### Quick quality check

* MultiQC: most reads should have adapters trimmed after round 1, and the repeat element mapping rate (STAR `repeats`) is typically high for eCLIP; the genome (`genome`) unique mapping rate of the remaining reads should be high (a low rate indicates a wrong genome or a poor library).
* `results/mapped/*.readnum.txt`: number of usable reads after PCR duplicate removal (ENCODE recommends roughly 1 million or more for IP samples).
* `results/idr/*/*.idr.out.png`: replicate rank plot; reproducible replicates show a clear diagonal.

### Troubleshooting

* *"Following files not found"* at start: file names in `reads/` do not match the `sample` column.
* *IDR fails*: it needs a reasonable number of peaks in both replicates (hundreds); very shallow or failed libraries do not have these.
* *STAR runs out of memory*: increase `mem_mb` of the STAR rules (`star_index_genome`, `star_repeats`, `star_genome` in `workflow/rules/`) or use fewer parallel jobs.
* Failed jobs write their errors to `logs/<tool>/<sample>.log`.

---

## Validation

The workflow was validated on the ENCODE experiment [ENCSR202BFN](https://www.encodeproject.org/experiments/ENCSR202BFN/) (U2AF2 eCLIP in HepG2, paired-end): the numbers of usable reads and peaks are within about 2% of ENCODE, 89-95% of the peaks overlap, and fold changes and p-values of overlapping peaks correlate with r = 0.94-0.99. Details, the exact steps and the comparison script are in [validation/](validation/README.md).

---

## Reference

Van Nostrand, E.L., Pratt, G.A., Shishkin, A.A. et al. Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP). *Nat Methods* 13, 508–514 (2016).

Yeo lab eCLIP pipeline: https://github.com/YeoLab/eclip
