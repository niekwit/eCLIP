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
* Genome and annotation are from GENCODE (hg38: v29, hg19: v19, mm10: vM25), which is what CLIPper's built-in annotations are based on (see below).

### Why GENCODE v29 for hg38?

CLIPper does not accept an arbitrary GTF: it ships its own pre-built annotation files for a fixed set of `--species` values, and the genome/GTF this workflow downloads has to match one of them exactly. For hg38, CLIPper also has a newer `GRCh38_v40` (GENCODE v40) option, but `v29e` is used by default because it is the exact annotation ENCODE's own eCLIP pipeline uses — this is what the [validation](validation/README.md) results are checked against. Set `gencode_release: "40"` to use `GRCh38_v40` instead (see below). For hg19 and mm10, `hg19` and `mm10v25` are the newest annotations CLIPper provides for those genomes (no newer or ENCODE-specific variant exists there), and mm39 is not supported at all because CLIPper has no mm39 annotation.

How much does this matter? Genome coordinates are identical across GENCODE releases (same GRCh38 assembly); only the annotation on top of it changes:

| | v29 (2018) | v40 (2022) | v50 (2024, latest) |
| --- | --- | --- | --- |
| Protein-coding genes | 19,940 | 19,988 | 20,107 |
| Protein-coding transcripts | 83,129 | 87,814 | 278,455 |
| lncRNA genes | 7,635 | 17,748 | 34,866 |
| Total genes | 58,721 | 61,544 | 78,733 |

* **Protein-coding genes:** essentially unchanged (+0.8% from v29 to v50). For an mRNA-binding protein, the gene set CLIPper works with is nearly the same in v29 as in the latest release.
* **Protein-coding transcript isoforms:** v29 to v40 barely changed (+6%), but v40 to v50 more than tripled (GENCODE's newer releases add many more long-read-supported alternative splice isoforms per gene). CLIPper classifies peaks by transcript region (exon/intron/UTR/proximal vs. distal intron), so a much richer isoform set can shift which region a peak is assigned to, without necessarily changing whether a peak is called there.
* **lncRNAs:** more than 4x more lncRNA genes are annotated now than in v29. This is the real gap: for an RBP that binds lncRNAs, v29e is missing roughly half of today's annotated lncRNA loci.

In short: staying on v29e for ENCODE parity costs almost nothing for protein-coding mRNA analyses, but is a genuine limitation for lncRNA-focused studies or fine isoform-level peak assignment.

### Using a newer GENCODE release

Set `gencode_release` in `config/config.yaml` to any human GENCODE release number (e.g. `"50"` for the latest) to use it instead of v29. Setting it to `"40"` uses CLIPper's other built-in annotation (`GRCh38_v40`) directly, no build step needed. For any other release, CLIPper has no built-in annotation, so one is built automatically from the downloaded GTF (`rule build_clipper_annotation`, `workflow/scripts/build_clipper_annotation.py`):

* one representative transcript per gene (the one with the largest genomic span, the same rule CLIPper's own annotation-building code uses) provides the gene coordinates and the `mrna_length`/`premrna_length` values CLIPper needs for its statistical test;
* exons of every transcript of a gene are merged into a non-overlapping list, used for splice-aware read assignment.

This replicates the two files CLIPper's ["Supporting additional species"](https://github.com/YeoLab/clipper/wiki/Supporting-additional-species) wiki page describes, and is passed to CLIPper with `--datadir` rather than modifying the installed package. It only runs when `gencode_release` is set to something other than `"29"` or `"40"`, and only builds these two files, not the full built-in CLIPper data directory (no `_genes.bed`, intron/UTR/poly-A region files etc., which CLIPper's core peak caller does not require, only its separate `clip_analysis` annotation tool does).

Using a non-default `gencode_release` is **not validated against ENCODE** the way v29e is (see [validation/](validation/README.md)): peak positions and counts will differ, both because of the annotation itself and because CLIPper's per-gene statistical test uses the picked representative transcript's length, which can differ between releases.

---

## Repeat elements and transposable elements

As in the ENCODE pipeline, reads from repeat elements are **removed** before the genome analysis, and by default transposable element (TE) derived reads are not analysed separately:

1. **Repeat element filter.** After adapter trimming, reads are mapped with STAR to a repeat element reference (`resources/{genome}_repeat_elements.fa`), which is built on the first run: the curated consensus sequences of the species from [Dfam](https://www.dfam.org) (for human ~1,400 families, e.g. LINE-1 (L1HS, L1PA, L1M, ...), Alu, SVA, LTR/ERV, DNA transposons and several small non-coding RNAs) plus the rDNA repeating unit (NCBI). Mapping is end-to-end, and a read may map to up to 30 places (`--outFilterMultimapNmax 30`), so reads from multi-copy elements are caught. **Reads that map are discarded by default; only the reads that do not map continue to the genome mapping** (unless `te_repeats.enabled` is set, see below). The STAR log of this step (`results/star/repeats/`, in the MultiQC report) gives the fraction of reads that were removed.
2. **Genome mapping.** The remaining reads are mapped to the genome, **keeping only reads with a unique alignment** (`--outFilterMultimapNmax 1`, as ENCODE). Reads from multi-copy TEs that were not caught by the repeat filter are lost in this step.
3. **Peak calling.** CLIPper, input normalisation and IDR use only these unique, PCR-duplicate removed reads. There is no masking or annotation of these peaks with TEs (RepeatMasker), the only region filter is the ENCODE eCLIP blacklist.

What this means for TEs by default:

* Reads that resemble the consensus sequence of a TE family are removed and are not counted anywhere by the main pipeline.
* Reads from TE copies that are diverged from the consensus (typically old families such as L2, MIR or old L1M) and that map to a single genomic position are retained, so peaks in TEs can appear in the results. Peaks in TE-derived sequence therefore cover only part of the TE-derived signal.
* Dfam consensus sequences are used instead of the RepBase sequences of ENCODE (RepBase is not freely available), so the reads that are removed are not exactly the same as in ENCODE.

### Locus-resolved TE binding analysis (optional)

Setting `te_repeats.enabled: True` in the config adds an opt-in analysis that answers a question the main pipeline cannot: *which individual TE copy is bound*, not just whether TE-derived reads exist. It re-uses the reads discarded in step 1 above (the repeat-mapped BAM is kept instead of deleted) rather than duplicating the repeat filter:

1. Reads that mapped to the repeat consensus reference (step 1) are re-aligned to the full genome with the same uniqueness requirement as the main genome mapping (`--outFilterMultimapNmax 1`). A read fully internal to a repeat copy still multi-maps genome-wide and is dropped here, exactly as in step 2 above; only reads with a unique anchor in the genome — typically a TE-to-flanking-sequence readthrough junction — survive, which is what makes assigning a specific genomic TE locus possible.
2. Surviving reads are PCR-duplicate removed with the same method as the main pipeline, then intersected with individual TE copies from UCSC RepeatMasker (downloaded automatically) to assign each read to one locus. `te_repeats.require_family_match` (default on) additionally requires that locus's RepeatMasker family to agree with the family the read hit in step 1, as a best-effort cross-check (Dfam and RepeatMasker family names mostly, but not always, agree).
3. Per IP sample: Fisher's exact test (BH-corrected) of read counts per locus and per family, IP vs size-matched input, normalised against each library's total confidently-assigned read count (not whole-library size) — `results/te_repeats/{sample}.locus_enrichment.tsv` and `.family_enrichment.tsv`. Replicate pairs additionally get a Spearman correlation of log2FC (`results/te_repeats/{s1}_vs_{s2}.reproducibility.txt`), used instead of IDR since IDR is designed for genomic-interval peak calling, not a fixed set of categorical loci/families.

This roughly doubles genome-mapping and dedup work (a second alignment pass on the repeat-mapped reads) and keeps the repeat-mapped BAM instead of deleting it (~1 GB per library in this project's ENCODE validation runs), so it is off by default.

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

| sample           | control          | adapter   |
| ---------------- | ----------------- | --------- |
| RBFOX2_1         | RBFOX2_input_1    | InvRil19  |
| RBFOX2_2         | RBFOX2_input_2    | InvRil19  |
| RBFOX2_input_1   |                    | InvRil19  |
| RBFOX2_input_2   |                    | InvRil19  |

Paired-end example (extra columns with the inline barcode IDs; `NIL` for barcode-less input):

| sample           | control          | barcode_a | barcode_b |
| ---------------- | ----------------- | --------- | --------- |
| RBFOX2_1         | RBFOX2_input_1    | A01       | B06       |
| RBFOX2_2         | RBFOX2_input_1    | C01       | D8f       |
| RBFOX2_input_1   |                    | NIL       | NIL       |

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

The workflow was validated on two ENCODE eCLIP experiments (U2AF2, paired-end): [ENCSR202BFN](https://www.encodeproject.org/experiments/ENCSR202BFN/) in HepG2 and [ENCSR893RAV](https://www.encodeproject.org/experiments/ENCSR893RAV/) in K562. For both, the numbers of usable reads and peaks are within about 2% of ENCODE, 89-96% of the peaks overlap, and fold changes and p-values of overlapping peaks correlate with r = 0.93-0.99. Details, the exact steps and the comparison script are in [validation/](validation/README.md).

---

## Reference

Van Nostrand, E.L., Pratt, G.A., Shishkin, A.A. et al. Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP). *Nat Methods* 13, 508–514 (2016).

Yeo lab eCLIP pipeline: https://github.com/YeoLab/eclip
