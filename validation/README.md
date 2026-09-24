# Validation on ENCODE data

The workflow was run on a real ENCODE eCLIP experiment and the results were compared with the peaks and alignments that ENCODE deposited for the same experiment.

**Experiment:** [ENCSR202BFN](https://www.encodeproject.org/experiments/ENCSR202BFN/), eCLIP of U2AF2 in HepG2 cells (paired-end, 2 biological replicates, hg38), with the size-matched input experiment [ENCSR049GND](https://www.encodeproject.org/experiments/ENCSR049GND/). The ENCODE alignments of this experiment were derived from these FASTQ files with fastqc, cutadapt, STAR and `barcode_collapse_pe.py`, i.e. the pipeline that this workflow reproduces.

## Summary

| Comparison with ENCODE                                   | Result                                                                                       |
| -------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| Usable reads (unique, PCR-duplicate removed), replicate 1 | 6,989,114 (ENCODE: 6,995,208, -0.09%)                                                        |
| Usable reads, replicate 2                                | 4,735,584 (ENCODE: 4,815,229, -1.7%)                                                         |
| Peaks per replicate (after input normalisation)          | 256,450 and 183,979 (ENCODE: 255,064 and 185,196)                                            |
| Significant peaks (-log10 p ≥ 3 and log2 fold change ≥ 3) | 17,741 and 15,394 (ENCODE: 17,406 and 15,759)                                                |
| IDR reproducible peaks                                   | 10,543 (ENCODE: 10,732)                                                                      |
| Peak overlap (same strand)                               | 89-95% of the peaks of either analysis overlap a peak of the other                          |
| Agreement of values in overlapping peaks                 | log2 fold change r = 0.94-0.99, -log10(p-value) r = 0.94-0.96 (Spearman 0.97-0.99)           |

The peaks that do not overlap are weak calls (median -log10 p-value 0.2-0.5, only 0.5-1% pass the significance cutoffs) in both directions. 82% of the ENCODE reproducible peaks that are not in the IDR set of the workflow overlap a significant peak of at least one of its replicates, so these are borderline IDR calls.

## Data

| Library                          | ENCODE files (R1, R2)           | Inline barcodes | Reads (pairs) |
| -------------------------------- | ------------------------------- | --------------- | ------------- |
| IP, replicate 1 (`U2AF2_1`)      | ENCFF565HKG, ENCFF559BLA        | A01, B06        | 10.5 M        |
| IP, replicate 2 (`U2AF2_2`)      | ENCFF912MAK, ENCFF440IJI        | A03, G07        | 8.8 M         |
| Size-matched input (`U2AF2_input_1`) | ENCFF156ZDE, ENCFF939YLN    | none (`NIL`)    | 19.5 M        |

ENCODE files used for the comparison: peaks of replicate 1 (ENCFF067JAD), replicate 2 (ENCFF243RQR) and IDR reproducible peaks (ENCFF721PWF), all GRCh38 narrowPeak files; alignments ENCFF358STL (replicate 1) and ENCFF033XVX (replicate 2).

**The FASTQ files on the ENCODE portal are already demultiplexed**: the inline barcodes are removed from read 1 and the UMI is at the start of the read name (`@CAAAA:HWI-D00611:...`), with both barcodes of a library in the same file. The workflow therefore ran with `demultiplexed: True` (see [config/README.md](../config/README.md)).

## How it was run

```bash
mkdir -p reads encode_peaks
# download the FASTQ files and check the MD5 checksums given by ENCODE
for acc in ENCFF565HKG ENCFF559BLA ENCFF912MAK ENCFF440IJI ENCFF156ZDE ENCFF939YLN; do
    curl -sL -o $acc.fastq.gz https://www.encodeproject.org/files/$acc/@@download/$acc.fastq.gz
done
# name them {sample}_R1_001.fastq.gz / {sample}_R2_001.fastq.gz in reads/ (table above)

# ENCODE peaks
for acc in ENCFF067JAD ENCFF243RQR ENCFF721PWF; do
    curl -sL -o encode_peaks/$acc.bed.gz https://www.encodeproject.org/files/$acc/@@download/$acc.bed.gz
done
```

`config/samples.csv`:

```csv
sample,control,barcode_a,barcode_b
U2AF2_1,U2AF2_input_1,A01,B06
U2AF2_2,U2AF2_input_1,A03,G07
U2AF2_input_1,,NIL,NIL
```

`config/config.yaml`: default config with `genome: hg38` and `demultiplexed: True`. Then:

```bash
snakemake --use-conda --cores 44
```

The whole workflow, including the download of the hg38 genome and the build of the STAR indices, took about 5 hours on 44 cores. CLIPper is by far the slowest step (about 3 hours per replicate with 20 CPUs each, running in parallel), the STAR index of hg38 takes about 35 minutes.

Comparison of the peaks:

```bash
python validation/compare_encode_peaks.py --results results --encode encode_peaks \
    --rep1 ENCFF067JAD --rep2 ENCFF243RQR --idr ENCFF721PWF --sample1 U2AF2_1 --sample2 U2AF2_2
```

## Results

### Mapping

| Step (STAR, repeat elements and genome)                 | Replicate 1 | Replicate 2 | Input      |
| ------------------------------------------------------- | ----------- | ----------- | ---------- |
| Input read pairs                                        | 10,536,763  | 8,831,751   | 19,536,752 |
| Trimmed reads mapped to repeat elements (removed)       | 17.0%       | 18.0%       | 43.1%      |
| Usable reads (unique, PCR-duplicate removed, read 2)    | 6,989,114   | 4,735,584   | 9,034,720  |

Alignment level, compared with the ENCODE BAM files (read 2): for replicate 1, 98% of the unique alignments (UMI, chromosome, position and strand) are the same. For replicate 2, 94.2% of the reads that both analyses kept have the same chromosome, position, strand and CIGAR.

### Peaks

| Peak set                          | This workflow | ENCODE  | Workflow peaks overlapping ENCODE | ENCODE peaks overlapping workflow |
| --------------------------------- | ------------- | ------- | --------------------------------- | --------------------------------- |
| Replicate 1, all peaks            | 256,450       | 255,064 | 242,179 (94.4%)                   | 242,171 (94.9%)                   |
| Replicate 1, significant peaks    | 17,741        | 17,406  | 16,609 (93.6%)                    | 16,590 (95.3%)                    |
| Replicate 2, all peaks            | 183,979       | 185,196 | 167,958 (91.3%)                   | 168,136 (90.8%)                   |
| Replicate 2, significant peaks    | 15,394        | 15,759  | 14,068 (91.4%)                    | 14,094 (89.4%)                    |
| IDR reproducible peaks            | 10,543        | 10,732  | 9,806 (93.0%)                     | 9,809 (91.4%)                     |

| Peak set                          | Overlapping pairs | log2FC Pearson | log2FC Spearman | -log10(p) Pearson | -log10(p) Spearman |
| --------------------------------- | ----------------- | -------------- | --------------- | ----------------- | ------------------ |
| Replicate 1, all peaks            | 242,179           | 0.988          | 0.987           | 0.962             | 0.986              |
| Replicate 1, significant peaks    | 16,609            | 0.949          | 0.941           | 0.952             | 0.981              |
| Replicate 2, all peaks            | 167,958           | 0.980          | 0.977           | 0.942             | 0.972              |
| Replicate 2, significant peaks    | 14,068            | 0.936          | 0.925           | 0.950             | 0.969              |
| IDR reproducible peaks            | 9,806             | 0.945          | 0.941           | 0.955             | 0.967              |

Overlap is on the same strand, and correlations are of the best overlapping ENCODE peak of each workflow peak. Remaining differences are expected: the workflow uses newer versions of the tools (STAR 2.7.11b, cutadapt 5.1, umi_tools/IDR versions, CLIPper from the current YeoLab repository), a Dfam based repeat element reference instead of RepBase, and the random choice of the retained read between PCR duplicates.

## Issues found with this data set

Running real data found two problems that the simulated test data did not show. Both are fixed in the workflow:

1. **Reads from the ENCODE portal are already demultiplexed.** The workflow assumed raw reads and would have tried to demultiplex them. This is now supported with `demultiplexed: True` (each library is processed as a whole, the barcode IDs in `samples.csv` still define the adapters).
2. **Trimming with barcodes that contain random bases (N).** Barcodes such as A03/G07 (used for replicate 2) contain 4 random bases. In the second round of adapter trimming (minimum overlap of 5, as in the ENCODE SOP), adapter chunks that start with `NNNN` match the last 5 bases of any read that ends with the next base of the chunk, so 5 real bases were removed from ~58% of the reads. Compared with ENCODE's alignments, only 42% of the reads had the same alignment and 3.8% fewer reads were kept. Second round trimming now only uses the adapter chunks without N (round 1 is unchanged, and barcodes without N such as A01/B06/C01/D8f are not affected): 94% of the reads now have the same alignment as ENCODE and the number of usable reads is within 1.7%.
