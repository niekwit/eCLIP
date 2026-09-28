# Validation on ENCODE data

The workflow was run on two real ENCODE eCLIP experiments and the results were compared with the peaks and alignments that ENCODE deposited for the same experiments.

| Experiment | Target | Cell line | Replicates | ENCODE peaks compared |
| ---------- | ------ | --------- | ---------- | --------------------- |
| [ENCSR202BFN](https://www.encodeproject.org/experiments/ENCSR202BFN/) | U2AF2 | HepG2 | 2 | ENCFF067JAD, ENCFF243RQR, ENCFF721PWF (IDR) |
| [ENCSR893RAV](https://www.encodeproject.org/experiments/ENCSR893RAV/) | U2AF2 | K562  | 2 | ENCFF182HYE, ENCFF405ETO, ENCFF290DFO (IDR) |

Both are paired-end, hg38, with a size-matched input experiment ([ENCSR049GND](https://www.encodeproject.org/experiments/ENCSR049GND/) and [ENCSR871RTL](https://www.encodeproject.org/experiments/ENCSR871RTL/) respectively). The ENCODE alignments and peaks of both experiments were derived from these FASTQ files with fastqc, cutadapt, STAR, `barcode_collapse_pe.py`, CLIPper and the input normalisation/IDR pipeline, i.e. the pipeline that this workflow reproduces.

## Summary

|                                                            | ENCSR202BFN (HepG2)                             | ENCSR893RAV (K562)                              |
| ---------------------------------------------------------- | ------------------------------------------------ | ------------------------------------------------ |
| Usable reads, replicate 1 (unique, PCR-duplicate removed)   | 6,989,114 (ENCODE: 6,995,208, -0.09%)            | 4,654,692 (ENCODE: 4,673,279, -0.40%)            |
| Usable reads, replicate 2                                   | 4,735,584 (ENCODE: 4,815,229, -1.7%)             | 6,853,134 (ENCODE: 6,882,780, -0.43%)            |
| Identical alignments (position, strand, CIGAR) of shared reads | 94-98%                                         | 99.8% (both replicates)                          |
| Peaks per replicate (after input normalisation)              | 256,450 and 183,979 (ENCODE: 255,064 and 185,196) | 165,923 and 241,065 (ENCODE: 163,877 and 237,177) |
| Significant peaks (-log10 p ≥ 3 and log2 fold change ≥ 3)    | 17,741 and 15,394 (ENCODE: 17,406 and 15,759)     | 4,815 and 4,900 (ENCODE: 4,941 and 5,185)         |
| IDR reproducible peaks                                       | 10,543 (ENCODE: 10,732)                          | 3,534 (ENCODE: 3,641)                            |
| Peak overlap (same strand)                                   | 89-95% of the peaks of either analysis overlap a peak of the other | 91-96%                     |
| Agreement of values in overlapping peaks                     | log2FC r = 0.94-0.99, -log10(p) r = 0.97-0.98 (Spearman 0.93-0.99) | log2FC r = 0.96-0.99, -log10(p) r = 0.93-0.98 (Spearman 0.93-0.99) |

The peaks that do not overlap are weak calls in both directions (median -log10 p-value close to 0, well under 1% pass the significance cutoffs). For ENCSR202BFN, 82% of the ENCODE reproducible peaks that are not in the IDR set of the workflow overlap a significant peak of at least one of its replicates, so these are borderline IDR calls.

The second experiment (K562) agrees with ENCODE at least as well as the first (HepG2) on every measure, including alignments (99.8% vs. 94-98% identical), which shows that the fixes made after the first run (see below) generalise and are not specific to one experiment or one set of barcodes.

## Data

### ENCSR202BFN (HepG2)

| Library                          | ENCODE files (R1, R2)           | Inline barcodes | Reads (pairs) |
| -------------------------------- | ------------------------------- | --------------- | ------------- |
| IP, replicate 1 (`U2AF2_1`)      | ENCFF565HKG, ENCFF559BLA        | A01, B06        | 10.5 M        |
| IP, replicate 2 (`U2AF2_2`)      | ENCFF912MAK, ENCFF440IJI        | A03, G07        | 8.8 M         |
| Size-matched input (`U2AF2_input_1`) | ENCFF156ZDE, ENCFF939YLN    | none (`NIL`)    | 19.5 M        |

Alignments used for the read-level comparison: ENCFF358STL (replicate 1), ENCFF033XVX (replicate 2).

### ENCSR893RAV (K562)

| Library                          | ENCODE files (R1, R2)           | Inline barcodes | Reads (pairs) |
| -------------------------------- | ------------------------------- | --------------- | ------------- |
| IP, replicate 1 (`U2AF2_1`)      | ENCFF567UMG, ENCFF886IXX        | A01, B06        | 11.8 M        |
| IP, replicate 2 (`U2AF2_2`)      | ENCFF270BLT, ENCFF708UWK        | C01, D8f        | 13.9 M        |
| Size-matched input (`U2AF2_input_1`) | ENCFF325JZX, ENCFF303MYY    | none (`NIL`)    | 12.9 M        |

Alignments used for the read-level comparison: ENCFF835KXL (replicate 1), ENCFF936JSP (replicate 2).

**The FASTQ files on the ENCODE portal are already demultiplexed** for both experiments: the inline barcodes are removed from read 1 and the UMI is at the start of the read name (`@CAAAA:HWI-D00611:...`), with both barcodes of a library in the same file. The workflow therefore ran with `demultiplexed: True` (see [config/README.md](../config/README.md)).

## How it was run

```bash
mkdir -p reads encode_peaks
# download the FASTQ files of the experiment and its input, and check the MD5 checksums given by ENCODE
for acc in <accessions of the tables above>; do
    curl -sL -o $acc.fastq.gz https://www.encodeproject.org/files/$acc/@@download/$acc.fastq.gz
done
# name them {sample}_R1_001.fastq.gz / {sample}_R2_001.fastq.gz in reads/ (tables above)

# ENCODE peaks (replicate 1, replicate 2, IDR), for the comparison below
for acc in <peak accessions of the summary table>; do
    curl -sL -o encode_peaks/$acc.bed.gz https://www.encodeproject.org/files/$acc/@@download/$acc.bed.gz
done
```

`config/samples.csv` (barcodes from the data tables above):

```csv
sample,control,barcode_a,barcode_b
U2AF2_1,U2AF2_input_1,A01,B06
U2AF2_2,U2AF2_input_1,<barcode_a>,<barcode_b>
U2AF2_input_1,,NIL,NIL
```

`config/config.yaml`: default config with `genome: hg38` and `demultiplexed: True`. Then:

```bash
snakemake --use-conda --cores 44
```

The whole workflow, including the download of the hg38 genome and the build of the STAR indices (shared between the two experiments, downloaded/built only once), took about 5 hours per experiment on 44 cores. CLIPper is by far the slowest step (about 3 hours per replicate with 20 CPUs each, running in parallel), the STAR index of hg38 takes about 35 minutes.

Comparison of the peaks:

```bash
python validation/compare_encode_peaks.py --results results --encode encode_peaks \
    --rep1 <accession> --rep2 <accession> --idr <accession> --sample1 U2AF2_1 --sample2 U2AF2_2
```

## Results

### ENCSR202BFN (HepG2)

#### Mapping

| Step (STAR, repeat elements and genome)                 | Replicate 1 | Replicate 2 | Input      |
| ------------------------------------------------------- | ----------- | ----------- | ---------- |
| Input read pairs                                        | 10,536,763  | 8,831,751   | 19,536,752 |
| Trimmed reads mapped to repeat elements (removed)       | 17.0%       | 18.0%       | 43.1%      |
| Usable reads (unique, PCR-duplicate removed, read 2)    | 6,989,114   | 4,735,584   | 9,034,720  |

Alignment level, compared with the ENCODE BAM files (read 2): for replicate 1, 98% of the unique alignments (UMI, chromosome, position and strand) are the same. For replicate 2, 94.2% of the reads that both analyses kept have the same chromosome, position, strand and CIGAR.

#### Peaks

| Peak set                          | This workflow | ENCODE  | Workflow peaks overlapping ENCODE | ENCODE peaks overlapping workflow |
| --------------------------------- | ------------- | ------- | --------------------------------- | --------------------------------- |
| Replicate 1, all peaks            | 256,450       | 255,064 | 242,179 (94.4%)                   | 242,171 (94.9%)                   |
| Replicate 1, significant peaks    | 17,741        | 17,406  | 16,609 (93.6%)                    | 16,590 (95.3%)                    |
| Replicate 2, all peaks            | 183,979       | 185,196 | 167,958 (91.3%)                   | 168,136 (90.8%)                   |
| Replicate 2, significant peaks    | 15,394        | 15,759  | 14,068 (91.4%)                    | 14,094 (89.4%)                    |
| IDR reproducible peaks            | 10,543        | 10,732  | 9,806 (93.0%)                     | 9,809 (91.4%)                     |

| Peak set                          | Overlapping pairs | log2FC Pearson | log2FC Spearman | -log10(p) Pearson | -log10(p) Spearman |
| --------------------------------- | ----------------- | -------------- | --------------- | ----------------- | ------------------ |
| Replicate 1, all peaks            | 242,179           | 0.989          | 0.988           | 0.978             | 0.987              |
| Replicate 1, significant peaks    | 16,609            | 0.958          | 0.948           | 0.967             | 0.986              |
| Replicate 2, all peaks            | 167,958           | 0.982          | 0.978           | 0.970             | 0.974              |
| Replicate 2, significant peaks    | 14,068            | 0.944          | 0.931           | 0.973             | 0.973              |
| IDR reproducible peaks            | 9,806             | 0.953          | 0.947           | 0.973             | 0.971              |

### ENCSR893RAV (K562)

#### Mapping

| Step (STAR, repeat elements and genome)                 | Replicate 1 | Replicate 2 | Input      |
| ------------------------------------------------------- | ----------- | ----------- | ---------- |
| Input read pairs                                        | 11,771,658  | 13,872,514  | 12,925,787 |
| Trimmed reads mapped to repeat elements (removed)       | 36.1%       | 23.4%       | 76.1%      |
| Usable reads (unique, PCR-duplicate removed, read 2)    | 4,654,692   | 6,853,134   | 1,833,273  |

Alignment level, compared with the ENCODE BAM files (read 2): 99.8% of the reads that both analyses kept have the same chromosome, position, strand and CIGAR, for both replicates.

#### Peaks

| Peak set                          | This workflow | ENCODE  | Workflow peaks overlapping ENCODE | ENCODE peaks overlapping workflow |
| --------------------------------- | ------------- | ------- | --------------------------------- | --------------------------------- |
| Replicate 1, all peaks            | 165,923       | 163,877 | 157,098 (94.7%)                   | 157,163 (95.9%)                   |
| Replicate 1, significant peaks    | 4,815         | 4,941   | 4,637 (96.3%)                     | 4,651 (94.1%)                     |
| Replicate 2, all peaks            | 241,065       | 237,177 | 227,914 (94.5%)                   | 227,916 (96.1%)                   |
| Replicate 2, significant peaks    | 4,900         | 5,185   | 4,690 (95.7%)                     | 4,726 (91.1%)                     |
| IDR reproducible peaks            | 3,534         | 3,641   | 3,374 (95.5%)                     | 3,384 (92.9%)                     |

| Peak set                          | Overlapping pairs | log2FC Pearson | log2FC Spearman | -log10(p) Pearson | -log10(p) Spearman |
| --------------------------------- | ----------------- | -------------- | --------------- | ----------------- | ------------------ |
| Replicate 1, all peaks            | 157,098           | 0.989          | 0.986           | 0.984             | 0.932               |
| Replicate 1, significant peaks    | 4,637             | 0.958          | 0.954           | 0.966             | 0.985               |
| Replicate 2, all peaks            | 227,914           | 0.990          | 0.987           | 0.974             | 0.972               |
| Replicate 2, significant peaks    | 4,690             | 0.963          | 0.954           | 0.979             | 0.985               |
| IDR reproducible peaks            | 3,374             | 0.955          | 0.953           | 0.975             | 0.982               |

The K562 input library is dominated by repeat/rRNA reads (76.1% removed by the repeat filter), leaving only 1.8 M usable reads for input normalisation; the peak numbers and agreement with ENCODE are nonetheless very close.

Overlap is on the same strand, and correlations are of the ENCODE peak with the largest overlap for each workflow peak. Remaining differences are expected: the workflow uses newer versions of the tools (STAR 2.7.11b, cutadapt 5.1, umi_tools/IDR versions, CLIPper from the current YeoLab repository), a Dfam based repeat element reference instead of RepBase, and the random choice of the retained read between PCR duplicates.

## Issues found with this data set

Running real data (ENCSR202BFN) found two problems that the simulated test data did not show. Both are fixed in the workflow, and both are confirmed fixed by the second experiment (ENCSR893RAV), which was run afterwards with no further changes and matched ENCODE at least as well:

1. **Reads from the ENCODE portal are already demultiplexed.** The workflow assumed raw reads and would have tried to demultiplex them. This is now supported with `demultiplexed: True` (each library is processed as a whole, the barcode IDs in `samples.csv` still define the adapters).
2. **Trimming with barcodes that contain random bases (N).** Barcodes such as A03/G07 (used for replicate 2 of ENCSR202BFN) contain 4 random bases. In the second round of adapter trimming (minimum overlap of 5, as in the ENCODE SOP), adapter chunks that start with `NNNN` match the last 5 bases of any read that ends with the next base of the chunk, so 5 real bases were removed from ~58% of the reads. Compared with ENCODE's alignments, only 42% of the reads had the same alignment and 3.8% fewer reads were kept. Second round trimming now only uses the adapter chunks without N (round 1 is unchanged, and barcodes without N such as A01/B06/C01/D8f are not affected): 94% of the reads now have the same alignment as ENCODE and the number of usable reads is within 1.7%. ENCSR893RAV uses barcode combinations without N for both replicates (A01/B06 and C01/D8f), and reaches 99.8% identical alignments, consistent with this explanation.
