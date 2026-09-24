# Samples.csv

Use this file to describe all libraries (IP/eCLIP libraries **and** size-matched input (SMInput) libraries), one library per row.

`sample`: sample name that matches the read file name(s) in `reads/` without extension (see below). Only alphanumeric characters and underscores are allowed. IP sample names must end with `_` followed by the replicate number (e.g. `RBFOX2_1`, `RBFOX2_2`). The part before this replicate number is the condition; replicates of the same condition are compared by IDR.

`control`: name of the size-matched input library (which has its own row) that is used for input normalisation of this IP sample. Leave empty for the size-matched input libraries themselves. Replicates can share the same input.

`adapter`: _(single-end only, optional)_ 3' adapter set that was used for the library: `InvRil19` (default, or set by `cutadapt: se_adapter` in `config.yaml`), or `InvRNA1` to `InvRNA8`.

`barcode_a`, `barcode_b`: _(paired-end only, required)_ ID of the inline barcode of the two barcodes that were ligated to the library: `A01`, `B06`, `C01`, `D8f`, `A03`, `G07`, `A04`, `F05` or `NIL` for libraries without barcode (size-matched input).

## Read files

Single-end reads (auto-detected):

```
reads/{sample}.fastq.gz
```

Paired-end reads (auto-detected):

```
reads/{sample}_R1_001.fastq.gz
reads/{sample}_R2_001.fastq.gz
```

## Examples

Single-end:

| sample         | control        | adapter  |
| -------------- | -------------- | -------- |
| RBFOX2_1       | RBFOX2_input_1 | InvRil19 |
| RBFOX2_2       | RBFOX2_input_2 | InvRil19 |
| RBFOX2_input_1 |                | InvRil19 |
| RBFOX2_input_2 |                | InvRil19 |

Paired-end (each IP library carries two inline barcodes, the size-matched input none):

| sample         | control        | barcode_a | barcode_b |
| -------------- | -------------- | --------- | --------- |
| RBFOX2_1       | RBFOX2_input_1 | A01       | B06       |
| RBFOX2_2       | RBFOX2_input_1 | C01       | D8f       |
| RBFOX2_input_1 |                | NIL       | NIL       |

# config.yaml

All settings have the ENCODE eCLIP pipeline (eCLIP-seq Processing Pipeline v2.2) values as default. Use Python style booleans (`True`/`False`).
