#!/usr/bin/env bash
set -euo pipefail

exec > "${snakemake_log[0]}" 2>&1

FASTA=${snakemake_input[fasta]}
FAI=${snakemake_input[fai]}
GTF=${snakemake_params[gtf]}
GENOME_DIR=${snakemake_output[0]}
THREADS=${snakemake[threads]}

# Settings as recommended in the STAR manual:
# --genomeSAindexNbases: min(14, log2(GenomeLength)/2 - 1)
# --genomeChrBinNbits: min(18, log2(max(GenomeLength/NumberOfReferences, ReadLength)))
# The latter is important for the repeat element reference, which has a lot of small sequences
read -r SA_INDEX_NBASES CHR_BIN_NBITS < <(gawk '
    {len += $2; n += 1}
    END {
        sa = log(len) / log(2) / 2 - 1
        if (sa > 14) sa = 14
        avg = len / n
        if (avg < 100) avg = 100
        bits = log(avg) / log(2)
        if (bits > 18) bits = 18
        printf "%d %d\n", sa, bits
    }' "$FAI")

mkdir -p "$GENOME_DIR"

ARGS=(
    --runMode genomeGenerate
    --genomeDir "$GENOME_DIR"
    --genomeFastaFiles "$FASTA"
    --runThreadN "$THREADS"
    --genomeSAindexNbases "$SA_INDEX_NBASES"
    --genomeChrBinNbits "$CHR_BIN_NBITS"
    --outFileNamePrefix "$GENOME_DIR/"
    --outTmpDir "$GENOME_DIR/_STARtmp"
)

if [ -n "$GTF" ]; then
    ARGS+=(--sjdbGTFfile "$GTF")
fi

STAR "${ARGS[@]}"
