#!/usr/bin/env bash
set -euo pipefail

exec > "${snakemake_log[0]}" 2>&1

BAM=${snakemake_input[bam]}
CHROM_SIZES=${snakemake_input[cs]}
POS_BW=${snakemake_output[pos]}
NEG_BW=${snakemake_output[neg]}
POS_STRAND=${snakemake_params[pos_strand]}
NEG_STRAND=${snakemake_params[neg_strand]}

# RPM normalisation (as makebigwigfiles.py of the Yeo lab): scale by 1 / (million mapped reads)
N=$(samtools view -c -F 4 "$BAM")
SCALE=$(gawk -v n="$N" 'BEGIN {printf "%.12f", 1000000 / n}')
NEG_SCALE=$(gawk -v n="$N" 'BEGIN {printf "%.12f", -1000000 / n}')

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

genomeCoverageBed -ibam "$BAM" -bg -strand "$POS_STRAND" -scale "$SCALE" -g "$CHROM_SIZES" -du -split |
    LC_ALL=C sort -k1,1 -k2,2n > "$TMP/pos.bg"
bedGraphToBigWig "$TMP/pos.bg" "$CHROM_SIZES" "$POS_BW"

genomeCoverageBed -ibam "$BAM" -bg -strand "$NEG_STRAND" -scale "$NEG_SCALE" -g "$CHROM_SIZES" -du -split |
    LC_ALL=C sort -k1,1 -k2,2n > "$TMP/neg.bg"
bedGraphToBigWig "$TMP/neg.bg" "$CHROM_SIZES" "$NEG_BW"
