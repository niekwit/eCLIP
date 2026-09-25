#!/usr/bin/env python3
"""
Compares the peaks of this workflow with peaks of an ENCODE eCLIP experiment (see README.md).

Usage:
  python compare_encode_peaks.py --results results --encode encode_peaks \\
      --rep1 ENCFF067JAD --rep2 ENCFF243RQR --idr ENCFF721PWF \\
      --sample1 U2AF2_1 --sample2 U2AF2_2

Requires pandas, numpy and bedtools (in PATH, or use --bedtools).
The ENCODE files (<accession>.bed.gz, narrowPeak) must be in the --encode directory.
"""

import argparse
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--results", default="results", help="results directory of the workflow")
parser.add_argument("--encode", default="encode_peaks", help="directory with the ENCODE peak files")
parser.add_argument("--rep1", required=True, help="ENCODE accession of the peaks of replicate 1")
parser.add_argument("--rep2", required=True, help="ENCODE accession of the peaks of replicate 2")
parser.add_argument("--idr", required=True, help="ENCODE accession of the IDR reproducible peaks")
parser.add_argument("--sample1", required=True, help="name of replicate 1 in the workflow")
parser.add_argument("--sample2", required=True, help="name of replicate 2 in the workflow")
parser.add_argument("--bedtools", default="bedtools")
parser.add_argument("--l10p", type=float, default=3, help="-log10(p-value) cutoff of significant peaks (ENCODE: 3)")
parser.add_argument("--l2fc", type=float, default=3, help="log2 fold change cutoff of significant peaks (ENCODE: 3)")
args = parser.parse_args()

TMP = tempfile.mkdtemp()


def read_encode(accession):
    """ENCODE narrowPeak: signalValue is the log2 fold change, pValue the -log10(p-value)"""
    df = pd.read_csv(os.path.join(args.encode, f"{accession}.bed.gz"), sep="\t", header=None)
    df.columns = ["chrom", "start", "end", "name", "score", "strand", "l2fc", "l10p", "q", "peak"]
    return df[["chrom", "start", "end", "strand", "l10p", "l2fc"]]


def read_mine(path):
    """Workflow BED: chromosome, start, end, -log10(p-value), log2 fold change, strand"""
    df = pd.read_csv(path, sep="\t", header=None)
    df.columns = ["chrom", "start", "end", "l10p", "l2fc", "strand"]
    return df[["chrom", "start", "end", "strand", "l10p", "l2fc"]]


def significant(df):
    return df[(df.l10p >= args.l10p) & (df.l2fc >= args.l2fc)].reset_index(drop=True)


def percentage(n, total):
    """Percentage as string, n/a if there is nothing to compare (empty peak set)"""
    return f"{100 * n / total:.1f}%" if total else "n/a"


def overlaps(a, b, tag):
    """Overlapping peaks (same strand) of a and b, as index pairs with the length of the overlap"""
    for name, df in (("a", a), ("b", b)):
        out = df.assign(i=range(len(df)))[["chrom", "start", "end", "i", "l2fc", "strand"]]
        out.sort_values(["chrom", "start"]).to_csv(f"{TMP}/{tag}_{name}.bed", sep="\t", header=False, index=False)
    out = subprocess.run(
        [args.bedtools, "intersect", "-s", "-wo", "-a", f"{TMP}/{tag}_a.bed", "-b", f"{TMP}/{tag}_b.bed"],
        capture_output=True, text=True, check=True,
    ).stdout
    rows = [line.split("\t") for line in out.strip().split("\n") if line]
    return pd.DataFrame(
        {
            "ia": [int(r[3]) for r in rows],
            "ib": [int(r[9]) for r in rows],
            "overlap": [int(r[12]) for r in rows],  # last column of -wo: number of overlapping bases
        }
    )


def compare(mine, enc, label):
    ov = overlaps(mine, enc, label.replace(" ", "_"))
    summary = {
        "peak set": label,
        "this workflow": len(mine),
        "ENCODE": len(enc),
        "workflow peaks overlapping ENCODE": f"{ov.ia.nunique()} ({percentage(ov.ia.nunique(), len(mine))})",
        "ENCODE peaks overlapping workflow": f"{ov.ib.nunique()} ({percentage(ov.ib.nunique(), len(enc))})",
    }
    # ENCODE peak with the largest overlap for each workflow peak
    best = ov.sort_values("overlap", ascending=False).drop_duplicates("ia")
    m, e = mine.iloc[best.ia.values].reset_index(drop=True), enc.iloc[best.ib.values].reset_index(drop=True)
    corr = {"peak set": label, "overlapping pairs": len(best)}
    for col, name in (("l2fc", "log2FC"), ("l10p", "-log10(p)")):
        if len(best) > 1:
            corr[f"{name} Pearson"] = round(float(np.corrcoef(m[col], e[col])[0, 1]), 3)
            corr[f"{name} Spearman"] = round(float(pd.Series(m[col].values).corr(pd.Series(e[col].values), method="spearman")), 3)
        else:
            corr[f"{name} Pearson"] = corr[f"{name} Spearman"] = "n/a"
    return summary, corr, ov


def not_overlapping(df, ov, index_column):
    return df[~df.index.isin(ov[index_column].unique())]


summaries, correlations = [], []
for label, sample, accession in (("replicate 1", args.sample1, args.rep1), ("replicate 2", args.sample2, args.rep2)):
    enc = read_encode(accession)
    mine = read_mine(os.path.join(args.results, "peaks", sample, f"{sample}.peaks.bed"))
    for kind, m, e in (("all peaks", mine, enc), ("significant peaks", significant(mine), significant(enc))):
        summary, corr, ov = compare(m, e, f"{label}, {kind}")
        summaries.append(summary)
        correlations.append(corr)
        if kind == "all peaks":
            only_enc = not_overlapping(enc, ov, "ib")
            only_mine = not_overlapping(mine, ov, "ia")
            print(
                f"{label}: peaks without overlap: ENCODE {len(only_enc)} (median -log10(p) {only_enc.l10p.median():.2f}, "
                f"{percentage(len(significant(only_enc)), len(only_enc))} significant), workflow {len(only_mine)} "
                f"(median -log10(p) {only_mine.l10p.median():.2f}, {percentage(len(significant(only_mine)), len(only_mine))} significant)"
            )

mine_idr = read_mine(os.path.join(args.results, "idr", f"{args.sample1}_vs_{args.sample2}", f"{args.sample1}_vs_{args.sample2}.reproducible_peaks.bed"))
summary, corr, ov = compare(mine_idr, read_encode(args.idr), "IDR reproducible peaks")
summaries.append(summary)
correlations.append(corr)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
print()
print(pd.DataFrame(summaries).to_string(index=False))
print()
print(pd.DataFrame(correlations).to_string(index=False))
