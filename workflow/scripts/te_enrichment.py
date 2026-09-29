"""
IP vs size-matched input enrichment of locus-resolved TE binding (see te_locus_assignment.py):
Fisher's exact test per TE locus/family, BH-corrected, normalised against the total pool of
confidently locus-assigned reads of each library (not whole-library read counts) -- this isolates
family/locus-specific enrichment from any bulk shift in overall repeat content between IP and
input. The per-locus table is written as-is; the per-family table is the same computation on the
locus table rolled up (grouped) by rep_family, giving a coarser, less sparse summary view.
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

sys.stderr = open(snakemake.log[0], "w")

ip = pd.read_csv(snakemake.input.ip_reads, sep="\t")
inp = pd.read_csv(snakemake.input.input_reads, sep="\t")
min_reads = snakemake.params.min_reads
ip_total, input_total = len(ip), len(inp)

print(f"IP: {ip_total} locus-assigned reads, input: {input_total} locus-assigned reads")


def enrich(ip, inp, group_col, meta_cols):
    counts = pd.concat(
        [
            ip.groupby(group_col).size().rename("ip_count"),
            inp.groupby(group_col).size().rename("input_count"),
        ],
        axis=1,
    ).fillna(0).astype(int)

    if meta_cols:
        meta = pd.concat([ip, inp]).drop_duplicates(group_col).set_index(group_col)[meta_cols]
        counts = counts.join(meta)

    counts = counts[counts["ip_count"] + counts["input_count"] >= min_reads].copy()

    # Symmetric pseudocount (+1) for RPM/log2FC to avoid log(0); Fisher's exact uses raw counts
    counts["ip_rpm"] = (counts["ip_count"] + 1) / ip_total * 1e6
    counts["input_rpm"] = (counts["input_count"] + 1) / input_total * 1e6
    counts["log2fc"] = np.log2(counts["ip_rpm"] / counts["input_rpm"])
    counts["pvalue"] = [
        fisher_exact([[a, ip_total - a], [c, input_total - c]])[1]
        for a, c in zip(counts["ip_count"], counts["input_count"])
    ]
    counts["qvalue"] = (
        multipletests(counts["pvalue"], method="fdr_bh")[1] if len(counts) else []
    )

    return counts.drop(columns=["ip_rpm", "input_rpm"]).sort_values("qvalue")


locus = enrich(ip, inp, "locus_id", ["chrom", "start", "end", "rep_name", "rep_family", "rep_class"])
locus.to_csv(snakemake.output.locus, sep="\t")

family = enrich(ip, inp, "rep_family", [])
family.to_csv(snakemake.output.family, sep="\t")
