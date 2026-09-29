"""
Reproducibility of locus-resolved TE binding enrichment (te_enrichment.py) between two IP
replicates of the same condition: Spearman correlation of log2FC, over loci/families present in
both replicates' enrichment tables (i.e. that passed each replicate's own min_reads filter).

This is used instead of IDR: IDR is designed for genomic-interval peak calling, not enrichment
scores over a fixed set of categorical loci/families.
"""

import sys

import pandas as pd
from scipy.stats import spearmanr

sys.stderr = open(snakemake.log[0], "w")


def reproducibility(path1, path2, id_col):
    s1 = pd.read_csv(path1, sep="\t").set_index(id_col)["log2fc"]
    s2 = pd.read_csv(path2, sep="\t").set_index(id_col)["log2fc"]
    joined = pd.concat([s1, s2], axis=1, join="inner")
    joined.columns = ["log2fc_1", "log2fc_2"]
    if len(joined) < 2:
        return len(joined), float("nan"), float("nan")
    rho, p = spearmanr(joined["log2fc_1"], joined["log2fc_2"])
    return len(joined), rho, p


with open(snakemake.output[0], "w") as f:
    for level, id_col, path1, path2 in [
        ("locus", "locus_id", snakemake.input.locus1, snakemake.input.locus2),
        ("family", "rep_family", snakemake.input.family1, snakemake.input.family2),
    ]:
        n, rho, p = reproducibility(path1, path2, id_col)
        f.write(
            f"{level}-level: n={n} shared {id_col}s, Spearman rho={rho:.3g}, p={p:.3g}\n"
        )
