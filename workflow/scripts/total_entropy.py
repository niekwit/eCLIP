"""
Calculates total entropy (relative information content) of the significant peaks of a sample,
as in the ENCODE eCLIP pipeline (calculate_entropy.py)
"""

import sys
import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

FULL_COLUMNS = [
    "chrom",
    "start",
    "end",
    "peak",
    "ip_count",
    "input_count",
    "pvalue",
    "chivalue",
    "chitype",
    "isenriched",
    "l10p",
    "l2fc",
]

with open(snakemake.input.ip_num) as f:
    ip_mapped = int(f.readline().strip())
with open(snakemake.input.input_num) as f:
    input_mapped = int(f.readline().strip())

peaks = pd.read_csv(snakemake.input.full, sep="\t", names=FULL_COLUMNS)
peaks = peaks[
    (peaks["l10p"] >= snakemake.params.l10p) & (peaks["l2fc"] >= snakemake.params.l2fc)
]

p_ip = peaks["ip_count"] / ip_mapped
p_input = peaks["input_count"] / input_mapped
total_entropy = float((p_ip * np.log2(p_ip / p_input)).sum())

with open(snakemake.output[0], "w") as f:
    f.write(f"{total_entropy}\n")
