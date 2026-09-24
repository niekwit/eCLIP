"""
Converts BED file with input normalised peaks to narrowPeak file and to a BED file that can be
converted to bigBed (as in ENCODE eCLIP pipeline).

Input BED: chrom, start, end, -log10(p-value), log2 fold change, strand
"""

import os
import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

COLUMNS = ["chrom", "start", "end", "pvalue", "l2fc", "strand"]

peaks = pd.read_csv(snakemake.input[0], sep="\t", names=COLUMNS)
db = snakemake.params.db
l10p = snakemake.params.l10p
l2fc = snakemake.params.l2fc

# narrowPeak: score is 1000 for significant peaks and 200 for others (only for coloring in genome browser)
name = os.path.basename(snakemake.input[0])
header = (
    f'track type=narrowPeak visibility=3 db={db} name="{name}" '
    f'description="{name} input-normalized peaks"'
)
narrowpeak = pd.DataFrame(
    {
        "chrom": peaks["chrom"],
        "start": peaks["start"],
        "end": peaks["end"],
        "name": ".",
        "score": ((peaks["pvalue"] >= l10p) & (peaks["l2fc"] >= l2fc)).map(
            {True: 1000, False: 200}
        ),
        "strand": peaks["strand"],
        "signalValue": peaks["l2fc"],
        "pValue": peaks["pvalue"],
        "qValue": -1,
        "peak": -1,
    }
)
with open(snakemake.output.narrowpeak, "w") as f:
    f.write(f"{header}\n")
narrowpeak.to_csv(
    snakemake.output.narrowpeak, sep="\t", header=False, index=False, mode="a"
)

# bigBed: p-value and log2 fold change are stored in name
fixed = pd.DataFrame(
    {
        "chrom": peaks["chrom"],
        "start": peaks["start"],
        "end": peaks["end"],
        "name": peaks["pvalue"].astype(str) + "|" + peaks["l2fc"].astype(str),
        "score": 0,
        "strand": peaks["strand"],
    }
)
fixed.to_csv(snakemake.output.fixed, sep="\t", header=False, index=False)
