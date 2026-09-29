"""
Creates an RPM normalised, strand-specific bigWig file from a BAM file, as in the ENCODE eCLIP
pipeline (makebigwigfiles.py of the Yeo lab), reimplemented with bedtools genomecov and UCSC tools.

Reads on `params.strand` are written to `output`, scaled by params.sign * 1e6 / (mapped reads):
the negative strand is written with a negative sign, so that genome browsers show it below zero.
"""

import os
import subprocess
import tempfile
from pathlib import Path

import pysam

log = open(snakemake.log[0], "w")

bam = snakemake.input.bam
chrom_sizes = snakemake.input.cs
output = snakemake.output[0]
strand = snakemake.params.strand
sign = snakemake.params.sign

# Number of mapped reads, for RPM (reads per million) normalisation (same count as
# `samtools view -c -F 4`, taken from the BAM index instead of scanning the whole file)
with pysam.AlignmentFile(bam, "rb") as f:
    mapped_reads = f.mapped
scale = sign * 1_000_000 / mapped_reads

with tempfile.TemporaryDirectory() as tmpdir:
    bedgraph = Path(tmpdir) / "coverage.bg"
    sorted_bedgraph = Path(tmpdir) / "coverage.sorted.bg"

    with open(bedgraph, "w") as f:
        subprocess.run(
            [
                "genomeCoverageBed",
                "-ibam",
                bam,
                "-bg",
                "-strand",
                strand,
                "-scale",
                f"{scale:.12f}",
                "-g",
                chrom_sizes,
                "-du",
                "-split",
            ],
            stdout=f,
            stderr=log,
            check=True,
        )

    with open(sorted_bedgraph, "w") as f:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", bedgraph],
            stdout=f,
            stderr=log,
            env={**os.environ, "LC_ALL": "C"},
            check=True,
        )

    subprocess.run(
        ["bedGraphToBigWig", sorted_bedgraph, chrom_sizes, output],
        stdout=log,
        stderr=log,
        check=True,
    )
