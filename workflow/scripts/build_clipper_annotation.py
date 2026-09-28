"""
Builds a custom CLIPper annotation from a GTF file, for genome/annotation versions that are not
built in to CLIPper (see CLIPper's wiki "Supporting additional species":
https://github.com/YeoLab/clipper/wiki/Supporting-additional-species, and the actual code in
clipper/src/utils.py: build_transcript_data_gtf() and get_exon_bed(), which this replicates).

CLIPper needs, per gene:
  - one representative transcript: the one with the largest genomic span (as CLIPper's own
    build_transcript_data_gtf() picks "the longest gene from a group of transcripts"). Its
    coordinates, spliced length (mrna_length) and genomic span (premrna_length) are written to
    the *.AS.STRUCTURE.COMPILED.gff file that CLIPper reads via --species/--datadir.
  - all its exons, merged with those of every other transcript of the same gene into a
    non-overlapping list, written to regions/*_exons.bed.
"""

import re
import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

ATTR_RE = re.compile(r'(\w+) "([^"]+)"')


def merge_intervals(starts_ends):
    """Merges overlapping/touching (start, end) intervals (0-based, half open), sorted by start"""
    merged = []
    for start, end in sorted(starts_ends):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


exons = []
with open(snakemake.input.gtf) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields[2] != "exon":
            continue
        attrs = dict(ATTR_RE.findall(fields[8]))
        exons.append(
            (
                fields[0],  # chrom
                int(fields[3]) - 1,  # start (0-based)
                int(fields[4]),  # end
                fields[6],  # strand
                attrs["gene_id"],
                attrs["transcript_id"],
            )
        )

exons = pd.DataFrame(exons, columns=["chrom", "start", "end", "strand", "gene_id", "transcript_id"])
exons["length"] = exons["end"] - exons["start"]
print(f"{len(exons)} exons, {exons.transcript_id.nunique()} transcripts, {exons.gene_id.nunique()} genes")

# One representative transcript per gene: the one with the largest genomic span, as CLIPper's own
# build_transcript_data_gtf() picks (not necessarily the one with the most spliced sequence)
transcripts = exons.groupby("transcript_id").agg(
    chrom=("chrom", "first"),
    strand=("strand", "first"),
    gene_id=("gene_id", "first"),
    start=("start", "min"),
    end=("end", "max"),
    mrna_length=("length", "sum"),
)
transcripts["span"] = transcripts["end"] - transcripts["start"]
representative = transcripts.loc[transcripts.groupby("gene_id")["span"].idxmax()]

# *.AS.STRUCTURE.COMPILED.gff: one gene per line, 1-based inclusive coordinates (GFF convention)
with open(snakemake.output.gff, "w") as f:
    for gene_id, t in representative.set_index("gene_id").iterrows():
        attrs = f"gene_id={gene_id};mrna_length={t.mrna_length};premrna_length={t.span}"
        f.write(
            f"{t.chrom}\tAS_STRUCTURE\tgene\t{t.start + 1}\t{t.end}\t.\t{t.strand}\t.\t{attrs}\n"
        )
print(f"wrote {len(representative)} genes to {snakemake.output.gff}")

# regions/*_exons.bed: exons of every transcript of a gene, merged into a non-overlapping list
with open(snakemake.output.exons, "w") as f:
    n = 0
    for (gene_id, strand), group in exons.groupby(["gene_id", "strand"]):
        for start, end in merge_intervals(list(zip(group.start, group.end))):
            f.write(f"{group.chrom.iloc[0]}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")
            n += 1
print(f"wrote {n} merged exons to {snakemake.output.exons}")
