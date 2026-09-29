"""
Assigns each uniquely-mapped, deduplicated TE-candidate read (see workflow/rules/te_repeats.smk)
to an individual TE locus (RepeatMasker copy). A read that is fully internal to a repeat copy
multi-maps genome-wide and was already discarded by the unique-mapping re-alignment step earlier
in this chain; only reads with a unique anchor (e.g. a TE-to-flank readthrough junction) survive
to reach this step, so overlapping a RepeatMasker locus here is real locus-level evidence.

A read overlapping more than one RepeatMasker entry (nested/adjacent elements) is assigned to the
one with the largest overlap (same largest-overlap resolution as validation/compare_encode_peaks.py
uses for peak-to-ENCODE-peak matching).

If `require_family_match` is set, the assignment is dropped unless the RepeatMasker family at the
locus agrees with the family the read hit during the repeat element pre-filter step (Dfam consensus
family, read from the RNAME of results/star/repeats/{unit}.Aligned.out.bam via the per-unit
candidate_family lookup) -- both vocabularies mostly overlap but are not guaranteed identical, so
this is a best-effort cross-check, not a strict identity.
"""

import subprocess
import sys
import tempfile

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

COLUMNS = [
    "a_chrom", "a_start", "a_end", "qname", "mapq", "a_strand",
    "chrom", "start", "end", "rep_name", "rep_family", "strand", "rep_class", "milli_div",
    "overlap",
]  # fmt: skip

with tempfile.NamedTemporaryFile(suffix=".tsv") as tmp:
    with open(tmp.name, "w") as f:
        subprocess.run(
            [
                "bedtools",
                "intersect",
                "-abam",
                snakemake.input.bam,
                "-b",
                snakemake.input.loci,
                "-bed",
                "-wo",
            ],  # fmt: skip
            stdout=f,
            stderr=sys.stderr,
            check=True,
        )
    hits = pd.read_csv(tmp.name, sep="\t", header=None, names=COLUMNS)

print(f"{len(hits)} read-locus overlaps ({hits['qname'].nunique()} unique reads)")

# Multi-overlapping reads: keep the largest-overlap assignment only
hits = hits.loc[hits.groupby("qname")["overlap"].idxmax()]

if snakemake.params.require_family_match:
    lookup = pd.concat(
        pd.read_csv(path, sep="\t", header=None, names=["qname", "family"])
        for path in snakemake.input.family_lookup
    )
    family = dict(zip(lookup["qname"], lookup["family"]))
    hits["consensus_family"] = hits["qname"].map(family)
    matched = hits["consensus_family"].str.casefold() == hits["rep_name"].str.casefold()
    print(f"{matched.sum()} / {len(hits)} reads pass the family-match cross-check")
    hits = hits[matched]

hits["locus_id"] = (
    hits["chrom"] + ":" + hits["start"].astype(str) + "-" + hits["end"].astype(str)
    + ":" + hits["rep_name"]
)  # fmt: skip

hits[
    [
        "qname",
        "locus_id",
        "chrom",
        "start",
        "end",
        "rep_name",
        "rep_family",
        "rep_class",
    ]
].to_csv(snakemake.output[0], sep="\t", index=False)
