"""
Downloads the UCSC RepeatMasker table (rmsk.txt.gz) and writes a BED of individual TE copies
(genomic coordinates + family), for the locus-resolved TE binding analysis (workflow/rules/
te_repeats.smk). This is the individual-copy counterpart of the family-level Dfam consensus
reference (get_repeat_elements.py) used for the repeat element pre-filter step: rmsk has no
consensus sequences, but does have per-copy genomic coordinates, which Dfam's API does not.

Only genuine transposable element classes are kept (LINE/SINE/LTR/DNA/RC/Retroposon), excluding
non-TE repeat classes the same table also lists (rRNA, tRNA, snRNA, scRNA, Satellite, Simple_repeat,
Low_complexity, Unknown, ...).
"""

import gzip
import sys
import tempfile
import urllib.request

TE_CLASSES = {"LINE", "SINE", "LTR", "DNA", "RC", "Retroposon"}

# rmsk.txt.gz columns (UCSC schema): bin, swScore, milliDiv, milliDel, milliIns, genoName,
# genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id
GENO_NAME, GENO_START, GENO_END, STRAND, REP_NAME, REP_CLASS, REP_FAMILY, MILLI_DIV = (
    5,
    6,
    7,
    9,
    10,
    11,
    12,
    2,
)

if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    n_total = 0
    n_kept = 0
    with tempfile.NamedTemporaryFile(suffix=".txt.gz") as raw:
        urllib.request.urlretrieve(snakemake.params.url, raw.name)

        with gzip.open(raw.name, "rt") as f_in, open(snakemake.output[0], "w") as f_out:
            for line in f_in:
                n_total += 1
                fields = line.rstrip("\n").split("\t")
                rep_class = fields[REP_CLASS].split("/")[
                    0
                ]  # e.g. "DNA/hAT-Charlie" -> "DNA"
                if rep_class not in TE_CLASSES:
                    continue
                n_kept += 1
                strand = "-" if fields[STRAND] == "C" else fields[STRAND]
                f_out.write(
                    "\t".join(
                        [
                            fields[GENO_NAME],
                            fields[GENO_START],
                            fields[GENO_END],
                            fields[REP_NAME],
                            fields[REP_FAMILY],
                            strand,
                            rep_class,
                            fields[MILLI_DIV],
                        ]
                    )
                    + "\n"
                )

    print(f"{n_kept} / {n_total} RepeatMasker entries are TE loci ({TE_CLASSES})")
