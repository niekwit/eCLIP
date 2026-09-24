#!/usr/bin/env python3
"""
Creates the test data of this directory (see README.md):

  resources/GRCh38.primary_assembly.genome.fa   chr21, real sequence only around the test genes (N elsewhere)
  resources/gencode.v29.annotation.gtf          GENCODE v29 annotation of the test genes
  resources/ENCFF269URO.bed                     ENCODE eCLIP blacklist (hg38)
  reads/*.fastq.gz                              simulated single-end eCLIP reads
  clipper_genes.txt                             --gene arguments for CLIPper (clipper: extra in config/config.yaml)

Usage: python make_test_data.py [download_dir]   (only Python 3 standard library is needed)
"""

import gzip
import json
import os
import random
import re
import sys
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DOWNLOADS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, ".downloads")

CHR21_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz"
GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz"
BLACKLIST_URL = "https://www.encodeproject.org/files/ENCFF269URO/@@download/ENCFF269URO.bed.gz"
ALU_URL = "https://www.dfam.org/api/families/DF000000002/sequence?format=fasta"  # AluY consensus

SEED = 7
MIN_TRANSCRIPT_LENGTH = 400  # nt of exons of the longest transcript of a gene
MAX_GENE_END = 20_000_000  # only genes at the start of chr21 (keeps the genome file small)
FLANK = 2000  # bp of real sequence on each side of a gene
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATATCTCGTATGCCGTCTTCTGCTTG"
READ_LENGTH = 75
UMI_LENGTH = 10

random.seed(SEED)
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def download(url, name):
    os.makedirs(DOWNLOADS, exist_ok=True)
    path = os.path.join(DOWNLOADS, name)
    if not os.path.isfile(path):
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, path)
    return path


# --- genome and annotation ---
genome = "".join(
    line.strip().upper()
    for line in gzip.open(download(CHR21_URL, "chr21.fa.gz"), "rt")
    if not line.startswith(">")
)
gtf = list(gzip.open(download(GTF_URL, "gencode.v29.annotation.gtf.gz"), "rt"))


def attributes(field):
    return dict(re.findall(r'(\w+) "([^"]+)"', field))


# longest transcript of each protein coding/lncRNA gene on chr21: (length, strand, exons)
exons = {}
for line in gtf:
    if line.startswith("#"):
        continue
    f = line.rstrip().split("\t")
    if f[0] != "chr21" or f[2] != "exon":
        continue
    a = attributes(f[8])
    if a.get("gene_type") not in ("protein_coding", "lncRNA"):
        continue
    exons.setdefault(a["transcript_id"], (a["gene_id"], f[6], []))[2].append(
        (int(f[3]) - 1, int(f[4]))
    )

longest = {}
for gene_id, strand, ex in exons.values():
    length = sum(end - start for start, end in ex)
    if gene_id not in longest or length > longest[gene_id][0]:
        longest[gene_id] = (length, strand, sorted(ex))

genes = {
    g: v
    for g, v in longest.items()
    if v[0] >= MIN_TRANSCRIPT_LENGTH and max(e for s, e in v[2]) < MAX_GENE_END
}
print(f"{len(genes)} test genes")

# gene loci
loci = []
for line in gtf:
    if line.startswith("#"):
        continue
    f = line.rstrip().split("\t")
    if f[0] == "chr21" and f[2] == "gene" and attributes(f[8])["gene_id"] in genes:
        loci.append((int(f[3]) - 1, int(f[4])))
genome_end = max(end for start, end in loci) + FLANK

# Genome: real sequence around the genes, N elsewhere (coordinates are unchanged, which is
# required for the annotation that is built in to CLIPper)
masked = ["N"] * genome_end
for start, end in loci:
    for i in range(max(0, start - FLANK), min(genome_end, end + FLANK)):
        masked[i] = genome[i]
masked = "".join(masked)

os.makedirs(os.path.join(HERE, "resources"), exist_ok=True)
os.makedirs(os.path.join(HERE, "reads"), exist_ok=True)

with open(os.path.join(HERE, "resources", "GRCh38.primary_assembly.genome.fa"), "w") as f:
    f.write(">chr21\n")
    for i in range(0, len(masked), 60):
        f.write(masked[i : i + 60] + "\n")

with open(os.path.join(HERE, "resources", "gencode.v29.annotation.gtf"), "w") as f:
    for line in gtf:
        if line.startswith("#"):
            f.write(line)
        else:
            fields = line.split("\t")
            if fields[0] == "chr21" and attributes(fields[8])["gene_id"] in genes:
                f.write(line)

with gzip.open(download(BLACKLIST_URL, "ENCFF269URO.bed.gz"), "rt") as f_in:
    with open(os.path.join(HERE, "resources", "ENCFF269URO.bed"), "w") as f_out:
        f_out.write(f_in.read())

with open(os.path.join(HERE, "clipper_genes.txt"), "w") as f:
    f.write(" ".join(f"--gene {g.split('.')[0]}" for g in sorted(genes)) + "\n")

# --- reads ---
transcripts = []
for gene_id, (length, strand, ex) in genes.items():
    seq = "".join(genome[start:end] for start, end in ex)
    if strand == "-":
        seq = revcomp(seq)
    if "N" not in seq:
        transcripts.append((gene_id, seq))

alu_response = urllib.request.urlopen(ALU_URL, timeout=60).read().decode()
try:
    alu_response = json.loads(alu_response)["body"]
except ValueError:
    pass
alu = "".join(x for x in alu_response.split("\n") if not x.startswith(">")).upper()


def umi():
    return "".join(random.choice("ACGT") for _ in range(UMI_LENGTH))


def make_read(insert):
    """UMI + RNA fragment (sense) + 3' adapter"""
    return (umi() + insert + ADAPTER)[:READ_LENGTH]


def simulate(genes_with_sites, n_background, n_repeat, reads_per_site=(25, 70)):
    reads = []
    # binding sites: reads that start around the site (5' end pile up), 30% PCR duplicates (same UMI)
    for gene_id, seq, sites in genes_with_sites:
        for center, strength in sites:
            for _ in range(int(random.randint(*reads_per_site) * strength)):
                start = max(0, min(int(random.gauss(center, 6)), len(seq) - 30))
                read = make_read(seq[start : start + random.randint(25, 55)])
                reads.append(read)
                if random.random() < 0.3:
                    reads.extend([read] * random.randint(1, 3))
    # background: uniform over the transcripts
    for _ in range(n_background):
        gene_id, seq, sites = random.choice(genes_with_sites)
        start = random.randint(0, len(seq) - 40)
        reads.append(make_read(seq[start : start + random.randint(25, 55)]))
    # Alu repeat reads (removed by the repeat element filter)
    for _ in range(n_repeat):
        start = random.randint(0, len(alu) - 40)
        reads.append(make_read(alu[start : start + random.randint(25, 55)]))
    random.shuffle(reads)
    return reads


def write_fastq(name, reads):
    with gzip.open(os.path.join(HERE, "reads", f"{name}.fastq.gz"), "wt") as f:
        for i, read in enumerate(reads):
            f.write(f"@{name}.{i} 1:N:0:ATCACG\n{read}\n+\n{'I' * len(read)}\n")


# Binding sites shared by the replicates (3-6 per gene)
shared = []
for gene_id, seq in transcripts:
    sites = [
        (random.randint(80, max(81, len(seq) - 80)), random.choice([1.0, 1.0, 0.7]))
        for _ in range(random.randint(3, 6))
    ]
    shared.append((gene_id, seq, sites))

# IP replicates: 90% of the shared sites, plus a replicate specific site in 20% of the genes
for replicate in (1, 2):
    ip = []
    for gene_id, seq, sites in shared:
        sites = [x for x in sites if random.random() < 0.9]
        if random.random() < 0.2:
            sites.append((random.randint(80, max(81, len(seq) - 80)), 0.8))
        ip.append((gene_id, seq, sites))
    write_fastq(f"RBFOX2_{replicate}", simulate(ip, n_background=3000, n_repeat=1500))

# Size-matched input: background only (shared by both replicates)
input_ = [(gene_id, seq, []) for gene_id, seq, sites in shared]
write_fastq("RBFOX2_input_1", simulate(input_, n_background=9000, n_repeat=3000))
print("Done")
