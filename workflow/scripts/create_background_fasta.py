import random
import convert_bed2fasta as utils
from Bio.Seq import Seq

bed = utils.load_bed(snakemake.input["bed"], snakemake.output["out"]) # bed file with crosslink position
fasta = utils.load_fasta(snakemake.input["fasta"]) # whole genome fasta
flank = int(snakemake.params["flank"])

# number of regions in bed file
n = len(bed)

# get chromosome names from fasta dict
chr_names = list(fasta.keys())

# get value lengths of each key in fasta dict
chr_lengths = [len(v) for v in fasta.values()]

# create dict with chromosome names as keys and chromosome lengths as values
chr_seq = dict(zip(chr_names, chr_lengths))

# dict to store random regions
random_regions = {}


def random_seq():
    """Returns random genomic sequence with no Ns in sequence"""
    while True:
        # randomly select chromosome
        chr = random.choice(chr_names)
        
        # randomly select strand
        strand = random.choice(["+", "-"])
        
        # randomly select start position
        start = random.randint(0, chr_seq[chr] - (2 * flank + 1)) # +1 to account for crosslink
        
        # calculate end position
        end = start + (2 * flank + 1)
        
        # extract sequence
        seq = fasta[chr][start:end]
        
        # reverse complement if on negative strand
        if strand == "-":
            seq = seq.reverse_complement()
        
        # check if sequence contains any Ns, if so pick another sequence
        if "N" not in seq:
            return seq, chr, start, strand

# create n random regions
for i in range(1, n + 1):
    # get random sequence, chromosome, start position and strand
    seq,chr,start,strand = random_seq()
    
    # store sequence in dict
    random_regions[f"random_{i}_{chr}_{start}_{strand}"] = seq

# write random regions to fasta file
utils.write_dict2fasta(random_regions, snakemake.output["out"])

