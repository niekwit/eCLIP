import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def load_bed(bed, out):
    """Load bed file as pandas df"""
    try:
        # load bed file as pandas df and add header names
        bed = pd.read_csv(bed, sep="\t", header=None)
        bed.columns = ["chr", "start", "end", "name", "score", "strand", "score"]
        
        # name column is actually Pureclip state, change to unique crosslink names
        bed["name"] = [f"crosslink_{i}" for i in range(1, len(bed) + 1)]
        
        return bed
    
    except pd.errors.EmptyDataError: # if bed file is empty (no peaks)
        print(f"Bed file ({bed}) is empty")
        print("Writing empty fasta file")
        with open(out, "w") as f:
            f.write("")
        sys.exit(0)


def load_fasta(fasta):
    """Load fasta file and return dict with chr as key and sequence as value"""
    chr_seq = {}
    for chr_ in SeqIO.parse(fasta, "fasta"):
        chr_seq[chr_.id] = chr_.seq
    return chr_seq


def write_dict2fasta(d, out):
    """Write dict to fasta file"""
    with open(out, "w") as f:
        for name, seq in d.items():
            f.write(f">{name}\n{seq}\n")


def main(bed, chr_seq):
    # dict to store sequence per region
    region_seq = {}

    # for each region in bed file, extract sequence and write to fasta file
    for row in bed.itertuples(index=False):
        chr = row.chr.replace("chr", "") # convert back to Ensembl format
        start = row.start - int(snakemake.params["flank"])
        end = row.end + int(snakemake.params["flank"])
        name = row.name
        strand = row.strand

        # extract sequence
        seq = chr_seq[chr][start:end]
        
        # reverse complement if on negative strand
        if strand == "-":
            seq = seq.reverse_complement()
        
        # store sequence in dict
        region_seq[f"{name}_{chr}_{start}_{strand}"] = seq

    # write to fasta file
    write_dict2fasta(region_seq, snakemake.output["out"])
           
if __name__ == "__main__":
    bed = load_bed(snakemake.input["bed"], snakemake.output["out"])
    chr_seq = load_fasta(snakemake.input["fasta"])
    
    main(bed, chr_seq)
    
    