import subprocess
import pysam
import pandas as pd
import os

threads = str(snakemake.threads)

def get_reads_fastq(fastq):
    """Get total number of reads in fastq file"""
    reads = int(subprocess.check_output(f"zcat {fastq} | wc -l", shell=True).decode().strip()) / 4
    return reads


def pysam_view(bam):
    """Get total number of reads in bam file"""
    reads = int(pysam.view("-c", "-@", threads, bam))
    return reads


def get_reads_from_star_log(log, line):
    """Get number of reads from STAR log"""
    # read entire log file 
    with open(log) as f:
        lines = f.readlines()
        
    # extract relevant line
    line = [x for x in lines if line in x][0]
    
    # extract number of reads
    reads = int(line.split("\t")[1].strip())
    
    return reads

# create df to store results
df = pd.DataFrame(columns = ["sample", "total_reads", "too_short_or_DLE", "too_shorted_mapped", "multimapping_reads", "pcr_duplicates", "usable_reads"])

# for each sample, calculate number of usable reads
for r1,log,ddup in zip(snakemake.input["r1"],
                       snakemake.input["log"],
                       snakemake.input["ddup"],
                              ):
    
    # get sample name
    sample = os.path.basename(r1).replace("_R1_001.fastq.gz", "")
    
    # get total number of reads
    total_reads = get_reads_fastq(r1)
    
    # get number of reads that have too short r1/r2 or have double ligation events
    star_input_reads = get_reads_from_star_log(log, "Number of input reads")
    too_short_or_DLE = total_reads - star_input_reads
        
    # get number of reads that are too short 
    too_short_mapped = get_reads_from_star_log(log, "Number of reads unmapped: too short")
    
    # get number of multi mapping reads    
    multimapping_reads_1 = get_reads_from_star_log(log, "Number of reads mapped to multiple loci")
    multimapping_reads_2 = get_reads_from_star_log(log, "Number of reads mapped to too many loci")
    multimapping_reads = multimapping_reads_1 + multimapping_reads_2
    
    
    # get number of PCR duplicates (number of mapped reads - number of reads after deduplication)
    unique_mapped = get_reads_from_star_log(log, "Uniquely mapped reads number")
    post_ddup_reads = pysam_view(ddup) # usable reads
    pcr_duplicates = unique_mapped - post_ddup_reads
    
    # write results to df
    df.loc[len(df)] = [sample, total_reads, too_short_or_DLE, 
                    too_short_mapped, multimapping_reads, 
                    pcr_duplicates, post_ddup_reads]

# write df to csv (plot bar graph in R)
df.to_csv(snakemake.output[0], index=False)
    
