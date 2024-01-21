import pandas as pd
#import csv
import sys
import os


def import_samples():
    csv = pd.read_csv("config/samples.csv")
    SAMPLES = csv["sample"]
    
    # check if sample names match file names
    counter = 0
    not_found = []
    for sample in SAMPLES:
        r1= f"reads/{sample}_R1_001.fastq.gz"
        r2= f"reads/{sample}_R2_001.fastq.gz"
        if not os.path.isfile(r1):
            not_found.append(r1)
            counter += 1
        if not os.path.isfile(r2):
            not_found.append(r2)
            counter += 1
    if counter != 0:
        not_found = "\n".join(not_found)
        raise ValueError(f"ERROR: {counter} files not found:\n{not_found}")
        
    return SAMPLES
       

def match_samples_with_input():
    """Match each sample to its input file for input normalization step"""
    # read samples.csv
    csv = pd.read_csv("config/samples.csv")
    
    # get non-input samples (i.e. ip samples)
    samples = csv["sample"].to_list()
    ip_samples = [sample for sample in samples if not "input" in sample]
    
    # match each ip sample with its input sample
    factors = csv["factor"].unique().tolist() # antibodies
    factors = [x for x in factors if not "input" in x] # remove input
    treatments = csv["treatment"].unique().tolist() # treatments
    
    if len(treatments) == 1:
        input_samples = [f"{s.split('_')[0]}_input" for s in ip_samples]
    else:
        csv_ip = csv[csv["sample"].isin(ip_samples)]
        genotypes = csv_ip["genotype"].tolist()
        treatments = csv_ip["treatment"].tolist()
        input_samples = [f"{g}_input_{t}" for g,t in zip(genotypes,treatments)]
    '''
    # check if these input files exist
    counter = 0
    for sample in input_samples:
        R2_bam = f"results/mapped/{sample}/{sample}_R2.bam"
        if not os.path.isfile(R2_bam):
            print(f"ERROR: {R2_bam} not found!")
            counter += 1
        if counter != 0:
            print(f"ERROR: {counter} input files not found!")
            print("Make sure to follow the naming convention for all files...")
            raise ValueError(f"ERROR: {counter} input files not found!")  
    '''
    return ip_samples, input_samples
    

        