import pandas as pd
import sys
import os

class Barcodes:
    '''Class for storing barcode/umi information
    '''
    def __init__(self, config):
        pass


def import_samples():
    try:
        csv = pd.read_csv("config/samples.csv")
        SAMPLES = csv["sample"]
        
        # check if sample names match file names
        for sample in SAMPLES:
            r1= f"reads/{sample}_R1_001.fastq.gz"
            r2= f"reads/{sample}_R2_001.fastq.gz"
            if not os.path.isfile(r1):
                print(f"ERROR: {r1} not found!")
                sys.exit(1)
            if not os.path.isfile(r2):
                print(f"ERROR: {r2} not found!")
                sys.exit(1)
        
        return SAMPLES

    except FileNotFoundError:
        print("ERROR: config/samples.csv not found!")
        sys.exit(1)
        


        