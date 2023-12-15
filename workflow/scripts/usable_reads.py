import pysam
import pandas as pd



#total_reads = number of reads in unprocessed fastq file
#too_short = difference between total_reads and number of reads after all trimming steps
#multimapping_reads = number of reads that map to multiple locations (after alignment)

#pcr_duplicates = difference between total_reads and number of usable_reads
#usable_reads = number of reads after deduplication and filtering
