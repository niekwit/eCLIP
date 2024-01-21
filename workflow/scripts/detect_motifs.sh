#!/usr/bin/env bash

# keep script going even if a command within it fails
set +e

LOG=${snakemake_log[0]}
FASTA=${snakemake_input[fasta]}

findMotifs.pl \
$FASTA \
fasta \
${snakemake_params[dr]} \
-fasta ${snakemake_input[background]} \
-rna \
-p ${snakemake[threads]} > $LOG 2>&1

EXITCODE=$?

# if no peaks were in bed file, findMotifs.pl exit code is still 0 !
# check log file for no data error message
# if no data error message is present, create empty output file and exit with 0

if [ $EXITCODE -eq 0 ]
then
    # check if findMotifs.pl failed due to no data in input file
    CHECK=$(grep "There is no data in your input file" $LOG | wc -l)

    if [ $CHECK -eq 1 ]
    then
        echo "WARNING: findMotifs.pl for $FASTA failed, most likely due to no regions in input file..." | tee 2>> $LOG
        echo "Please check log to confirm..." | tee 2>> $LOG
        echo "Creating empty output file..." | tee 2>> $LOG
        mkdir -p ${snakemake_params[dr]}
        touch ${snakemake_output[html]}
        exit 0
    else
        exit 0
    fi
else 
    echo "ERROR: findMotifs.pl failed with exit code $EXITCODE..." | tee 2>> $LOG
    echo "Please check log..." | tee 2>> $LOG
    exit $EXITCODE
fi

