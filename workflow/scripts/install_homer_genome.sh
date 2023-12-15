#!/usr/bin/env bash

# install promoter and enhancer annotations for HOMER
GENOME=${snakemake_params[genome]}

if [ $GENOME == "hg19" ] || [ $GENOME == "hg38" ]; then
    perl $CONDA_PREFIX/share/homer/configureHomer.pl -install human-p > ${snakemake_log[0]} 2>&1
elif [ $GENOME == "mm9" ] || [ $GENOME == "mm10" ]; then
    perl $CONDA_PREFIX/share/homer/configureHomer.pl -install mouse-p > ${snakemake_log[0]} 2>&1
else
    echo "Choosen genome (${GENOME}) not (yet) supported by HOMER for peak annotation"
    exit 1
fi

# install remaining genome annotations for HOMER (append stdout to log file)
perl $CONDA_PREFIX/share/homer/configureHomer.pl -install ${snakemake_params[genome]} >> ${snakemake_log[0]} 2>&1

