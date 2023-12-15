#!/usr/bin/env bash

# fix issue with not finding winextract (line 63 in compute_CLmotif_scores.sh)
#SCRIPT_PATH=$(which compute_CLmotif_scores.sh)
#sed -i 's,if \[ ! -e "$WINEXTRACT" \]; then,if ! winextract --version \>/dev/null 2\>\&1; then,' $SCRIPT_PATH

export BEDTOOLS=$(which bedtools)        # if not specified, PATH is searched
export FIMO=$(which fimo)              # if not specified, PATH is searched
export WINEXTRACT=$(which winextract)    # built together with PureCLIP

# run compute_CLmotif_scores.sh
compute_CLmotif_scores.sh ${snakemake_input[fasta]} ${snakemake_input[bam]} ${snakemake_input[xml]} ${snakemake_input[txt]} ${snakemake_output[0]} > ${snakemake_log[0]} 2>&1


