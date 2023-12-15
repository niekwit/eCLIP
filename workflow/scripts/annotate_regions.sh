#!/usr/bin/env bash

# prepare bed file for HOMER (chromosome,start,stop,region id, not used, strand)
BED=${snakemake_input[bed]}
BED_HOMER=${snakemake_output[bed]}

echo $BED
echo $BED_HOMER

awk 'BEGIN{OFS=FS="\t"} {print "chr"$1,$2,$3,"region_"NR,".",$6}' $BED > $BED_HOMER 2> ${snakemake_log[0]}

# annotate regions with HOMER and perform GO analysis
ANNOTATED=${snakemake_output[txt]}

annotatePeaks.pl $BED_HOMER ${snakemake_params[genome]} -go ${snakemake_output[go]} > $ANNOTATED 2>> ${snakemake_log[0]}


