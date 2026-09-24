# ENCODE reproducible peaks (Yeo lab merge_peaks pipeline): replicates of the same condition are compared using the
# Irreproducible Discovery Rate (IDR) framework, with peaks ranked by entropy (information content) over the input
rule compress_peaks_full:
    input:
        full="results/peaks/{sample}/{sample}.normed.bed.full",
        script=yeolab_script(
            "compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl"
        ),
    output:
        bed="results/idr/entropy/{sample}.compressed.bed",
        full="results/idr/entropy/{sample}.compressed.bed.full",
    log:
        "logs/idr/compress_peaks/{sample}.log",
    conda:
        "../envs/peaks.yaml"
    threads: 1
    resources:
        runtime=config["resources"]["idr"]["time"],
    shell:
        "perl {input.script} {input.full} {output.bed} {output.full} > {log} 2>&1"


rule information_content:
    input:
        unpack(peaks_input),
        full="results/idr/entropy/{sample}.compressed.bed.full",
        script=yeolab_script("make_informationcontent_from_peaks.pl"),
    output:
        entropy="results/idr/entropy/{sample}.entropy.full",
        excess_reads="results/idr/entropy/{sample}.entropy.excess_reads",
    log:
        "logs/idr/information_content/{sample}.log",
    conda:
        "../envs/peaks.yaml"
    threads: 1
    resources:
        runtime=config["resources"]["idr"]["time"],
        mem_mb=16000,
    shell:
        "perl {input.script} "
        "{input.full} "
        "{input.ip_num} "
        "{input.input_num} "
        "{output.entropy} "
        "{output.excess_reads} "
        "> {log} 2>&1"


rule entropy_bed:
    input:
        full="results/idr/entropy/{sample}.entropy.full",
        script=yeolab_script("full_to_bed.py"),
    output:
        "results/idr/entropy/{sample}.entropy.bed",
    log:
        "logs/idr/entropy_bed/{sample}.log",
    conda:
        "../envs/peaks.yaml"
    threads: 1
    resources:
        runtime=15,
    shell:
        # Columns 7-10 are added because IDR >= 2.0.3 expects narrowPeak-like BED files
        # (ranking is done on column 5: entropy)
        "python {input.script} --input {input.full} --output {output}.tmp > {log} 2>&1 && "
        "awk 'BEGIN {{OFS=\"\\t\"}} {{print $0, 0, 0, 0, -1}}' {output}.tmp > {output} 2>> {log} && "
        "rm {output}.tmp"


rule idr:
    input:
        rep1="results/idr/entropy/{s1}.entropy.bed",
        rep2="results/idr/entropy/{s2}.entropy.bed",
    output:
        idr="results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.idr.out",
        plot="results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.idr.out.png",
    log:
        "logs/idr/idr/{s1}_vs_{s2}.log",
    conda:
        "../envs/idr.yaml"
    threads: 1
    resources:
        runtime=config["resources"]["idr"]["time"],
    shell:
        "idr "
        "--samples {input.rep1} {input.rep2} "
        "--input-file-type bed "
        "--rank 5 "
        "--peak-merge-method max "
        "--plot "
        "-o {output.idr} "
        "> {log} 2>&1"


rule parse_idr_peaks:
    input:
        idr="results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.idr.out",
        entropy1="results/idr/entropy/{s1}.entropy.full",
        entropy2="results/idr/entropy/{s2}.entropy.full",
        script=yeolab_script("parse_idr_peaks.pl"),
    output:
        "results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.idr.out.bed",
    log:
        "logs/idr/parse_idr_peaks/{s1}_vs_{s2}.log",
    conda:
        "../envs/peaks.yaml"
    threads: 1
    resources:
        runtime=config["resources"]["idr"]["time"],
    shell:
        "perl {input.script} "
        "{input.idr} "
        "{input.entropy1} "
        "{input.entropy2} "
        "{output} "
        "> {log} 2>&1"


rule input_normalisation_idr:
    # Normalises the reads of each replicate over the input in the IDR peak regions
    input:
        unpack(peaks_input),
        peaks="results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.idr.out.bed",
        script=yeolab_script("overlap_peakfi_with_bam.pl"),
    output:
        bed="results/idr/{s1}_vs_{s2}/{sample}.idr_peaks.normed.bed",
        full="results/idr/{s1}_vs_{s2}/{sample}.idr_peaks.normed.bed.full",
    log:
        "logs/idr/input_normalisation/{s1}_vs_{s2}/{sample}.log",
    conda:
        "../envs/peaks.yaml"
    threads: config["resources"]["idr"]["cpu"]
    resources:
        runtime=config["resources"]["peaks"]["time"],
        mem_mb=16000,
    shell:
        "perl {input.script} "
        "{input.ip_bam} "
        "{input.input_bam} "
        "{input.peaks} "
        "{input.ip_num} "
        "{input.input_num} "
        "{output.bed} "
        "> {log} 2>&1"


rule reproducible_peaks:
    input:
        rep1="results/idr/{s1}_vs_{s2}/{s1}.idr_peaks.normed.bed.full",
        rep2="results/idr/{s1}_vs_{s2}/{s2}.idr_peaks.normed.bed.full",
        entropy1="results/idr/entropy/{s1}.entropy.full",
        entropy2="results/idr/entropy/{s2}.entropy.full",
        idr="results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.idr.out",
        script=yeolab_script("get_reproducing_peaks.pl"),
    output:
        bed=temp(
            "results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.reproducible_peaks.unsorted.bed"
        ),
        custombed="results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.reproducible_peaks.custombed",
        rep1_full="results/idr/{s1}_vs_{s2}/{s1}.reproducible_peaks.full",
        rep2_full="results/idr/{s1}_vs_{s2}/{s2}.reproducible_peaks.full",
    log:
        "logs/idr/reproducible_peaks/{s1}_vs_{s2}.log",
    conda:
        "../envs/peaks.yaml"
    threads: 1
    resources:
        runtime=config["resources"]["idr"]["time"],
    shell:
        "perl {input.script} "
        "{input.rep1} "
        "{input.rep2} "
        "{output.rep1_full} "
        "{output.rep2_full} "
        "{output.bed} "
        "{output.custombed} "
        "{input.entropy1} "
        "{input.entropy2} "
        "{input.idr} "
        "> {log} 2>&1"


rule sort_reproducible_peaks:
    input:
        "results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.reproducible_peaks.unsorted.bed",
    output:
        "results/idr/{s1}_vs_{s2}/{s1}_vs_{s2}.reproducible_peaks.bed",
    log:
        "logs/idr/sort_reproducible_peaks/{s1}_vs_{s2}.log",
    conda:
        "../envs/peaks.yaml"
    threads: 1
    resources:
        runtime=15,
    shell:
        "sort -k1,1 -k2,2n {input} > {output} 2> {log}"
