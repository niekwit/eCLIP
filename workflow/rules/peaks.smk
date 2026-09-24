# ENCODE: call peak clusters with CLIPper, normalise reads in peaks over the size-matched input (Fisher
# exact/chi-square test and log2 fold change), compress overlapping peaks and remove blacklisted regions
rule clipper:
    input:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
    output:
        "results/clipper/{sample}.peakClusters.bed",
    params:
        species=resources.clipper_species,
        extra=config["clipper"]["extra"],
    threads: config["resources"]["clipper"]["cpu"]
    resources:
        runtime=config["resources"]["clipper"]["time"],
    log:
        "logs/clipper/{sample}.log",
    conda:
        "../envs/clipper.yaml"
    shell:
        "clipper "
        "--species {params.species} "
        "--bam {input.bam} "
        "--outfile {output} "
        "--processors {threads} "
        "{params.extra} "
        "> {log} 2>&1"


rule input_normalisation:
    input:
        unpack(peaks_input),
        peaks="results/clipper/{sample}.peakClusters.bed",
        script=yeolab_script("overlap_peakfi_with_bam.pl"),
    output:
        bed="results/peaks/{sample}/{sample}.normed.bed",
        full="results/peaks/{sample}/{sample}.normed.bed.full",
    threads: config["resources"]["peaks"]["cpu"]
    resources:
        runtime=config["resources"]["peaks"]["time"],
        mem_mb=16000,
    log:
        "logs/input_normalisation/{sample}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        "perl {input.script} "
        "{input.ip_bam} "
        "{input.input_bam} "
        "{input.peaks} "
        "{input.ip_num} "
        "{input.input_num} "
        "{output.bed} "
        "> {log} 2>&1"


rule compress_peaks:
    input:
        bed="results/peaks/{sample}/{sample}.normed.bed",
        script=yeolab_script(
            "compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl"
        ),
    output:
        "results/peaks/{sample}/{sample}.normed.compressed.bed",
    threads: 1
    resources:
        runtime=config["resources"]["peaks"]["time"],
    log:
        "logs/compress_peaks/{sample}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        "perl {input.script} {input.bed} {output} > {log} 2>&1"


rule sort_peaks:
    input:
        "results/{prefix}.bed",
    output:
        temp("results/{prefix}.sorted.bed"),
    wildcard_constraints:
        prefix=".+compressed",
    threads: 1
    resources:
        runtime=15,
    log:
        "logs/sort_peaks/{prefix}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        "sort -k1,1 -k2,2n {input} > {output} 2> {log}"


rule remove_blacklisted_regions:
    input:
        bed="results/peaks/{sample}/{sample}.normed.compressed.sorted.bed",
        blacklist=resources.blacklist,
    output:
        "results/peaks/{sample}/{sample}.peaks.bed",
    params:
        # ENCODE removes peaks on the same strand as the (BED6) eCLIP blacklist regions
        strand="-s" if resources.blacklist_stranded else "",
    threads: 1
    resources:
        runtime=15,
    log:
        "logs/blacklist/{sample}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        "bedtools intersect -v {params.strand} -a {input.bed} -b {input.blacklist} > {output} 2> {log}"


rule narrowpeak:
    # narrowPeak file and BED file that can be converted to bigBed
    input:
        "results/{prefix}.bed",
    output:
        narrowpeak="results/{prefix}.narrowPeak",
        fixed=temp("results/{prefix}.fixed.bed"),
    params:
        db=genome,
        l10p=config["peaks"]["l10p"],
        l2fc=config["peaks"]["l2fc"],
    wildcard_constraints:
        prefix=".+peaks",
    threads: 1
    resources:
        runtime=15,
    log:
        "logs/narrowpeak/{prefix}.log",
    conda:
        "../envs/peaks.yaml"
    script:
        "../scripts/peak_files.py"


rule bigbed:
    input:
        bed="results/{prefix}.fixed.bed",
        cs=resources.chrom_sizes,
    output:
        "results/{prefix}.bb",
    wildcard_constraints:
        prefix=".+peaks",
    threads: 1
    resources:
        runtime=15,
    log:
        "logs/bigbed/{prefix}.log",
    conda:
        "../envs/peaks.yaml"
    shell:
        "bedToBigBed {input.bed} {input.cs} {output} > {log} 2>&1"


rule total_entropy:
    input:
        unpack(peaks_input),
        full="results/peaks/{sample}/{sample}.normed.bed.full",
    output:
        "results/peaks/{sample}/{sample}.total_entropy.txt",
    params:
        l10p=config["peaks"]["l10p"],
        l2fc=config["peaks"]["l2fc"],
    threads: 1
    resources:
        runtime=15,
    log:
        "logs/total_entropy/{sample}.log",
    conda:
        "../envs/peaks.yaml"
    script:
        "../scripts/total_entropy.py"
