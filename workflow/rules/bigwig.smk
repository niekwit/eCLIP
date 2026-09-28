# ENCODE: RPM normalised, strand specific read density files (Yeo lab makebigwigfiles.py,
# reimplemented with bedtools genomecov and UCSC tools)
rule bigwig_pos:
    input:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
        cs=resources.chrom_sizes,
    output:
        "results/bigwig/{sample}.norm.pos.bw",
    log:
        "logs/bigwig/{sample}.pos.log",
    conda:
        "../envs/bigwig.yaml"
    threads: config["resources"]["bigwig"]["cpu"]
    resources:
        runtime=config["resources"]["bigwig"]["time"],
    params:
        # Single-end reads are not reversed (direction f), paired-end read 2 is (direction r)
        strand="-" if PAIRED_END else "+",
        sign=1,
    script:
        "../scripts/bigwig.py"


use rule bigwig_pos as bigwig_neg with:
    output:
        "results/bigwig/{sample}.norm.neg.bw",
    log:
        "logs/bigwig/{sample}.neg.log",
    params:
        strand="+" if PAIRED_END else "-",
        sign=-1,
