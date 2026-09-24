# ENCODE: RPM normalised, strand specific read density files (Yeo lab makebigwigfiles.py,
# reimplemented with bedtools genomecov and UCSC tools)
rule bigwig:
    input:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
        cs=resources.chrom_sizes,
    output:
        pos="results/bigwig/{sample}.norm.pos.bw",
        neg="results/bigwig/{sample}.norm.neg.bw",
    log:
        "logs/bigwig/{sample}.log",
    conda:
        "../envs/bigwig.yaml"
    threads: config["resources"]["bigwig"]["cpu"]
    resources:
        runtime=config["resources"]["bigwig"]["time"],
    params:
        # Single-end reads are not reversed (direction f), paired-end read 2 is (direction r)
        pos_strand="-" if PAIRED_END else "+",
        neg_strand="+" if PAIRED_END else "-",
    script:
        "../scripts/bigwig.sh"
