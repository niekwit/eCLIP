rule bigwig:
    input:
        bam="results/mapped/{sample}/{sample}_R2.bam",
        bai="results/mapped/{sample}/{sample}_R2.bam.bai",
    output:
        "results/bigwig/{sample}.bw",
    params:
        binSize=config["bigwig"]["binSize"],
        normalise=config["bigwig"]["normalise"],
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/bigwig/{sample}.log",
    conda:
        "../envs/deeptools.yaml",
    shell:
        "bamCoverage " 
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.binSize} "
        "--normalizeUsing {params.normalise} "
        "--numberOfProcessors {threads} "
        "> {log} 2>&1"
    