rule index_sorted_bam:
    input:
        "results/mapped/{sample}/{sample}_sorted.bam",
    output:
        "results/mapped/{sample}/{sample}_sorted.bam.bai",
    log:
        "logs/samtools/index_sorted_{sample}.log",
    params:
        extra="",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    wrapper:
        "v3.3.3/bio/samtools/index"


rule deduplication:
    input:
        bam="results/mapped/{sample}/{sample}_sorted.bam",
        bai="results/mapped/{sample}/{sample}_sorted.bam.bai",
    output:
        "results/mapped/{sample}/{sample}_sorted.dedup.bam",
    params:
        stats_prefix="results/mapped/{sample}/{sample}_sorted.dedup.stats"
    log:
        "logs/dedup/{sample}.log",
    threads: config["resources"]["trim"]["cpu"],
    resources:
        runtime=config["resources"]["trim"]["time"],
    conda:
        "../envs/trim.yaml",
    shell:
        "umi_tools dedup "
        "--output-stats={params.stats_prefix} "
        "--paired "
        "-I {input.bam} "
        "-S {output} "
        "> {log} 2>&1"

