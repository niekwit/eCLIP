rule calculate_reads_distribution:
    input:
        r1=expand("reads/{sample}_R1_001.fastq.gz", sample=SAMPLES), # raw data
        log=expand("results/mapped/{sample}/{sample}_Log.final.out", sample=SAMPLES), # STAR log
        ddup=expand("results/mapped/{sample}/{sample}_sorted.dedup.bam", sample=SAMPLES), # reads after barcode collapse = usable reads
    output:
        "results/samtools/usable_reads.csv"
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/usable_reads.yaml"
    log:
        "logs/qc/reads_distribution.log"
    script:
        "../scripts/calculate_reads_distribution.py"


rule plot_reads_distribution:
    input:
        "results/samtools/usable_reads.csv"
    output:
        "results/samtools/usable_reads.pdf"
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/usable_reads.yaml"
    log:
        "logs/qc/plot_reads_distribution.log"
    script:
        "../scripts/plot_reads_distribution.R"
    