rule get_usable_reads:
    input:
        "{sample}.bam" # PLACEHOLDER, SHOULD BE MANY BAM FILES
    output:
        "results/samtools/usable_reads.csv"
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/usable_reads.yaml"
    script:
        "../scripts/usable_reads.py"


rule plot_usable_reads:
    input:
        csv = "results/samtools/usable_reads.csv"
    output:
        "results/samtools/usable_reads.pdf"
    conda:
        "../envs/usable_reads.yaml"
    script:
        "../scripts/plot_usable_reads.R"
    