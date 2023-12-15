rule fastqc_raw_data:
    input:
        "reads/{sample}_{end}.fastq.gz"
    output:
        html="results/qc/fastqc_raw/{sample}_{end}.html",
        zip="results/qc/fastqc_raw/{sample}_{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/qc/fastqc_raw_data/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        mem_mb = 1024
    wrapper:
        "v3.1.0/bio/fastqc"


rule multiqc_raw_data:
    input:
        expand("results/qc/fastqc_raw/{sample}_{end}_fastqc.zip", sample=SAMPLES, end=["R1_001","R2_001"]),
    output:
        r="results/qc/multiqc_raw_data/multiqc.html",
        d=directory("results/qc/multiqc_raw_data"),
        t="results/qc/multiqc_raw_data/multiqc_data/multiqc_general_stats.txt",
    params:
        extra=""
    log:
        "logs/qc/multiqc/multiqc_raw_data.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "multiqc " 
        "--force "
        "--outdir {output.d} "
        "-n multiqc.html "
        "{params.extra} "
        "{input} "
        "> {log} 2>&1"


rule fastqc_trim_adapters:
    input:
        "results/cutadapt1/{sample}_{end}.fastq.gz"
    output:
        html="results/qc/fastqc1/{sample}_{end}.html",
        zip="results/qc/fastqc1/{sample}_{end}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/qc/fastqc1/{sample}_{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        mem_mb=1024
    wrapper:
        "v3.1.0/bio/fastqc"


rule fastqc_trim_double_ligation_events:
    input:
        "results/cutadapt2/{sample}_{end}.fastq.gz"
    output:
        html="results/qc/fastqc2/{sample}_{end}.html",
        zip="results/qc/fastqc2/{sample}_{end}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
    params:
        extra="--quiet"
    log:
        "logs/qc/fastqc2/{sample}_{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        mem_mb=1024
    wrapper:
        "v3.1.0/bio/fastqc"


rule multiqc_trimmed_reads:
    input:
        expand("results/qc/fastqc{exp}/{sample}_{end}_fastqc.zip", sample=SAMPLES, end=["R1_001","R2_001"],exp=[1,2]),
    output:
        r="results/qc/multiqc_trimmed_reads/multiqc.html",
        d=directory("results/qc/multiqc_trimmed_reads/"),
        t="results/qc/multiqc_trimmed_reads/multiqc_data/multiqc_general_stats.txt",
    params:
        extra=""
    log:
        "logs/qc/multiqc/multiqc_trimmed_reads.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "multiqc " 
        "--force "
        "--outdir {output.d} "
        "-n multiqc.html "
        "{params.extra} "
        "{input} "
        "> {log} 2>&1"


rule samtools_flagstat:
    input:
        "results/mapped/{sample}/{sample}_sorted.dedup.bam",
    output:
        "results/mapped/{sample}/{sample}_sorted.dedup.bam.flagstat",
    log:
        "logs/samtools/flagstat_{sample}.log",
    threads: config["resources"]["samtools"]["cpu"],
    resources:
        runtime=config["resources"]["samtools"]["time"],
    wrapper:
        "v3.1.0/bio/samtools/flagstat"


rule multiqc_dedup_flagstat:
    input:
        expand("results/mapped/{sample}/{sample}_sorted.dedup.bam.flagstat", sample=SAMPLES),
    output:
        r="results/qc/multiqc_dedup_flagstat/multiqc.html",
        d=directory("results/qc/multiqc_dedup_flagstat/"),
        t="results/qc/multiqc_dedup_flagstat/multiqc_data/multiqc_general_stats.txt",
    params:
        extra=""
    log:
        "logs/qc/multiqc/multiqc_dedup_flagstat.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "multiqc " 
        "--force "
        "--outdir {output.d} "
        "-n multiqc.html "
        "{params.extra} "
        "{input} "
        "> {log} 2>&1"