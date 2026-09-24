if PAIRED_END:

    rule fastqc_raw:
        input:
            "reads/{sample}_R{end}_001.fastq.gz",
        output:
            html="results/qc/fastqc/raw/{sample}_R{end}.html",
            zip="results/qc/fastqc/raw/{sample}_R{end}_fastqc.zip",
        log:
            "logs/fastqc/raw/{sample}_R{end}.log",
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb=2048,
        params:
            extra="--quiet",
        wrapper:
            f"{wrapper_version}/bio/fastqc"

else:

    rule fastqc_raw:
        input:
            "reads/{sample}.fastq.gz",
        output:
            html="results/qc/fastqc/raw/{sample}.html",
            zip="results/qc/fastqc/raw/{sample}_fastqc.zip",
        log:
            "logs/fastqc/raw/{sample}.log",
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb=2048,
        params:
            extra="--quiet",
        wrapper:
            f"{wrapper_version}/bio/fastqc"


rule fastqc_trimmed:
    # ENCODE runs FastQC after both rounds of cutadapt
    input:
        "results/trimmed/round{round}/{unit}.r{end}.fq.gz",
    output:
        html="results/qc/fastqc/trim{round}/{unit}.r{end}.html",
        zip="results/qc/fastqc/trim{round}/{unit}.r{end}_fastqc.zip",
    log:
        "logs/fastqc/trim{round}/{unit}.r{end}.log",
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        runtime=config["resources"]["fastqc"]["time"],
        mem_mb=2048,
    params:
        extra="--quiet",
    wrapper:
        f"{wrapper_version}/bio/fastqc"


rule multiqc:
    input:
        multiqc_input(),
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="MultiQC",
        ),
    log:
        "logs/multiqc/multiqc.log",
    conda:
        "../envs/mapping.yaml"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        runtime=config["resources"]["fastqc"]["time"],
        mem_mb=2048,
    params:
        dir=lambda wildcards, output: os.path.dirname(output[0]),
        # Sample names are prefixed with their directory (raw/trim1/trim2, round1/round2, repeats/genome)
        extra="--dirs --dirs-depth 1",
    shell:
        "multiqc "
        "--force "
        "--outdir {params.dir} "
        "-n multiqc.html "
        "{params.extra} "
        "{input} "
        "> {log} 2>&1"
