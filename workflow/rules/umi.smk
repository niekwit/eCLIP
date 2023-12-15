rule umi_tools: # adds barcode sequence to read name
    input:
        r1="reads/{sample}_R1_001.fastq.gz",
        r2="reads/{sample}_R2_001.fastq.gz",
    output:
        r1="results/umi_tools/{sample}_R1_001.fastq.gz",
        r2="results/umi_tools/{sample}_R2_001.fastq.gz",
    params:
        bc_pattern=bc_pattern,
    conda:
        "../envs/trim.yaml",
    threads: config["resources"]["trim"]["cpu"]
    log: 
        "logs/umi_tools/{sample}.log",
    shell:
        "umi_tools extract "
        "--bc-pattern={params.bc_pattern} "
        "-I {input.r2} " # barcode is 3' end of read2, so read2 and read1 are swapped
        "--stdout {output.r2} "
        "--read2-in={input.r1} "
        "--read2-out={output.r1} "
        "--log2stderr True "
        "2> {log}"




