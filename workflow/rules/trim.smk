rule trim_adapters: # not the same r1/r2 as in umi_tools
    input:
        r1="results/umi_tools/{sample}_R1_001.fastq.gz",
        r2="results/umi_tools/{sample}_R2_001.fastq.gz",
    output:
        r1="results/cutadapt1/{sample}_R1_001.fastq.gz",
        r2="results/cutadapt1/{sample}_R2_001.fastq.gz",
    conda:
        "../envs/trim.yaml",
    threads: config["resources"]["trim"]["cpu"]
    log: 
        "logs/cutadapt/{sample}_adapter_trim.log",
    shell:
        "cutadapt "
        "--cores {threads} "
        "--match-read-wildcards "
        "--times 1 " # default value, remove?
        "-e 0.1 " # default value, remove?
        "--quality-cutoff 6,6 "
        "-U 3 " # removes n bases from r2
        "-m 18 " # remove reads shorter than n bp from r1
        "-g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT " # 5' adapter sequencess r1
        "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC " # 3' adapter sequence r1
        "-A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT " # 3' adapter sequences r2
        "-o {output.r1} "
        "-p {output.r2} "
        "{input.r1} {input.r2} "
        "> {log} 2>&1"


rule trim_double_ligation_events:
    input:
        r1="results/cutadapt1/{sample}_R1_001.fastq.gz",
        r2="results/cutadapt1/{sample}_R2_001.fastq.gz",
    output:
        r1="results/cutadapt2/{sample}_R1_001.fastq.gz",
        r2="results/cutadapt2/{sample}_R2_001.fastq.gz",
    conda:
        "../envs/trim.yaml",
    threads: config["resources"]["trim"]["cpu"]
    log: 
        "logs/cutadapt/{sample}_double_ligation_trim.log",
    shell:
        "cutadapt "
        "--cores {threads} "
        "--match-read-wildcards "
        "--times 1 " # default value, remove?
        "-e 0.1 " # default value, remove?
        "-O 5 " # minimum overlap length
        "--quality-cutoff 6 "
        "-m 18 " # remove reads shorter than n bp from r1
        "-A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT " # 3' adapter sequences r2
        "-o {output.r1} "
        "-p {output.r2} "
        "{input.r1} {input.r2} "
        "> {log} 2>&1"


