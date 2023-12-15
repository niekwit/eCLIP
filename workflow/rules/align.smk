'''
rule hisat2_index:
    input:
        fasta=resources.fasta
    output:
        directory("resources/index"),
    params:
        prefix=f"resources/index/index"
    log:
        "logs/hisat2_index/index.log"
    threads: config["resources"]["index"]["cpu"]
    resources: 
        runtime=config["resources"]["index"]["time"]
    wrapper:
        "v3.1.0/bio/hisat2/index"


rule mapping:
    input:
        reads=["results/cutadapt2/{sample}_R1.fastq.gz", "results/cutadapt2/{sample}_R2.fastq.gz"],
        idx="resources/index",
    output:
        "results/mapped/{sample}.bam",
    log:
        "logs/hisat2_align/{sample}.log",
    params:
        extra="",
    threads: config["resources"]["mapping"]["cpu"]
    resources: 
        runtime=config["resources"]["mapping"]["time"]
    wrapper:
        "v3.1.0/bio/hisat2/align"


rule remove_multi_mapping_reads:
    input:
        "results/mapped/{sample}.bam",
    output:
        "results/mapped/{sample}_unique.bam",
    log:
        "logs/samtools/view_{sample}.log",
    params:
        extra="",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view {params.extra} -h {input} | grep -P '(NH:i:1|^@)' | samtools view -Sb - > {output}"


rule sort:
    input:
        "results/mapped/{sample}_unique.bam",
    output:
        "results/mapped/{sample}_sorted.bam",
    log:
        "logs/samtools/sort_{sample}.log",
    params:
        extra="",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    wrapper:
        "v3.1.0/bio/samtools/sort"
'''
rule get_readlength:  
    input:
        "results/qc/multiqc_raw_data/multiqc_data/multiqc_general_stats.txt",
    output:
        "results/qc/readlength.txt",
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/mapping/readlength.log"
    script:
        "../scripts/get_readlength.sh"


rule genome_index:
    input:
        fa=resources.fasta,
        gtf=resources.gtf,
        rl="results/qc/readlength.txt",
    output:
        directory("resources/genome_index/"),
    threads: config["resources"]["star_index"]["cpu"],
    resources:
        runtime=config["resources"]["star_index"]["time"],
    conda:
        "../envs/mapping.yaml",
    log:
        "logs/genome_index/genome_index.log"
    shell:
        "mkdir -p {output} ; "
        "STAR "
        "--runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fa} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang $(cat {input.rl}) "
        "> {log} 2>&1"


rule mapping:
    input:
        r1="results/cutadapt2/{sample}_R1_001.fastq.gz", 
        r2="results/cutadapt2/{sample}_R2_001.fastq.gz",
        idx="resources/genome_index/",
    output:
        "results/mapped/{sample}/{sample}_Aligned.out.bam",
    params:
        extra=config["STAR"]["extra"]
    threads: config["resources"]["mapping"]["cpu"],
    resources:
        runtime=config["resources"]["mapping"]["time"],
    conda:
        "../envs/mapping.yaml",
    log:
        "logs/mapping/{sample}.log"
    shell:
        "rm -rf temp_{wildcards.sample}/ ;"
        "STAR "
        "--runThreadN {threads} "
        "--runMode alignReads "
        "--genomeDir {input.idx} "
        "--readFilesCommand zcat "
        "--readFilesIn {input.r1} {input.r2} "
        "--outSAMunmapped Within "
        "--outFilterMultimapNmax 1 " # removes multimapping reads when set to 1
        "--outFilterMultimapScoreRange 1 " # the score range below the maximum score for multimapping alignments
        "--outFileNamePrefix results/mapped/{wildcards.sample}/{wildcards.sample}_ "
        "--outSAMattributes All "
        "--outSAMtype BAM Unsorted "
        "--outFilterType BySJout " # reduces the number of ”spurious” junctions
        "--outReadsUnmapped Fastx "
        "--outFilterScoreMin 10 "
        "--outSAMattrRGline ID:{wildcards.sample} "
        "--outStd Log "
        "--alignEndsType EndToEnd "
        "--outBAMcompression 10 "
        "--outSAMmode Full "
        "--outTmpDir temp_{wildcards.sample}/ "
        "> {log} 2>&1"


rule sort:
    input:
        "results/mapped/{sample}/{sample}_Aligned.out.bam",
    output:
        "results/mapped/{sample}/{sample}_sorted.bam",
    log:
        "logs/samtools/sort_{sample}.log",
    params:
        extra="-m 4G", # sort by name to ensure read pairs are adjacent
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    wrapper:
        "v3.1.0/bio/samtools/sort"


