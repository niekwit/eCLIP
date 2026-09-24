# ENCODE: sort reads (to reduce randomness), map to repeat elements with STAR and keep the reads that do not map
rule sort_trimmed_fastq:
    input:
        "results/trimmed/round2/{unit}.r{end}.fq.gz",
    output:
        temp("results/trimmed/sorted/{unit}.r{end}.fq.gz"),
    log:
        "logs/fastq_sort/trimmed/{unit}.r{end}.log",
    conda:
        "../envs/mapping.yaml"
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"],
    shell:
        "pigz -dc {input} | "
        "fastq-sort --id 2> {log} | "
        "pigz -p {threads} > {output}"


rule star_repeats:
    input:
        idx=resources.star_repeat_index,
        reads=lambda wildcards: expand(
            "results/trimmed/sorted/{unit}.r{end}.fq.gz",
            unit=wildcards.unit,
            end=ENDS,
        ),
    output:
        bam=temp("results/star/repeats/{unit}.Aligned.out.bam"),
        log_final="results/star/repeats/{unit}.Log.final.out",
        unmapped=temp(
            expand("results/star/repeats/{{unit}}.Unmapped.out.mate{end}", end=ENDS)
        ),
    log:
        "logs/star/repeats/{unit}.log",
    conda:
        "../envs/mapping.yaml"
    threads: config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["mapping"]["time"],
        mem_mb=32000,
    params:
        prefix="results/star/repeats/{unit}.",
        extra=config["star"]["repeats_extra"],
    shell:
        "STAR "
        "--runMode alignReads "
        "--runThreadN {threads} "
        "--genomeDir {input.idx} "
        "--genomeLoad NoSharedMemory "
        "--readFilesCommand zcat "
        "--readFilesIn {input.reads} "
        "--alignEndsType EndToEnd "
        "--outSAMunmapped Within "
        "--outFilterMultimapNmax 30 "
        "--outFilterMultimapScoreRange 1 "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype BAM Unsorted "
        "--outFilterType BySJout "
        "--outBAMcompression 10 "
        "--outReadsUnmapped Fastx "
        "--outFilterScoreMin 10 "
        "--outSAMattrRGline ID:foo "
        "--outSAMattributes All "
        "--outSAMmode Full "
        "--outStd Log "
        "{params.extra} "
        "> {log} 2>&1"


rule sort_unmapped_fastq:
    # STAR does not always output mates in the same order
    input:
        "results/star/repeats/{unit}.Unmapped.out.mate{end}",
    output:
        temp("results/repeats/unmapped/{unit}.r{end}.fq.gz"),
    log:
        "logs/fastq_sort/unmapped/{unit}.r{end}.log",
    conda:
        "../envs/mapping.yaml"
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"],
    shell:
        "fastq-sort --id {input} 2> {log} | " "pigz -p {threads} > {output}"
