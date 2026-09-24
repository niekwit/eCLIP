# ENCODE: map reads that do not map to repeat elements to the genome with STAR (unique mappers only),
# remove PCR duplicates using the UMIs and (paired-end) merge the two inline barcodes and keep read 2
rule star_genome:
    input:
        unpack(star_input),
        idx=resources.star_index,
    output:
        bam=temp("results/star/genome/{unit}.Aligned.out.bam"),
        log_final="results/star/genome/{unit}.Log.final.out",
        unmapped=temp(
            expand("results/star/genome/{{unit}}.Unmapped.out.mate{end}", end=ENDS)
        ),
    log:
        "logs/star/genome/{unit}.log",
    conda:
        "../envs/mapping.yaml"
    threads: config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["mapping"]["time"],
        mem_mb=40000,
    params:
        prefix="results/star/genome/{unit}.",
        extra=config["star"]["genome_extra"],
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
        "--outFilterMultimapNmax 1 "
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


if PAIRED_END:

    rule namesort_pe:
        # Read pairs must be adjacent for barcodecollapsepe.py
        input:
            "results/star/genome/{unit}.Aligned.out.bam",
        output:
            temp("results/mapped/genome/{unit}.namesorted.bam"),
        log:
            "logs/samtools_sort/namesort/{unit}.log",
        conda:
            "../envs/mapping.yaml"
        threads: config["resources"]["samtools"]["cpu"]
        resources:
            runtime=config["resources"]["samtools"]["time"],
        shell:
            "samtools sort -n -@ {threads} -o {output} {input} 2> {log}"

    rule barcode_collapse_pe:
        # Randomer aware PCR duplicate removal (Yeo lab script)
        input:
            bam="results/mapped/genome/{unit}.namesorted.bam",
            script=yeolab_script("barcodecollapsepe.py"),
        output:
            bam=temp("results/mapped/genome/{unit}.rmdup.bam"),
            metrics="results/qc/barcode_collapse/{unit}.metrics",
        log:
            "logs/barcode_collapse/{unit}.log",
        conda:
            "../envs/yeolab_py2.yaml"
        threads: 1
        resources:
            runtime=config["resources"]["umi_tools"]["time"],
            mem_mb=32000,
        shell:
            "python {input.script} "
            "-b {input.bam} "
            "-o {output.bam} "
            "-m {output.metrics} "
            "> {log} 2>&1"

    rule sort_rmdup_pe:
        input:
            "results/mapped/genome/{unit}.rmdup.bam",
        output:
            temp("results/mapped/genome/{unit}.rmdup.sorted.bam"),
        log:
            "logs/samtools_sort/rmdup/{unit}.log",
        conda:
            "../envs/mapping.yaml"
        threads: config["resources"]["samtools"]["cpu"]
        resources:
            runtime=config["resources"]["samtools"]["time"],
        shell:
            "samtools sort -@ {threads} -o {output} {input} 2> {log}"

    rule merge_barcodes_pe:
        # Merges the technical replicates (the two inline barcodes)
        input:
            lambda wildcards: expand(
                "results/mapped/genome/{unit}.rmdup.sorted.bam",
                unit=[
                    f"{wildcards.sample}.{x}"
                    for x in sample_barcodes(wildcards.sample)
                ],
            ),
        output:
            temp("results/mapped/merged/{sample}.bam"),
        log:
            "logs/samtools_merge/{sample}.log",
        conda:
            "../envs/mapping.yaml"
        threads: config["resources"]["samtools"]["cpu"]
        resources:
            runtime=config["resources"]["samtools"]["time"],
        shell:
            "samtools merge -@ {threads} {output} {input} 2> {log}"

    rule select_read2_pe:
        # Only read 2 is used with the single stranded peak caller (final BAM file)
        input:
            "results/mapped/merged/{sample}.bam",
        output:
            "results/mapped/{sample}.bam",
        log:
            "logs/samtools_view/read2/{sample}.log",
        conda:
            "../envs/mapping.yaml"
        threads: config["resources"]["samtools"]["cpu"]
        resources:
            runtime=config["resources"]["samtools"]["time"],
        shell:
            "samtools view -f 128 -b -@ {threads} -o {output} {input} 2> {log}"

else:

    rule sort_genome_bam_se:
        # ENCODE sorts by name and then by position to make the order of reads deterministic
        input:
            "results/star/genome/{unit}.Aligned.out.bam",
        output:
            temp("results/mapped/genome/{unit}.sorted.bam"),
        log:
            "logs/samtools_sort/{unit}.log",
        conda:
            "../envs/mapping.yaml"
        threads: config["resources"]["samtools"]["cpu"]
        resources:
            runtime=config["resources"]["samtools"]["time"],
        shell:
            "samtools sort -n -u -@ {threads} {input} 2> {log} | "
            "samtools sort -@ {threads} -o {output} - 2>> {log}"

    rule umi_dedup_se:
        input:
            bam="results/mapped/genome/{unit}.sorted.bam",
            bai="results/mapped/genome/{unit}.sorted.bam.bai",
        output:
            bam=temp("results/mapped/dedup/{unit}.bam"),
            stats=(
                multiext(
                    "results/qc/umi_tools/{unit}",
                    "_edit_distance.tsv",
                    "_per_umi.tsv",
                    "_per_umi_per_position.tsv",
                )
                if config["umi_tools"]["dedup_stats"]
                else []
            ),
        log:
            "logs/umi_tools/dedup/{unit}.log",
        conda:
            "../envs/umi_tools.yaml"
        threads: config["resources"]["umi_tools"]["cpu"]
        resources:
            runtime=config["resources"]["umi_tools"]["time"],
            mem_mb=32000,
        params:
            stats=lambda wildcards: (
                f"--output-stats results/qc/umi_tools/{wildcards.unit}"
                if config["umi_tools"]["dedup_stats"]
                else ""
            ),
        shell:
            "umi_tools dedup "
            "--random-seed 1 "
            "-I {input.bam} "
            "--method unique "
            "{params.stats} "
            "-S {output.bam} "
            "> {log} 2>&1"

    rule sort_dedup_se:
        input:
            "results/mapped/dedup/{sample}.bam",
        output:
            "results/mapped/{sample}.bam",
        log:
            "logs/samtools_sort/dedup/{sample}.log",
        conda:
            "../envs/mapping.yaml"
        threads: config["resources"]["samtools"]["cpu"]
        resources:
            runtime=config["resources"]["samtools"]["time"],
        shell:
            "samtools sort -@ {threads} -o {output} {input} 2> {log}"


rule bam_index:
    input:
        "results/{path}.bam",
    output:
        "results/{path}.bam.bai",
    log:
        "logs/samtools_index/{path}.log",
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"],
    params:
        extra="",  # optional params string
    wrapper:
        f"{wrapper_version}/bio/samtools/index"


rule mapped_readnum:
    # Number of mapped reads of the final BAM file (used for input normalisation)
    input:
        "results/mapped/{sample}.bam",
    output:
        "results/mapped/{sample}.readnum.txt",
    log:
        "logs/samtools_view/readnum/{sample}.log",
    conda:
        "../envs/mapping.yaml"
    threads: 1
    resources:
        runtime=15,
    shell:
        "samtools view -c -F 4 {input} > {output} 2> {log}"
