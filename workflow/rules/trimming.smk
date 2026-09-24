# ENCODE: identify UMIs (single-end: umi_tools; paired-end: eclipdemux, which also demultiplexes the
# inline barcodes) and trim adapters twice with cutadapt (to remove double ligation events)
if PAIRED_END:

    rule demux_pe:
        input:
            r1="reads/{sample}_R1_001.fastq.gz",
            r2="reads/{sample}_R2_001.fastq.gz",
            barcodes="resources/yeolab/barcodes.fasta",
        output:
            dir=temp(directory("results/demux/{sample}")),
            metrics="results/qc/demux/{sample}.metrics",
        log:
            "logs/demux/{sample}.log",
        conda:
            "../envs/yeolab_py2.yaml"
        threads: 1
        resources:
            runtime=config["resources"]["trim"]["time"],
        params:
            barcode_a=lambda wildcards: csv.loc[
                csv["sample"] == wildcards.sample, "barcode_a"
            ].iloc[0],
            barcode_b=lambda wildcards: csv.loc[
                csv["sample"] == wildcards.sample, "barcode_b"
            ].iloc[0],
            length=config["umi"]["pe_length"],
        shell:
            # eclipdemux writes its output files to the working directory
            "R1=$(realpath {input.r1}); "
            "R2=$(realpath {input.r2}); "
            "BARCODES=$(realpath {input.barcodes}); "
            "METRICS=$(realpath -m {output.metrics}); "
            "LOG=$(realpath -m {log}); "
            "mkdir -p {output.dir} && cd {output.dir} && "
            "demux "
            "--metrics $METRICS "
            "--expectedbarcodeida {params.barcode_a} "
            "--expectedbarcodeidb {params.barcode_b} "
            "--fastq_1 $R1 "
            "--fastq_2 $R2 "
            "--newname {wildcards.sample} "
            "--dataset eclip "
            "--barcodesfile $BARCODES "
            "--length {params.length} "
            "> $LOG 2>&1"

    rule cutadapt_pe_round1:
        input:
            "results/demux/{sample}",
        output:
            r1=temp("results/trimmed/round1/{sample}.{barcode}.r1.fq.gz"),
            r2=temp("results/trimmed/round1/{sample}.{barcode}.r2.fq.gz"),
            qc="results/qc/cutadapt/round1/{sample}.{barcode}.txt",
        log:
            "logs/cutadapt/{sample}.{barcode}.round1.log",
        conda:
            "../envs/trimming.yaml"
        threads: config["resources"]["trim"]["cpu"]
        resources:
            runtime=config["resources"]["trim"]["time"],
        params:
            adapters=lambda wildcards: cutadapt_args(
                f"{wildcards.sample}.{wildcards.barcode}", 1
            ),
            error_rate=config["cutadapt"]["error_rate"],
            quality_cutoff=config["cutadapt"]["quality_cutoff"],
            min_length=config["cutadapt"]["min_length"],
        shell:
            "cutadapt "
            "--match-read-wildcards "
            "--times 1 "
            "-e {params.error_rate} "
            "-O 1 "
            "--quality-cutoff {params.quality_cutoff} "
            "-m {params.min_length} "
            "-j {threads} "
            "{params.adapters} "
            "-o {output.r1} "
            "-p {output.r2} "
            "{input}/eclip.{wildcards.sample}.{wildcards.barcode}.r1.fq.gz "
            "{input}/eclip.{wildcards.sample}.{wildcards.barcode}.r2.fq.gz "
            "> {output.qc} 2> {log}"

    rule cutadapt_pe_round2:
        input:
            r1="results/trimmed/round1/{sample}.{barcode}.r1.fq.gz",
            r2="results/trimmed/round1/{sample}.{barcode}.r2.fq.gz",
        output:
            r1=temp("results/trimmed/round2/{sample}.{barcode}.r1.fq.gz"),
            r2=temp("results/trimmed/round2/{sample}.{barcode}.r2.fq.gz"),
            qc="results/qc/cutadapt/round2/{sample}.{barcode}.txt",
        log:
            "logs/cutadapt/{sample}.{barcode}.round2.log",
        conda:
            "../envs/trimming.yaml"
        threads: config["resources"]["trim"]["cpu"]
        resources:
            runtime=config["resources"]["trim"]["time"],
        params:
            adapters=lambda wildcards: cutadapt_args(
                f"{wildcards.sample}.{wildcards.barcode}", 2
            ),
            error_rate=config["cutadapt"]["error_rate"],
            quality_cutoff=config["cutadapt"]["quality_cutoff"],
            min_length=config["cutadapt"]["min_length"],
        shell:
            "cutadapt "
            "--match-read-wildcards "
            "--times 1 "
            "-e {params.error_rate} "
            "-O 5 "
            "--quality-cutoff {params.quality_cutoff} "
            "-m {params.min_length} "
            "-j {threads} "
            "{params.adapters} "
            "-o {output.r1} "
            "-p {output.r2} "
            "{input.r1} "
            "{input.r2} "
            "> {output.qc} 2> {log}"

else:

    rule umi_extract_se:
        input:
            "reads/{unit}.fastq.gz",
        output:
            temp("results/umi/{unit}.r1.fq.gz"),
        log:
            "logs/umi_tools/extract/{unit}.log",
        conda:
            "../envs/umi_tools.yaml"
        threads: 1
        resources:
            runtime=config["resources"]["umi_tools"]["time"],
        params:
            pattern=lambda wildcards: "N" * config["umi"]["se_length"],
        shell:
            "umi_tools extract "
            "--random-seed 1 "
            "--bc-pattern {params.pattern} "
            "--log {log} "
            "--stdin {input} "
            "--stdout {output}"

    rule cutadapt_se_round1:
        input:
            "results/umi/{unit}.r1.fq.gz",
        output:
            fastq=temp("results/trimmed/round1/{unit}.r1.fq.gz"),
            qc="results/qc/cutadapt/round1/{unit}.txt",
        log:
            "logs/cutadapt/{unit}.round1.log",
        conda:
            "../envs/trimming.yaml"
        threads: config["resources"]["trim"]["cpu"]
        resources:
            runtime=config["resources"]["trim"]["time"],
        params:
            adapters=lambda wildcards: cutadapt_args(wildcards.unit, 1),
            error_rate=config["cutadapt"]["error_rate"],
            quality_cutoff=config["cutadapt"]["quality_cutoff"],
            min_length=config["cutadapt"]["min_length"],
        shell:
            "cutadapt "
            "--match-read-wildcards "
            "--times 1 "
            "-e {params.error_rate} "
            "-O 1 "
            "--quality-cutoff {params.quality_cutoff} "
            "-m {params.min_length} "
            "-j {threads} "
            "{params.adapters} "
            "-o {output.fastq} "
            "{input} "
            "> {output.qc} 2> {log}"

    rule cutadapt_se_round2:
        input:
            "results/trimmed/round1/{unit}.r1.fq.gz",
        output:
            fastq=temp("results/trimmed/round2/{unit}.r1.fq.gz"),
            qc="results/qc/cutadapt/round2/{unit}.txt",
        log:
            "logs/cutadapt/{unit}.round2.log",
        conda:
            "../envs/trimming.yaml"
        threads: config["resources"]["trim"]["cpu"]
        resources:
            runtime=config["resources"]["trim"]["time"],
        params:
            adapters=lambda wildcards: cutadapt_args(wildcards.unit, 2),
            error_rate=config["cutadapt"]["error_rate"],
            quality_cutoff=config["cutadapt"]["quality_cutoff"],
            min_length=config["cutadapt"]["min_length"],
        shell:
            "cutadapt "
            "--match-read-wildcards "
            "--times 1 "
            "-e {params.error_rate} "
            "-O 5 "
            "--quality-cutoff {params.quality_cutoff} "
            "-m {params.min_length} "
            "-j {threads} "
            "{params.adapters} "
            "-o {output.fastq} "
            "{input} "
            "> {output.qc} 2> {log}"
