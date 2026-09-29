# Locus-resolved TE binding analysis (optional, off by default -- see README.md).
#
# Reads that map to the repeat element consensus reference in the pre-filter step
# (results/star/repeats/{unit}.Aligned.out.bam, normally discarded) are the candidate pool: they
# are re-aligned to the full genome with the same uniqueness requirement as the main genome
# mapping step. A read fully internal to a repeat copy multi-maps genome-wide and is dropped;
# only reads with a unique anchor (e.g. a TE-to-flank readthrough junction) survive, which is what
# lets a specific genomic TE locus (RepeatMasker copy) be assigned with confidence.
if TE_REPEATS["enabled"]:

    if PAIRED_END:

        rule te_candidate_fastq:
            # Primary, mapped alignments only (drop unmapped/secondary/supplementary): fastq of
            # the reads that mapped to the repeat consensus reference, name-sorted so mates stay
            # paired for samtools fastq
            input:
                "results/star/repeats/{unit}.Aligned.out.bam",
            output:
                r1=temp("results/te_repeats/candidate_fastq/{unit}.r1.fq.gz"),
                r2=temp("results/te_repeats/candidate_fastq/{unit}.r2.fq.gz"),
            log:
                "logs/te_repeats/candidate_fastq/{unit}.log",
            conda:
                "../envs/mapping.yaml"
            threads: config["resources"]["samtools"]["cpu"]
            resources:
                runtime=config["resources"]["samtools"]["time"],
            shell:
                "samtools view -F 2308 -u {input} | "
                "samtools sort -n -@ {threads} -o - - 2> {log} | "
                "samtools fastq -n -1 {output.r1} -2 {output.r2} "
                "-0 /dev/null -s /dev/null - >> {log} 2>&1"

    else:

        rule te_candidate_fastq:
            input:
                "results/star/repeats/{unit}.Aligned.out.bam",
            output:
                r1=temp("results/te_repeats/candidate_fastq/{unit}.r1.fq.gz"),
            log:
                "logs/te_repeats/candidate_fastq/{unit}.log",
            conda:
                "../envs/mapping.yaml"
            threads: config["resources"]["samtools"]["cpu"]
            resources:
                runtime=config["resources"]["samtools"]["time"],
            shell:
                "samtools view -F 2308 -u {input} | "
                "samtools fastq -n -0 {output.r1} - 2> {log}"

    rule te_candidate_family:
        # Read name -> repeat consensus family (from the RNAME of the repeat-mapping BAM, e.g.
        # "DF000000001.4|MIR" -> "MIR"), for the optional family-match cross-check downstream
        input:
            "results/star/repeats/{unit}.Aligned.out.bam",
        output:
            temp("results/te_repeats/candidate_family/{unit}.tsv.gz"),
        log:
            "logs/te_repeats/candidate_family/{unit}.log",
        conda:
            "../envs/mapping.yaml"
        threads: 1
        resources:
            runtime=15,
        shell:
            "samtools view -F 2308 {input} 2> {log} | "
            """awk -F'\\t' 'BEGIN{{OFS="\\t"}}{{split($3,a,"|"); print $1,a[2]}}' | """
            "pigz > {output} 2>> {log}"

    use rule star_genome as star_te_candidates with:
        input:
            reads=expand(
                "results/te_repeats/candidate_fastq/{{unit}}.r{end}.fq.gz", end=ENDS
            ),
            idx=resources.star_index,
        output:
            bam=temp("results/star/te_candidates/{unit}.Aligned.out.bam"),
            log_final="results/star/te_candidates/{unit}.Log.final.out",
            unmapped=temp(
                expand(
                    "results/star/te_candidates/{{unit}}.Unmapped.out.mate{end}",
                    end=ENDS,
                )
            ),
        log:
            "logs/star/te_candidates/{unit}.log",
        params:
            prefix="results/star/te_candidates/{unit}.",
            extra=config["star"]["genome_extra"],

    if PAIRED_END:

        use rule namesort_pe as te_namesort_pe with:
            input:
                "results/star/te_candidates/{unit}.Aligned.out.bam",
            output:
                temp("results/te_repeats/mapped/{unit}.namesorted.bam"),
            log:
                "logs/te_repeats/samtools_sort/namesort/{unit}.log",

        use rule barcode_collapse_pe as te_barcode_collapse_pe with:
            input:
                bam="results/te_repeats/mapped/{unit}.namesorted.bam",
                script=yeolab_script("barcodecollapsepe.py"),
            output:
                bam=temp("results/te_repeats/mapped/{unit}.rmdup.bam"),
                metrics="results/qc/te_repeats/barcode_collapse/{unit}.metrics",
            log:
                "logs/te_repeats/barcode_collapse/{unit}.log",

        use rule sort_rmdup_pe as te_sort_rmdup_pe with:
            input:
                "results/te_repeats/mapped/{unit}.rmdup.bam",
            output:
                temp("results/te_repeats/mapped/{unit}.rmdup.sorted.bam"),
            log:
                "logs/te_repeats/samtools_sort/rmdup/{unit}.log",

        if not DEMULTIPLEXED:

            use rule merge_barcodes_pe as te_merge_barcodes_pe with:
                input:
                    lambda wildcards: expand(
                        "results/te_repeats/mapped/{unit}.rmdup.sorted.bam",
                        unit=[
                            f"{wildcards.sample}.{x}"
                            for x in sample_barcodes(wildcards.sample)
                        ],
                    ),
                output:
                    temp("results/te_repeats/mapped/merged/{sample}.bam"),
                log:
                    "logs/te_repeats/samtools_merge/{sample}.log",

        use rule select_read2_pe as te_select_read2_pe with:
            input:
                te_select_read2_input,
            output:
                "results/te_repeats/mapped/{sample}.bam",
            log:
                "logs/te_repeats/samtools_view/read2/{sample}.log",

    else:

        use rule sort_genome_bam_se as te_sort_genome_bam_se with:
            input:
                "results/star/te_candidates/{unit}.Aligned.out.bam",
            output:
                temp("results/te_repeats/mapped/{unit}.sorted.bam"),
            log:
                "logs/te_repeats/samtools_sort/{unit}.log",

        use rule umi_dedup_se as te_umi_dedup_se with:
            input:
                bam="results/te_repeats/mapped/{unit}.sorted.bam",
                bai="results/te_repeats/mapped/{unit}.sorted.bam.bai",
            output:
                bam=temp("results/te_repeats/mapped/dedup/{unit}.bam"),
                stats=(
                    multiext(
                        "results/qc/te_repeats/umi_tools/{unit}",
                        "_edit_distance.tsv",
                        "_per_umi.tsv",
                        "_per_umi_per_position.tsv",
                    )
                    if config["umi_tools"]["dedup_stats"]
                    else []
                ),
            log:
                "logs/te_repeats/umi_tools/dedup/{unit}.log",
            params:
                stats=lambda wildcards: (
                    f"--output-stats results/qc/te_repeats/umi_tools/{wildcards.unit}"
                    if config["umi_tools"]["dedup_stats"]
                    else ""
                ),

        use rule sort_dedup_se as te_sort_dedup_se with:
            input:
                "results/te_repeats/mapped/dedup/{sample}.bam",
            output:
                "results/te_repeats/mapped/{sample}.bam",
            log:
                "logs/te_repeats/samtools_sort/dedup/{sample}.log",

    rule te_locus_assignment:
        input:
            bam="results/te_repeats/mapped/{sample}.bam",
            bai="results/te_repeats/mapped/{sample}.bam.bai",
            loci=resources.repeatmasker_loci,
            family_lookup=lambda wildcards: expand(
                "results/te_repeats/candidate_family/{unit}.tsv.gz",
                unit=sample_units(wildcards.sample),
            ),
        output:
            "results/te_repeats/{sample}.locus_reads.tsv.gz",
        log:
            "logs/te_repeats/locus_assignment/{sample}.log",
        conda:
            "../envs/te_repeats.yaml"
        threads: 1
        resources:
            runtime=config["resources"]["peaks"]["time"],
            mem_mb=8000,
        params:
            require_family_match=TE_REPEATS["require_family_match"],
        script:
            "../scripts/te_locus_assignment.py"

    rule te_enrichment:
        input:
            unpack(
                lambda wildcards: {
                    "ip_reads": f"results/te_repeats/{wildcards.sample}.locus_reads.tsv.gz",
                    "input_reads": f"results/te_repeats/{control(wildcards.sample)}.locus_reads.tsv.gz",
                }
            ),
        output:
            locus="results/te_repeats/{sample}.locus_enrichment.tsv",
            family="results/te_repeats/{sample}.family_enrichment.tsv",
        log:
            "logs/te_repeats/enrichment/{sample}.log",
        conda:
            "../envs/te_repeats.yaml"
        threads: 1
        resources:
            runtime=config["resources"]["peaks"]["time"],
        params:
            min_reads=TE_REPEATS["min_reads"],
        script:
            "../scripts/te_enrichment.py"

    rule te_reproducibility:
        input:
            locus1="results/te_repeats/{s1}.locus_enrichment.tsv",
            family1="results/te_repeats/{s1}.family_enrichment.tsv",
            locus2="results/te_repeats/{s2}.locus_enrichment.tsv",
            family2="results/te_repeats/{s2}.family_enrichment.tsv",
        output:
            "results/te_repeats/{s1}_vs_{s2}.reproducibility.txt",
        log:
            "logs/te_repeats/reproducibility/{s1}_vs_{s2}.log",
        conda:
            "../envs/te_repeats.yaml"
        threads: 1
        resources:
            runtime=15,
        script:
            "../scripts/te_reproducibility.py"
