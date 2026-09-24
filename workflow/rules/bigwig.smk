# ENCODE: RPM normalised, strand specific read density files (Yeo lab makebigwigfiles.py,
# reimplemented with bedtools genomecov and UCSC tools)
rule bigwig:
    input:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
        cs=resources.chrom_sizes,
    output:
        pos="results/bigwig/{sample}.norm.pos.bw",
        neg="results/bigwig/{sample}.norm.neg.bw",
    params:
        # Single-end reads are not reversed (direction f), paired-end read 2 is (direction r)
        pos_strand="-" if PAIRED_END else "+",
        neg_strand="+" if PAIRED_END else "-",
    threads: config["resources"]["bigwig"]["cpu"]
    resources:
        runtime=config["resources"]["bigwig"]["time"],
    log:
        "logs/bigwig/{sample}.log",
    conda:
        "../envs/bigwig.yaml"
    shell:
        """
        N=$(samtools view -c -F 4 {input.bam})
        SCALE=$(gawk -v n=$N 'BEGIN {{printf "%.12f", 1000000 / n}}')
        NEG_SCALE=$(gawk -v n=$N 'BEGIN {{printf "%.12f", -1000000 / n}}')
        TMP=$(mktemp -d)
        trap 'rm -rf $TMP' EXIT

        genomeCoverageBed -ibam {input.bam} -bg -strand {params.pos_strand} -scale $SCALE -g {input.cs} -du -split 2> {log} | LC_ALL=C sort -k1,1 -k2,2n > $TMP/pos.bg
        bedGraphToBigWig $TMP/pos.bg {input.cs} {output.pos} >> {log} 2>&1

        genomeCoverageBed -ibam {input.bam} -bg -strand {params.neg_strand} -scale $NEG_SCALE -g {input.cs} -du -split 2>> {log} | LC_ALL=C sort -k1,1 -k2,2n > $TMP/neg.bg
        bedGraphToBigWig $TMP/neg.bg {input.cs} {output.neg} >> {log} 2>&1
        """
