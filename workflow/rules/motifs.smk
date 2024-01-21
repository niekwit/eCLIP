rule convert_bed2fasta:
    input:
        bed="results/pureclip/crosslink_sites/{ip_sample}_vs_{input_sample}.bed",
        fasta=resources.fasta
    output:
        out="results/motifs/fasta/{ip_sample}_vs_{input_sample}.fa"
    params:
        flank=config["motifs"]["crosslink_flank"]
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/pureclip.yaml"
    log:
        "logs/motifs/convert_bed2fasta/{ip_sample}_vs_{input_sample}.log"
    script:
        "../scripts/convert_bed2fasta.py"


rule create_background_fasta:
    input:
        bed="results/pureclip/crosslink_sites/{ip_sample}_vs_{input_sample}.bed",
        fasta=resources.fasta
    output:
        out="results/motifs/fasta/background_{ip_sample}_vs_{input_sample}.fa"
    params:
        flank=config["motifs"]["crosslink_flank"]
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/pureclip.yaml"
    log:
        "logs/motifs/create_background_fasta/{ip_sample}_vs_{input_sample}.log"
    script:
        "../scripts/create_background_fasta.py"


rule detect_motifs:
    input:
        fasta="results/motifs/fasta/{ip_sample}_vs_{input_sample}.fa",
        background="results/motifs/fasta/background_{ip_sample}_vs_{input_sample}.fa",
        homer="resources/homer_genome_installed",
    output:
        html="results/motifs/{ip_sample}_vs_{input_sample}/homerResults.html",
    params:
        dr=directory("results/motifs/{ip_sample}_vs_{input_sample}/"),
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources: 
        runtime=config["resources"]["deeptools"]["time"]
    conda:
        "../envs/pureclip.yaml"
    log:
        "logs/motifs/detect_motifs/{ip_sample}_vs_{input_sample}.log"
    script:
        "../scripts/detect_motifs.sh"

