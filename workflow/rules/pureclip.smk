rule get_read2:
    input:
        "results/mapped/{sample}/{sample}_sorted.dedup.bam",
    output:
        "results/mapped/{sample}/{sample}_R2.bam",
    log:
        "logs/samtools/read2_{sample}.log",
    params:
        extra="-f 130",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    wrapper:
        "v3.1.0/bio/samtools/view"


rule index_read2_bam:
    input:
        "results/mapped/{sample}/{sample}_R2.bam",
    output:
        "results/mapped/{sample}/{sample}_R2.bam.bai",
    log:
        "logs/samtools/index_R2_{sample}.log",
    params:
        extra="",
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    wrapper:
        "v3.1.0/bio/samtools/index"


rule get_common_CL_motifs_files:
    output:
        txt="resources/motifs.txt",
        xml="resources/motifs.xml",
    conda:
        "../envs/pureclip.yaml"
    log:
        "logs/pureclip/get_common_CL_motifs_files.log"
    shell:
        "wget -O {output.txt} https://raw.githubusercontent.com/skrakau/PureCLIP_data/master/common_CL-motifs/dreme.w10.k4.txt > {log} 2>&1; "
        "wget -O {output.xml} https://raw.githubusercontent.com/skrakau/PureCLIP_data/master/common_CL-motifs/dreme.w10.k4.xml "
        ">> {log} 2>&1" # append to same log as other file


rule common_CL_motifs:
    input:
        bam="results/mapped/{sample}/{sample}_R2.bam",
        fasta=resources.fasta,
        txt="resources/motifs.txt",
        xml="resources/motifs.xml",
    output:
        "results/pureclip/common_cl_motifs/fimo_clmotif_occurences_{sample}.bed",
    conda:
        "../envs/pureclip.yaml"
    log:
        "logs/pureclip/common_cl_motifs_{sample}.log"
    script:
        "../scripts/common_CL_motifs.sh"


rule pureclip:
    input:
        bam="results/mapped/{sample}/{sample}_R2.bam",
        bai="results/mapped/{sample}/{sample}_R2.bam.bai",
        fasta=resources.fasta,
        common_cl="results/pureclip/common_cl_motifs/fimo_clmotif_occurences_{sample}.bed",
    output:
        crosslink_sites="results/pureclip/crosslink_sites_{sample}.bed",
        binding_regions="results/pureclip/binding_regions_{sample}.bed",
        par="results/pureclip/par_{sample}.txt",
    threads: config["resources"]["pureclip"]["cpu"]
    resources: 
        runtime=config["resources"]["pureclip"]["time"]
    conda:
        "../envs/pureclip.yaml"
    log:
        "logs/pureclip/pureclip_{sample}.log"
    shell:
        "pureclip "
        "--nt {threads} "
        "--nta {threads} "
        "-iv '1;2;3;' " # chromosomes for model training (Ensembl format)
        "-i {input.bam} "
        "-bai {input.bai} "
        "-g {input.fasta} "
        "-o {output.crosslink_sites} "
        "--or {output.binding_regions} "
        "--par {output.par} "
        "-nim 4 "
        "-fis {input.common_cl} "
        "> {log} 2>&1"


rule annotate_regions: # includes GO analysis on nearest genes
    input:
        bed="results/pureclip/binding_regions_{sample}.bed",
        homer="resources/homer_genome_installed",
    output:
        bed="results/homer/binding_regions_{sample}_homer.bed",
        txt="results/homer/annotated_regions_{sample}.txt",
        go=directory("results/homer/GO_{sample}"),
    conda:
        "../envs/pureclip.yaml"
    log:
        "logs/homer/annotate_regions_{sample}.log"
    script:
        "../scripts/annotate_regions.sh"

