rule get_read2:
    input:
        "results/mapped/{sample}/{sample}_sorted.dedup.bam",
    output:
        "results/mapped/{sample}/{sample}_R2.bam",
    log:
        "logs/samtools/read2_{sample}.log",
    params:
        extra="-f 128", # second in pair flag
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    wrapper:
        "v3.3.3/bio/samtools/view"


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
        "v3.3.3/bio/samtools/index"


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


rule compute_common_CL_motifs:
    input:
        bam="results/mapped/{input_sample}/{input_sample}_sorted.dedup.bam", #pre processed bam file (after deduplication and sorting)
        fasta=resources.fasta,
        txt="resources/motifs.txt",
        xml="resources/motifs.xml",
    output:
        "results/pureclip/common_cl_motifs/fimo_clmotif_occurences_{input_sample}.bed",
    conda:
        "../envs/pureclip.yaml"
    log:
        "logs/pureclip/common_cl_motifs_{input_sample}.log"
    script:
        "../scripts/common_CL_motifs.sh"


rule crosslink_detection:
    input:
        bam="results/mapped/{ip_sample}/{ip_sample}_R2.bam",
        bai="results/mapped/{ip_sample}/{ip_sample}_R2.bam.bai",
        ibam="results/mapped/{input_sample}/{input_sample}_R2.bam",
        ibai="results/mapped/{input_sample}/{input_sample}_R2.bam.bai",
        fasta=resources.fasta,
        common_cl="results/pureclip/common_cl_motifs/fimo_clmotif_occurences_{input_sample}.bed",
    output:
        crosslink_sites="results/pureclip/crosslink_sites/{ip_sample}_vs_{input_sample}.bed",
        par="results/pureclip/parameters/{ip_sample}_vs_{input_sample}.txt",
    threads: config["resources"]["pureclip"]["cpu"]
    resources: 
        runtime=config["resources"]["pureclip"]["time"]
    conda:
        "../envs/pureclip.yaml"
    log:
        "logs/pureclip/pureclip_{ip_sample}_vs_{input_sample}.log"
    shell:
        "pureclip "
        "--nt {threads} "
        "--nta {threads} "
        "-iv '1;2;3;' " # chromosomes for model training (Ensembl format)
        "-i {input.bam} "
        "-bai {input.bai} "
        "-g {input.fasta} "
        "-o {output.crosslink_sites} "
        "--par {output.par} "
        "-nim 4 "
        "-fis {input.common_cl} "
        "-ibam {input.ibam} "
        "-ibai {input.ibai} "
        "> {log} 2>&1"

